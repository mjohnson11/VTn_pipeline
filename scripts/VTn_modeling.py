import pandas as pd
import numpy as np
from glob import glob
import statsmodels.formula.api as smf
import scipy.stats as sci_stats
import scipy.odr as sci_o
from statsmodels.stats.multitest import fdrcorrection as benjamini_hochberg

p1_pops = ['P1B04','P1C04','P1C05','P1D03','P1F05','P1G04']
p3_pops = ['P3C03','P3D03','P3F03','P3G02','P3G05','P3G06']
byrm_s = pd.read_csv('../../output/BYxRM_s.csv')
byrm_x = pd.read_csv('../../output/BYxRM_x.csv')
byrm_s = byrm_s.merge(byrm_x[['Sample', 'Fitness']], on='Sample', how='left')

def indicator(row, gen, pop):
    if row['Pop'] == pop and row['Gen']>=gen:
        return 1
    else:
        return 0
    
    
def expand_df(df):
    df['Gen'] = df['Sample'].apply(lambda s: int(s.split('_')[0][1:]))
    df['Pop'] = df['Sample'].apply(lambda s: s.split('_')[1].split('-')[0])
    df['Env'] = df['Sample'].apply(lambda s: s.split('-')[1])
    df['Cond'] = df.apply(lambda r: r['Pop'][:2]+'_'+r['Env'], axis=1)
    return df


def f(B, x): # linear function for odr
    return B[0]+B[1]*x


def call_slope(row, cond):
    if pd.isnull(row[cond+'_slope']): return np.nan
    if row[cond+'_p']<0.05:
        if row[cond+'_slope'] > 0.05:
            return '+'
        elif row[cond+'_slope'] < -0.05:
            return '-'
    return 'NS'

    
def do_modeling(in_s, in_x, outfile, outfile2, cbc_cutoff=5):

    vtn_s = expand_df(pd.read_csv(in_s))
    vtn_x = expand_df(pd.read_csv(in_x))
    g70_fits = {i[0]: i[1] for i in np.array(vtn_x[vtn_x.Gen==70][['Cond', 'Fitness']].groupby('Cond').mean().reset_index()[['Cond', 'Fitness']])}
    vtn_x['Fitness_sub_70'] = vtn_x.apply(lambda row: row['Fitness']-g70_fits[row['Cond']], axis=1)
    vtn_s = vtn_s.merge(vtn_x[['Sample', 'Fitness', 'Fitness_sub_70', 'Fitness_std']], on='Sample', how='left')
    
    # making a bunch of indicator variables that say: this is at or after timepoint X in population Y
    gens = [70, 1410, 2640, 5150, 7530, 10150]
    for pop in p1_pops+p3_pops:
        for gen in gens[:-1]: # no 10K indicators, that is fitting one point
            vtn_s['ind_'+pop+'_'+str(gen)] = vtn_s.apply(lambda r: indicator(r, gen, pop), axis=1)

    indicators = [i for i in vtn_s if i[:3]=='ind']
    all_edges = sorted(set(vtn_s['Edge']))
    conditions = ['P1_YPD_30C', 'P3_SC_37C', 'P1_SC_37C', 'P3_bad_SC_37C']
    
    # Getting ancestor 
    anc_s = {c: dict() for c in conditions}
    for cond in conditions:
        for edge in all_edges:
            td = vtn_s[(vtn_s.Edge==edge) & (vtn_s.Cond==cond) & (vtn_s.num_cbcs>=cbc_cutoff) & (vtn_s.Gen==70)]
            anc_s[cond][edge] = np.nanmean(td['s'])

    vtn_s['g70_s'] = vtn_s.apply(lambda row: anc_s[row.Cond][row.Edge], axis=1)
    vtn_s['s_sub_g70_s'] = vtn_s['s'] - vtn_s['g70_s'] 

    # Changes:
    # don't allow 10K indicators or any that only fit one point (e.g. a 7.5k indicator, but there is no 10K measurement for that pop)
    # don't allow more than one indicator param per population
    edge_results = {c: {'no_x': dict(), 'x': dict()} for c in conditions}
    s_var = 's_sub_g70_s'
    vtn_s['dummy'] = [0]*len(vtn_s)
    for cond in conditions:
        print(cond)
        cc = 0
        for edge in all_edges:
            cc += 1
            if cc % 10 == 0:
                print(cc)
            td = vtn_s[(vtn_s.Edge==edge) & (vtn_s.Cond==cond) & (vtn_s.num_cbcs>=cbc_cutoff) & (vtn_s['s_sub_g70_s'].notnull())]
            if len(td)>=20:
                for base_case in [['no_x', []], ['x', ['Fitness_sub_70']]]:
                    params = base_case[1]
                    if base_case[0] == 'x':
                        if cond == 'P3_bad_SC_37C':
                            break
                        results = [smf.ols(formula=s_var+' ~ Fitness_sub_70 -1', data=td).fit()]
                    else:
                        results = [smf.ols(formula=s_var+' ~ dummy -1', data=td).fit()]
                    ind_use = indicators
                    while True:
                        rec = {ind: smf.ols(formula=s_var+' ~ ' + ' + '.join([ind] + params) + '-1', data=td).fit() for ind in ind_use}
                        best = sorted(rec.keys(), key=lambda x: rec[x].bic)
                        if rec[best[0]].bic - results[-1].bic >= -2:
                            break
                        params.append(best[0])
                        results.append(rec[best[0]])
                        pops_w_params = set([ind.split('_')[1] for ind in params if ind!='Fitness_sub_70'])
                        # this implements the criteria in the comment at the top
                        ind_use = [ind for ind in ind_use if ind.split('_')[1] not in pops_w_params and len(td[td[ind]==1])>1] 
                        if len(ind_use) == 0:
                            break
                    edge_results[cond][base_case[0]][edge] = [params, results]

    model_fixer = {'x': 'FM', 'no_x': 'IM'}
    mat = []
    for cond in conditions:
        for base in ['x', 'no_x']:
            td = edge_results[cond][base]
            for edge in all_edges:
                if edge in td:
                    er = td[edge]
                    full_model = er[1][-1]
                    coeffs = dict(full_model.params)
                    pvals = dict(full_model.pvalues)
                    coeff_list = [c for c in coeffs]
                    # Using 1-full_model.ssr/full_model.centered_tss to get R2 because otherwise the fixed-intercept model is comparing our predictions to the sum of squares 
                    # of differences from the fixed intercept (mean gen 70 s) rather than from the mean (inflating R2)
                    mat.append([edge, cond, model_fixer[base], 1-full_model.ssr/full_model.centered_tss, full_model.llf, full_model.bic, ';'.join(coeff_list), 
                                ';'.join([str(coeffs[c]) for c in coeff_list]), ';'.join([str(pvals[c]) for c in coeff_list])])
                    if base == 'x':
                        x_model = er[1][0]
                        coeffs = dict(x_model.params)
                        pvals = dict(x_model.pvalues)
                        coeff_list = [c for c in coeffs]
                        mat.append([edge, cond, 'XM', 1-x_model.ssr/x_model.centered_tss, x_model.llf, x_model.bic, ';'.join(coeff_list), 
                                    ';'.join([str(coeffs[c]) for c in coeff_list]), ';'.join([str(pvals[c]) for c in coeff_list])])

    modeling = pd.DataFrame(mat, columns=['Edge', 'Cond', 'Model', 'R2', 'LLF', 'BIC', 'Params', 'Coeffs', 'Pvalues'])
    # see note above about using centered_tss for R2 calculation - this can mean an R2 below zero, which we will change to nan
    modeling['R2'] = np.clip(modeling['R2'], 0, 1)
    modeling.to_csv(outfile, index=False)
    
    # Reformatting modeling data
    modeling['Cmodel'] = modeling['Cond']+'_'+modeling['Model']
    vcols = ['R2', 'LLF', 'BIC', 'Params', 'Coeffs', 'Pvalues']
    dats = [modeling.pivot(index='Edge', columns='Cmodel', values=v).reset_index() for v in vcols]
    edge_models = dats[0]
    base_cols = [i for i in edge_models if i!='Edge']
    for i in range(1, len(dats)):
        edge_models = edge_models.merge(dats[i], on='Edge', how='outer', suffixes=('', '_'+vcols[i]))
    edge_models = edge_models.rename(columns={i: i+'_R2' for i in base_cols})
    
    # Gene annotation data etc.
    edge_info = pd.read_csv('../accessory_files/TP_data_by_edge.csv')
    e2g = {i[0]:i[1] for i in np.array(edge_info[['Edge', 'Gene.Use']])}

    mat = []
    for edge in set(vtn_s['Edge']):
        tmp = [edge]
        dfs = [byrm_s] + [vtn_s[vtn_s.Cond==cond] for cond in conditions]
        things = ['BYxRM'] + conditions
        i = 0
        for df in dfs:
            # For each condition
            # Filtering for >= 5 cBCs (or >= 3 cBCs for clone modeling)
            td = df[(pd.notnull(df['s'])) & (df['Edge']==edge) & (df['num_cbcs']>=cbc_cutoff)]
            if len(td) >= 20: # recording mean s and variance of s
                tmp += [np.mean(td['s']), np.var(td['s'])]
            else:
                tmp += [np.nan, np.nan]
            td = df[(pd.notnull(df['Fitness'])) & (pd.notnull(df['s'])) & (df['Edge']==edge) & (df['num_cbcs']>=cbc_cutoff)]
            if len(td) >= 20: # recording regression results
                lr = sci_stats.linregress(td['Fitness'], td['s'])
                tmp += [lr[0], lr[3], lr[2]**2] # slope, P, R^2
            else:
                #print(things[i], e2g[edge], edge)
                tmp += [np.nan, np.nan, np.nan]
            i += 1
        mat.append(tmp)
    cols = ['Edge']
    for c in ['BYxRM'] + conditions:
        cols += [c+'_s_mean', c+'_s_var', c+'_slope', c+'_p', c+'_x_R2']
    # turning it into a dataframe
    edge_stats = pd.DataFrame(mat, columns=cols)
    for cond in ['BYxRM'] + conditions:
        edge_stats[cond+'_call'] = edge_stats.apply(lambda row: call_slope(row, cond), axis=1)
        td = edge_stats[edge_stats[cond+'_p'].notnull()]
        edges = list(td['Edge'])
        corrected_ps = benjamini_hochberg(list(td[cond+'_p']))[1]
        p_dict = {edges[i]:corrected_ps[i] for i in range(len(edges))}                                                                                          
        edge_stats[cond+'_bh_p'] = edge_stats['Edge'].apply(lambda e: p_dict.get(e, np.nan))
        

    edge_short = edge_info[['Edge', 'chromosome', 'Type', 'Gene.Use', 'briefDescription', 'insertion_edge', 'phenotypeSummary', 'phenotypeSummary.nearby']].rename(columns={'Gene.Use': 'Gene_Use'})
    edge_stats = edge_stats.merge(edge_short, on='Edge', how='left') # adding Gene annotations etc.
    edge_stats = edge_stats.merge(edge_models, on='Edge', how='outer') # adding modeling data
    edge_stats.to_csv(outfile2, index=False)

    

do_modeling('../../output/VTn_s.csv', '../../output/VTn_x.csv', '../../output/VTn_modeling.csv', '../../output/VTn_modeling_wide.csv')
do_modeling('../../output/VTn_s_SHUFFLED.csv', '../../output/VTn_x.csv', '../../output/VTn_modeling_SHUFFLED.csv', '../../output/VTn_modeling_wide_SHUFFLED.csv')
do_modeling('../../output/VTn_s_simulated.csv', '../../output/VTn_x.csv', '../../output/VTn_modeling_simulated.csv', '../../output/VTn_modeling_wide_simulated.csv')
do_modeling('../../output/VTn_s_clones.csv', '../../output/VTn_x_clones.csv', '../../output/VTn_modeling_clones.csv', '../../output/VTn_modeling_wide_clones.csv', cbc_cutoff=3)