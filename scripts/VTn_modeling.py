import pandas as pd
import numpy as np
from glob import glob
import statsmodels.formula.api as smf

vtn_s = pd.read_csv('../../output/VTn_s.csv')
vtn_x = pd.read_csv('../../output/VTn_x.csv')
vtn_x['Cond'] = vtn_x['Pop'].str[:2] + '_' + vtn_x['Environment']
g70_fits = {i[0]: i[1] for i in np.array(vtn_x[vtn_x.Gen==70][['Cond', 'Fitness']].groupby('Cond').mean().reset_index()[['Cond', 'Fitness']])}
vtn_x['Fitness_sub_70'] = vtn_x.apply(lambda row: row['Fitness']-g70_fits[row['Cond']], axis=1)
vtn_s = vtn_s.merge(vtn_x[['Sample', 'Fitness', 'Freq_T0', 's_VLTE', 'Fitness_sub_70']], on='Sample', how='left')

def indicator(row, gen, pop):
    if row['Pop'] == pop and row['Gen']>=gen:
        return 1
    else:
        return 0
    
p1_pops = ['P1B04','P1C04','P1C05','P1D03','P1F05','P1G04']
p3_pops = ['P3C03','P3D03','P3F03','P3G02','P3G05','P3G06']

# making a bunch of indicator variables that say: this is at or after timepoint X in population Y
gens = [70, 1410, 2640, 5150, 7530, 10150]
for pop in p1_pops+p3_pops:
    for gen in gens[:-1]: # no 10K indicators, that is fitting one point
        vtn_s['ind_'+pop+'_'+str(gen)] = vtn_s.apply(lambda r: indicator(r, gen, pop), axis=1)
        
indicators = [i for i in vtn_s if i[:3]=='ind']
all_edges = sorted(set(vtn_s['Edge']))
conditions = ['P1_YPD_30C', 'P3_SC_37C', 'P1_SC_37C', 'P3_bad_SC_37C']

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
        td = vtn_s[(vtn_s.Edge==edge) & (vtn_s.Cond==cond) & (vtn_s.num_cbcs>=5) & (vtn_s['s_sub_g70_s'].notnull())]
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
                mat.append([edge, cond, model_fixer[base], 1-full_model.ssr/full_model.centered_tss, full_model.bic, ';'.join(coeff_list), 
                            ';'.join([str(coeffs[c]) for c in coeff_list]), ';'.join([str(pvals[c]) for c in coeff_list])])
                if base == 'x':
                    x_model = er[1][0]
                    coeffs = dict(x_model.params)
                    pvals = dict(x_model.pvalues)
                    coeff_list = [c for c in coeffs]
                    mat.append([edge, cond, 'XM', 1-x_model.ssr/x_model.centered_tss, x_model.bic, ';'.join(coeff_list), 
                                ';'.join([str(coeffs[c]) for c in coeff_list]), ';'.join([str(pvals[c]) for c in coeff_list])])
                    
modeling = pd.DataFrame(mat, columns=['Edge', 'Cond', 'Model', 'R2', 'BIC', 'Params', 'Coeffs', 'Pvalues'])
modeling.to_csv('../../output/VTn_modeling.csv', index=False)