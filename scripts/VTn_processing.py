import pandas as pd
import numpy as np
from glob import glob
import scipy.stats as sci_stats
from statsmodels.stats.multitest import fdrcorrection as benjamini_hochberg

conditions = ['P1_YPD_30C', 'P3_SC_37C', 'P1_SC_37C']

# Reading our data
sdat_V2 = {i.split('/')[-2]: pd.read_csv(i) for i in glob('/n/holyscratch01/desai_lab/mjohnson/V2Tn/output/s_estimation_V2/*/*_edge_s.csv') if 'rep' not in i}
sdat_V1 = {i.split('/')[-2]: pd.read_csv(i) for i in glob('/n/holyscratch01/desai_lab/mjohnson/VTn/output/s_estimation_V2/*/*_edge_s.csv') if 'rep' not in i}
# Formatting our data
tidy_one = []
use_cols = ['Edge', 'mean.s', 'std_bcs', 'std_neut', 'stderr.s', 'total.cbcs', 'pval', 
            'rep1.cbcs', 'rep1.s', 'rep1_std_bcs', 'rep1_std_neut', 'rep1.stderr.s', 
            'rep2.cbcs', 'rep2.s', 'rep2_std_bcs', 'rep2_std_neut', 'rep2.stderr.s']
form_cols = ['Edge', 's', 'std_bcs', 'std_neut', 'stderr', 'num_cbcs', 'pval', 
             'rep1_num_cbcs', 'rep1_s', 'rep1_std_bcs', 'rep1_std_neut', 'rep1_stderr', 
             'rep2_num_cbcs', 'rep2_s', 'rep2_std_bcs', 'rep2_std_neut', 'rep2_stderr']
for samp in sdat_V2:
    td = sdat_V2[samp]
    for row in np.array(td[use_cols]):
        #tidy_one.append([samp+'_SC_37C', 'SC_37C', samp, samp.split('_')[1], int(samp.split('_')[0][1:])]+list(row))
        tidy_one.append([samp+'-SC_37C']+list(row))
        
for samp in sdat_V1:
    #if 'P3' not in samp:
    td = sdat_V1[samp]
    for row in np.array(td[use_cols]):
        if samp.split('_')[1][:2] == 'P3':
            # in the first assay, a bad batch of SC media reagents was used. We still process that data but it is not included in the paper
            #tidy_one.append([samp+'_bad_SC_37C', 'bad_SC_37C', samp, samp.split('_')[1], int(samp.split('_')[0][1:])]+list(row))
            tidy_one.append([samp+'-bad_SC_37C']+list(row))
        else:
            #tidy_one.append([samp+'_YPD_30C', 'YPD_30C', samp, samp.split('_')[1], int(samp.split('_')[0][1:])]+list(row))
            tidy_one.append([samp+'-YPD_30C']+list(row))
        
vtn_s = pd.DataFrame(tidy_one, columns=['Sample']+form_cols)
#vtn_s['Cond'] = vtn_s['Pop'].str[:2] + '_' + vtn_s['Environment']

"""
# Getting the average s at gen 70 for each mutation in each condition
all_edges = sorted(set(vtn_s['Edge']))
conditions = ['P1_YPD_30C', 'P3_SC_37C', 'P1_SC_37C', 'P3_bad_SC_37C']
anc_s = {c: dict() for c in conditions}
for cond in conditions:
    for edge in all_edges:
        td = vtn_s[(vtn_s.Edge==edge) & (vtn_s.Cond==cond) & (vtn_s.num_cbcs>=5) & (vtn_s.Gen==70)]
        anc_s[cond][edge] = np.nanmean(td['s'])

vtn_s['g70_s'] = vtn_s.apply(lambda row: anc_s[row.Cond][row.Edge], axis=1)
vtn_s['s_sub_g70_s'] = vtn_s['s'] - vtn_s['g70_s'] 
"""
# Formatting our old data (Johnson et al 2019)
edge_info = pd.read_csv('../accessory_files/TP_data_by_edge.csv')
s = edge_info.melt(id_vars='Edge', value_vars=[i for i in edge_info if 'mean.s' in i], value_name='s')
std = edge_info.melt(id_vars='Edge', value_vars=[i for i in edge_info if 'stderr.s' in i and 'rep' not in i], value_name='stderr')
cbcs = edge_info.melt(id_vars='Edge', value_vars=[i for i in edge_info if 'total.cbcs' in i], value_name='num_cbcs')
pvals = edge_info.melt(id_vars='Edge', value_vars=[i for i in edge_info if 'pval' in i], value_name='pval')
dfs = [s, std, cbcs, pvals]
for df in dfs:
    df['Sample'] = df['variable'].str.split('.').str[0]

cols = ['Edge', 'Sample']
byrm_s = s[cols+['s']].merge(std[cols+['stderr']], on=['Edge', 'Sample'], how='inner').merge(cbcs[cols+['num_cbcs']], on=['Edge', 'Sample'], how='inner').merge(pvals[cols+['pval']], on=['Edge', 'Sample'], how='inner')

seg_fits = pd.read_csv('../accessory_files/Clones_For_Tn96_Experiment.csv')
byrm_x = seg_fits[seg_fits['segregant'].isin(set(byrm_s['Sample']))][['segregant', 'initial fitness, YPD 30C', 'std err']].rename(columns={'segregant': 'Sample', 'initial fitness, YPD 30C': 'Fitness', 'std err': 'Fitness_std'})

# Benjamini-Hochberg correction for p values
byrm_s['sig']=benjamini_hochberg(byrm_s['pval'])[0]
vtn_s['sig']=benjamini_hochberg(vtn_s['pval'])[0]

# Outputting tidy fitness effect (s) data
vtn_s.to_csv('../../output/VTn_s.csv', index=False)
byrm_s = byrm_s[pd.notnull(byrm_s['s'])] # didn't have to do this for vtn_s bc it has no nulls
byrm_s.to_csv('../../output/BYxRM_s.csv', index=False)
byrm_x.to_csv('../../output/BYxRM_x.csv', index=False)

# Randomly shuffling data
vtn_s['Cond'] = vtn_s.apply(lambda r: r['Sample'].split('_')[1][:2]+'_'+r['Sample'].split('-')[1], axis=1)
dfs = []
for i in vtn_s.groupby(['Cond', 'Edge']):
    td = i[1]
    td['Sample'] = list(td['Sample'].sample(frac=1)) # randomly shuffle the samples
    dfs.append(td)
pd.concat(dfs).to_csv('../../output/VTn_s_SHUFFLED.csv')

# simulating data
dfs = []
for cond in conditions:
    df = vtn_s[vtn_s['Cond']==cond][['Sample', 'Edge', 'stderr', 'num_cbcs']]
    for edge in set(list(df.Edge)):
        td = df[(df.Edge==edge) & (df.num_cbcs >= 5)].copy()
        # we pull s values from a normal distribution with standard dev equal to the mean standard error in our empirical data
        std_mean = np.mean(td.stderr)
        td['s'] = np.random.normal(loc=0, scale=std_mean, size=len(td))
        td['num_cbcs'] = 10
        dfs.append(td[['Sample', 'Edge', 's', 'num_cbcs']])
pd.concat(dfs).to_csv('../../output/VTn_s_simulated.csv', index=False)

s_cols = ['s', 'std_bcs', 'std_neut', 'stderr', 'num_cbcs']
r1 = vtn_s[['Sample', 'Edge']+['rep1_'+s for s in s_cols]].rename(columns={'rep1_'+s:s for s in s_cols})
r1['Clone'] = 'A'
r2 = vtn_s[['Sample', 'Edge']+['rep2_'+s for s in s_cols]].rename(columns={'rep2_'+s:s for s in s_cols})
r2['Clone'] = 'B'
vtn_s_clones = pd.concat([r1, r2])
vtn_s_clones['Sample'] = vtn_s_clones['Sample']+'-'+vtn_s_clones['Clone']
vtn_s_clones[['Sample', 'Edge']+s_cols].to_csv('../../output/VTn_s_clones.csv', index=False)

## CALCULATING DFE STATS ##
## BY x RM (old experiment)
# cutoff is 60 mutations for by x rm
byrm_svc = dict(byrm_s.Sample.value_counts())
byrm_s_for_dfe = byrm_s[byrm_s.Sample.apply(lambda s: byrm_svc[s]>=60)]
byrm_dfe = byrm_s_for_dfe[['Sample', 's']].groupby('Sample').mean().reset_index().rename(columns={'s': 'DFE_mean'})
byrm_dfe = byrm_dfe.merge(byrm_s_for_dfe[['Sample', 's']].groupby('Sample').std().reset_index().rename(columns={'s': 'DFE_std'}), on='Sample')
byrm_dfe = byrm_dfe.merge(byrm_s_for_dfe[['Sample', 's']].groupby('Sample').count().reset_index().rename(columns={'s': 's_count'}), on='Sample')
byrm_dfe.to_csv('../../output/byrm_dfe_stats.csv', index=False)

## VTn DFE ##
def get_dfe_stats(df, complete_dataset_cutoff=40):
    # Given a dataframe with columns: ['Sample', 'Edge', 's', 'std_bcs', 'std_neut', 'stderr', 'num_cbcs']
    # Return a dataframe with columns: ['Sample', 's_count', 'del_s_count', 'DFE_mean', 'DFE_std', 'DFE_skew', 'DFE_mean_std']

    # filtering the df
    df['Cond'] = df.apply(lambda r: r['Sample'].split('_')[1][:2]+'_'+r['Sample'].split('-')[1], axis=1)
    df = df[(df['num_cbcs']>=3) & (df['Cond'].isin(conditions))]
    # get the average effect for each edge
    edge_avgs = {c: {i[0]:i[1] for i in np.array(df[df['Cond']==c][['Edge', 's']].groupby('Edge').mean().reset_index()[['Edge', 's']])} for c in conditions}
    edge_avg_stds = {c: {i[0]:i[1] for i in np.array(df[df['Cond']==c][['Edge', 'std_bcs']].groupby('Edge').mean().reset_index()[['Edge', 'std_bcs']])} for c in conditions}

    # getting a list of mutations that are deleterious on average for each condition:
    del_edges = {c: [i for i in edge_avgs[c] if edge_avgs[c][i]<-0.05] for c in edge_avgs}
    df['broadly_deleterious'] = df.apply(lambda r: r['Edge'] in del_edges[r['Cond']], axis=1)
    # Samples must have at least 60 mutations with s measurements in one condition to be included in DFE-level analyses:
    svc = {cond: dict(df[df.Cond==cond].Sample.value_counts()) for cond in conditions}
    df_for_dfe = df[df.apply(lambda r: svc[r['Cond']][r['Sample']]>=60, axis=1)]
    # Going through the samples for each condition in order of the number of mutations they have measurements for until
    # we hit a cutoff of the set of shared_mutations
    shared_dataset = dict() # keys are conditions, values like [[samples], {edges}] where all edges are measured in all samples
    for cond in conditions:
        samples_in_order = sorted([i for i in svc[cond] if svc[cond][i]>=60], key=lambda x: -1*svc[cond][x])
        sample_counter = 0
        new_edges = set(df_for_dfe[df_for_dfe['Sample']==samples_in_order[sample_counter]]['Edge'])
        while len(new_edges) >= complete_dataset_cutoff:
            edges = new_edges
            sample_counter += 1
            if sample_counter == len(samples_in_order):
                break
            new_edges = edges.intersection(set(df_for_dfe[df_for_dfe['Sample']==samples_in_order[sample_counter]]['Edge']))
        shared_dataset[cond] = [samples_in_order[:sample_counter], edges]
        print('Shared for', cond, 'Samples:', len(samples_in_order[:sample_counter]), 'Edges:', len(edges))
    df_for_dfe['shared_dataset'] = df_for_dfe.apply(lambda r: r['Sample'] in shared_dataset[r['Cond']][0] and r['Edge'] in shared_dataset[r['Cond']][1], axis=1)
    mat = []
    for sample in set(df_for_dfe['Sample']):
        td = df_for_dfe[df_for_dfe['Sample']==sample]
        td_shared = td[td['shared_dataset']]
        s_list = list(td['s'])
        s_list_shared = list(td_shared['s'])
        s_list_filled = s_list + [edge_avgs[td['Cond'].iloc[0]][e] for e in edge_avgs[td['Cond'].iloc[0]] if e not in list(td['Edge'])]
        std_list_filled = list(td['std_bcs']) + [edge_avg_stds[td['Cond'].iloc[0]][e] for e in edge_avgs[td['Cond'].iloc[0]] if e not in list(td['Edge'])]
        tmp = [sample, len(td), len(td[td['broadly_deleterious']]), np.mean(s_list), np.mean(s_list_filled), np.nanmean(s_list_shared),
               len(s_list_shared), np.std(s_list), sci_stats.skew(s_list)]
        neut_std = td['std_neut'].iloc[0] # this is the same for all edges for this sample, so taking the first is ok
        dfe_mean_std = np.sqrt(np.sum(np.power(td['std_bcs'], 2))/(len(td)**2) + neut_std**2)
        shared_dfe_mean_std = np.sqrt(np.sum(np.power(td_shared['std_bcs'], 2))/(len(td)**2) + neut_std**2)
        filled_dfe_mean_std = np.sqrt(np.sum(np.power(np.array(std_list_filled), 2))/(len(td)**2) + neut_std**2)
        tmp += [dfe_mean_std, filled_dfe_mean_std, shared_dfe_mean_std]
        mat.append(tmp)
        
    return pd.DataFrame(mat, columns=['Sample', 's_count', 'del_s_count', 'DFE_mean', 'filled_DFE_mean', 'shared_DFE_mean', 'shared_DFE_count',
                                      'DFE_std', 'DFE_skew', 'DFE_mean_std', 'filled_DFE_mean_std', 'shared_DFE_mean_std'])


vtn_dfe = get_dfe_stats(vtn_s)
vtn_dfe.to_csv('../../output/vtn_dfe_stats.csv', index=False)

vtn_clones_dfe = get_dfe_stats(vtn_s_clones)
vtn_clones_dfe.to_csv('../../output/vtn_clones_dfe_stats.csv', index=False)