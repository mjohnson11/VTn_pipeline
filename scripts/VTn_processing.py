import pandas as pd
import numpy as np
from glob import glob

# Reading our data
sdat_V2 = {i.split('/')[-2]: pd.read_csv(i) for i in glob('/n/holyscratch01/desai_lab/mjohnson/V2Tn/output/s_estimation_V2/*/*_edge_s.csv') if 'rep' not in i}
sdat_V1 = {i.split('/')[-2]: pd.read_csv(i) for i in glob('/n/holyscratch01/desai_lab/mjohnson/VTn/output/s_estimation_V2/*/*_edge_s.csv') if 'rep' not in i}
# Formatting our data
tidy_one = []
for samp in sdat_V2:
    td = sdat_V2[samp]
    for row in np.array(td[['Edge', 'mean.s', 'stderr.s', 'total.cbcs']]):
        tidy_one.append([samp+'_SC_37C', 'SC_37C', samp, samp.split('_')[1], int(samp.split('_')[0][1:])]+list(row))
        
for samp in sdat_V1:
    #if 'P3' not in samp:
    td = sdat_V1[samp]
    for row in np.array(td[['Edge', 'mean.s', 'stderr.s', 'total.cbcs']]):
        if samp.split('_')[1][:2] == 'P3':
            # in the first assay, a bad batch of SC media reagents was used. We still process that data but it is not included in the paper
            tidy_one.append([samp+'_bad_SC_37C', 'bad_SC_37C', samp, samp.split('_')[1], int(samp.split('_')[0][1:])]+list(row))
        else:
            tidy_one.append([samp+'_YPD_30C', 'YPD_30C', samp, samp.split('_')[1], int(samp.split('_')[0][1:])]+list(row))
        
tidy_df = pd.DataFrame(tidy_one, columns=['Sample', 'Environment', 'Gen_pop', 'Pop', 'Gen', 'Edge', 's', 'stderr', 'num_cbcs'])
tidy_df['Cond'] = tidy_df['Pop'].str[:2] + '_' + tidy_df['Environment']

# Getting the average s at gen 70 for each mutation in each condition
all_edges = sorted(set(tidy_df['Edge']))
conditions = ['P1_YPD_30C', 'P3_SC_37C', 'P1_SC_37C', 'P3_bad_SC_37C']
anc_s = {c: dict() for c in conditions}
for cond in conditions:
    for edge in all_edges:
        td = tidy_df[(tidy_df.Edge==edge) & (tidy_df.Cond==cond) & (tidy_df.num_cbcs>=5) & (tidy_df.Gen==70)]
        anc_s[cond][edge] = np.nanmean(td['s'])

tidy_df['g70_s'] = tidy_df.apply(lambda row: anc_s[row.Cond][row.Edge], axis=1)
tidy_df['s_sub_g70_s'] = tidy_df['s'] - tidy_df['g70_s'] 

# Formatting our old data (Johnson et al 2019)
edge_info = pd.read_csv('../accessory_files/TP_data_by_edge.csv')
s = edge_info.melt(id_vars='Edge', value_vars=[i for i in edge_info if 'mean.s' in i], value_name='s')
std = edge_info.melt(id_vars='Edge', value_vars=[i for i in edge_info if 'stderr.s' in i and 'rep' not in i], value_name='stderr')
cbcs = edge_info.melt(id_vars='Edge', value_vars=[i for i in edge_info if 'total.cbcs' in i], value_name='num_cbcs')
dfs = [s, std, cbcs]
for df in dfs:
    df['Sample'] = df['variable'].str.split('.').str[0]

cols = ['Edge', 'Sample']
tidy_old_data = s[cols+['s']].merge(std[cols+['stderr']], on=['Edge', 'Sample'], how='inner').merge(cbcs[cols+['num_cbcs']], on=['Edge', 'Sample'], how='inner')

seg_fits = pd.read_csv('../accessory_files/Clones_For_Tn96_Experiment.csv')
simple_seg_fits = seg_fits[seg_fits['segregant'].isin(set(tidy_old_data['Sample']))][['segregant', 'initial fitness, YPD 30C']].rename(columns={'segregant': 'Sample', 'initial fitness, YPD 30C': 'Fitness'})

# Outputting tidy fitness effect (s) data
tidy_df.to_csv('../../output/VTn_s.csv', index=False)
tidy_old_data.to_csv('../../output/BYxRM_s.csv', index=False)
simple_seg_fits.to_csv('../../output/BYxRM_x.csv', index=False)



