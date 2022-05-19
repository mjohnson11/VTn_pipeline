import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as pl
import matplotlib.lines as lines
import seaborn as sns
from scipy import stats as sci_stats
from statsmodels.stats.multitest import fdrcorrection as benjamini_hochberg
import matplotlib.gridspec as gridspec
colors = sns.color_palette('colorblind')
pl.rcParams['font.family'] = 'sans-serif'
pl.rcParams['font.sans-serif'] = 'Noto Sans'

def expand_df(df):
    df['Gen'] = df['Sample'].apply(lambda s: int(s.split('_')[0][1:]))
    df['Pop'] = df['Sample'].apply(lambda s: s.split('_')[1].split('-')[0])
    df['Env'] = df['Sample'].apply(lambda s: s.split('-')[1])
    df['Cond'] = df.apply(lambda r: r['Pop'][:2]+'_'+r['Env'], axis=1)
    return df

conditions = ['P1_YPD_30C', 'P3_SC_37C', 'P1_SC_37C']
p1_pops = ['P1B04','P1C04','P1C05','P1D03','P1F05','P1G04']
p3_pops = ['P3C03','P3D03','P3F03','P3G02','P3G05','P3G06']
cond_pops = {'P1': p1_pops, 'P3': p3_pops}
color_map = {p1_pops[i]:colors[i] for i in range(len(p1_pops))}
color_map.update({p3_pops[i]:colors[i] for i in range(len(p3_pops))})
cond_to_title = {'P1_YPD_30C': 'YPD 30°C', 'P3_SC_37C': 'SC 37°C', 'P1_SC_37C': 'Evo. YPD 30°C, in SC 37°C',
                 'P3_bad_SC_37C': 'Evo. SC 37°C, in bad media'}

# This dictionary changes the recorded generation numbers to the correct generation numbers
# Since P3 only does 8 gens/day it is very different, the other differences are due to little recording errors
gen_fixer = {70: {'P1': 70, 'P2': 70, 'P3': 56, 'P4': 70},
             550: {'P1': 560, 'P2': 560, 'P3': 448, 'P4': 560},
             1410: {'P1': 1410, 'P2': 1410, 'P3': 1128, 'P4': 1410},
             2640: {'P1': 2640, 'P2': 2640, 'P3': 2104, 'P4': 2640},
             3630: {'P1': 3660, 'P2': 3660, 'P3': 2920, 'P4': 3660},
             5150: {'P1': 5170, 'P2': 5170, 'P3': 4128, 'P4': 5170},
             7530: {'P1': 7550, 'P2': 7560, 'P3': 6040, 'P4': 7560},
             10150: {'P1': 10190, 'P2': 10200, 'P3': 8096, 'P4': 10200}}

vtn_s = expand_df(pd.read_csv('../../output/VTn_s.csv'))
vtn_s_clones = expand_df(pd.read_csv('../../output/VTn_s_clones.csv'))
vtn_x = expand_df(pd.read_csv('../../output/VTn_x.csv'))
vtn_x_clones = expand_df(pd.read_csv('../../output/VTn_x_clones.csv'))
vtn_s = vtn_s.merge(vtn_x[['Sample', 'Fitness', 'Fitness_std']], on='Sample', how='left')
vtn_s_clones = vtn_s_clones.merge(vtn_x_clones[['Sample', 'Fitness', 'Fitness_std']], on='Sample', how='left')
vtn_s = vtn_s[vtn_s.Cond.isin(conditions)]
vtn_s_clones = vtn_s_clones[vtn_s_clones.Cond.isin(conditions)]


byrm_s = pd.read_csv('../../output/BYxRM_s.csv')
byrm_x = pd.read_csv('../../output/BYxRM_x.csv')
byrm_s = byrm_s.merge(byrm_x, on='Sample', how='left')

vtn_dfe = expand_df(pd.read_csv('../../output/vtn_dfe_stats.csv'))
vtn_clones_dfe = expand_df(pd.read_csv('../../output/vtn_clones_dfe_stats.csv'))
vtn_dfe = vtn_dfe.merge(vtn_x[['Sample', 'Fitness', 'Fitness_std']], on='Sample', how='left')
vtn_clones_dfe = vtn_clones_dfe.merge(vtn_x_clones[['Sample', 'Fitness', 'Fitness_std']], on='Sample', how='left')

byrm_dfe = pd.read_csv('../../output/byrm_dfe_stats.csv')
byrm_dfe = byrm_dfe.merge(byrm_x[['Sample', 'Fitness', 'Fitness_std']], on='Sample', how='left')


g70_fits = {i[0]: i[1] for i in np.array(vtn_x[vtn_x.Gen==70][['Cond', 'Fitness']].groupby('Cond').mean().reset_index()[['Cond', 'Fitness']])}
# Getting the average s at gen 70 for each mutation in each condition
all_edges = sorted(set(vtn_s['Edge']))
anc_s = {c: dict() for c in conditions}
for cond in conditions:
    for edge in all_edges:
        td = vtn_s[(vtn_s.Edge==edge) & (vtn_s.Cond==cond) & (vtn_s.num_cbcs>=5) & (vtn_s.Gen==70)]
        anc_s[cond][edge] = np.nanmean(td['s'])

vtn_s['g70_s'] = vtn_s.apply(lambda row: anc_s[row.Cond][row.Edge], axis=1)
vtn_s['s_sub_g70_s'] = vtn_s['s'] - vtn_s['g70_s']

g70_fits_clones = {i[0]: i[1] for i in np.array(vtn_x_clones[vtn_x_clones.Gen==70][['Cond', 'Fitness']].groupby('Cond').mean().reset_index()[['Cond', 'Fitness']])}
# clones now Getting the average s at gen 70 for each mutation in each condition
anc_s_clones = {c: dict() for c in conditions}
for cond in conditions:
    for edge in all_edges:
        td = vtn_s_clones[(vtn_s_clones.Edge==edge) & (vtn_s_clones.Cond==cond) & (vtn_s_clones.num_cbcs>=3) & (vtn_s_clones.Gen==70)]
        anc_s_clones[cond][edge] = np.nanmean(td['s'])

vtn_s_clones['g70_s'] = vtn_s_clones.apply(lambda row: anc_s_clones[row.Cond][row.Edge], axis=1)
vtn_s_clones['s_sub_g70_s'] = vtn_s_clones['s'] - vtn_s_clones['g70_s']

vtn_modeling = pd.read_csv('../../output/VTn_modeling_wide.csv')
vtn_modeling_clones = pd.read_csv('../../output/VTn_modeling_wide_clones.csv')
vtn_modeling_shuf = pd.read_csv('../../output/VTn_modeling_wide_SHUFFLED.csv')
vtn_modeling_sim = pd.read_csv('../../output/VTn_modeling_wide_simulated.csv')

e2g = {i[0]:i[1] for i in np.array(vtn_modeling[['Edge', 'Gene_Use']])}