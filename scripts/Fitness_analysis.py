import numpy as np
from scipy import stats as sci_stats
import pandas as pd
# this was installed like: pip install FlowCytometryTools
import FlowCytometryTools as fct
from glob import glob
from collections import defaultdict

####  SETTING UP PLATE LAYOUTS ETC  ####
def well_to_plate_well(well, plates):
    # From 384 well plate to plate name and 96 well plate well name
    row = ord(well[0])-65
    col = int(well[1:])
    plate = plates[row % 2][(col-1) %2]
    plate_row = chr((row-2)//2+66)
    plate_col = (col-1)//2+1
    return plate + '_' + plate_row + str(plate_col).zfill(2)


lets = [chr(65+i) for i in range(8)]
lets_long = [chr(65+i) for i in range(16)]
gens = [70, 1410, 2640, 5150, 7530, 10150]

P1_wells = ['P1B04', 'P1G04', 'P1C04', 'P1F05', 'P1D03', 'P1C05']
P3_wells = ['P3F03', 'P3G02', 'P3G06', 'P3G05', 'P3C03', 'P3D03']
# making plate layout dictionaries for where samples are in the 96 well plates
P1_plate = dict()
P3_plate = dict()
for col, gen in [(1, 1410), (11, 5150), (12, 10150)]:
    for i in range(6):
        P1_plate[lets[i]+str(col).zfill(2)] = 'G'+str(gen)+'_'+P1_wells[i]
        P3_plate[lets[i]+str(col).zfill(2)] = 'G'+str(gen)+'_'+P3_wells[i]

all_P1_clones = ['G'+str(g)+'_'+p+'_'+c for c in ['A', 'B'] for g in gens[::-1] for p in P1_wells]
all_P3_clones = ['G'+str(g)+'_'+p+'_'+c for c in ['A', 'B'] for g in gens[::-1] for p in P3_wells]
for i in range(len(all_P1_clones)):
    P1_plate[lets[i%8]+str(2+i//8).zfill(2)] = all_P1_clones[i]
    P3_plate[lets[i%8]+str(2+i//8).zfill(2)] = all_P3_clones[i]
    


# layouts for each assay (the ignore ones are a test I did where I preconditioned in YPD and assayed in SC, which I am not using)
pl_fa12 = [['VTnC_P1_1', 'VTnC_P1_2'], ['VTnC_P3_1', 'VTnC_P3_2']]
plate_layout = {'FA1': pl_fa12, 'FA2': pl_fa12, 'FA3': [['VTnC_P1_1', 'VTnC_P1_2'], ['ignore', 'ignore']]} 

    

#### GATES ####
p1_fortessa = fct.PolyGate([(6800, -1), (6500, 5200), (7000, 6000), (9000, 7700), (10000, 11000), (10000, -1)], ['FITC-A', 'PE-A'])
sc_fortessa_ssc = fct.PolyGate([(5000, 5200), (11000, 11500), (13000, 12000), (10000, -1)], ['FITC-A', 'SSC-A'])
cg = fct.ThresholdGate(7000, 'SSC-A', 'above')
sc_fortessa_ssc_fa3 = fct.PolyGate([(5300, 5200), (10000, 9500), (13000, 12000), (10000, -1)], ['FITC-A', 'SSC-A'])
cg_fa3 = fct.ThresholdGate(5750, 'SSC-A', 'above')


#### READING DATA ####
def get_plate_data(dir_base, plate_set, cellgate, plate_type='384'):
    # given a directory base like 'FA_timepoint_0/Specimen_001_', this will try to read in all the files corresponding to each well in a 96 or 384 well plate
    # and will return a dictionary like td['well_id'] = dataframe with facs info
    print('reading from', dir_base)
    td = dict()
    wells_missed = []
    if plate_type == '384':
        let_top, col_top = 16, 25
    elif plate_type == '96':
        let_top, col_top = 8, 13
    else:
        print('unrecognized plate type')
        return None
    for let in [chr(i+65) for i in range(let_top)]:
        for col in range(1, col_top):
            orig_well = let + str(col).zfill(2)
            if plate_type == '384':
                well = well_to_plate_well(orig_well, plate_set)
            else:
                well = orig_well
            flist = glob(dir_base + '*' + orig_well + '.fcs')
            try:
                assert len(flist) == 1
                # Reading in file and immediately gating on good cells
                samp = fct.FCMeasurement(ID=well, datafile=flist[0]).transform('tlog', channels=['FITC-A', 'SSC-A', 'PE-A'])
                if cellgate:
                    samp = samp.gate(cellgate)
                td[well] = samp
            except AssertionError:
                wells_missed.append(well)
                td[well] = None
    if len(wells_missed) > 0:
        print('Missed files for', len(wells_missed))
    return td
    

cellgates = {'FA1': None, 'FA2': cg, 'FA3': cg_fa3}
tps = [0,1,2,3]
dirs = ['../../data/FACS/VTn_FA'+f+'_T'+str(tp) for tp in tps for f in ['1', '2', '3']]
dir_d = {d.split('/')[-1]: d for d in dirs}
dat_d = dict()
for d in dir_d:
    dat_d[d] = get_plate_data(dir_d[d] + '/Specimen_001_', plate_layout[d.split('_')[-2]], cellgates[d.split('_')[-2]])
    
    

#### MEASURING REF COUNTS AND GETTING FITNESSES ####
def get_ref_counts(df, use_gate):
    ref_counts = df.gate(use_gate).shape[0]
    total_counts = df.shape[0]
    density = df.shape[0]/np.nanmax(df['Time'])
    if total_counts < 1000: # Excluding timepoints with less than 1000 reads
        freq = np.nan
        #print('Low counts for a sample...')
    else:
        freq = ref_counts/total_counts
    return ref_counts, total_counts-ref_counts, freq, density

## Getting blank refs
blanks_1 = ['G01', 'G11', 'G12', 'H11'] #excluding corners (H01 and H12)
blanks_2 = ['E01', 'F01', 'G01']
assays = [('P1_YPD', 'FA1', 'VTnC_P1', p1_fortessa, blanks_1, [0,1,2,3]), # name, assay code, prefix, gate, blank wells, tps for blanks
          ('P3_SC', 'FA2', 'VTnC_P3', sc_fortessa_ssc, blanks_1, [0,1,2,3]), 
          ('P1_SC', 'FA3', 'VTnC_P1', sc_fortessa_ssc_fa3, blanks_2, [1,2,3])] # just forgot to measure the blanks at T0

blank_ref_rec = defaultdict(list)
# Using P1 ref-only wells to find the not glowing % in YPD and P3 ref-only wells to do the same in SC 37C
for info in assays:
    for t in info[5]:
        for w in info[4]:
            for r in ['1', '2']:
                if dat_d['VTn_' + info[1] + '_T'+str(t)][info[2]+'_'+r+'_'+w]:
                    blank_ref_rec[info[0]].append(get_ref_counts(dat_d['VTn_' + info[1] + '_T'+str(t)][info[2]+'_'+r+'_'+w], info[3])[2])

def get_fit(row, tps):
    # Getting ref frequencies, excluding low count (nan) timepoints
    ref_freqs = np.array([row['Ref_Freq_T' + str(t)] for t in tps if pd.notnull(row['Ref_Freq_T' + str(t)])])
    times = np.array([t for t in tps if pd.notnull(row['Ref_Freq_T' + str(t)])])*10
    # excluding time intervals where ref or test is >95% in both timepoints
    use_tp_until = len(times)
    for t in range(1, len(times)):
        if (ref_freqs[t] > 0.95 and ref_freqs[t-1] > 0.95) or (ref_freqs[t] < 0.05 and ref_freqs[t-1] > 0.05):
            use_tp_until = t-1
            break
    if use_tp_until > 0:
        test_freqs = 1-ref_freqs
        # s = log slope of test freq / reference freq
        return sci_stats.linregress(times[:use_tp_until+1], np.log(test_freqs[:use_tp_until+1]/ref_freqs[:use_tp_until+1]))[0] 
    
final_d = dict()
tps_use = [0,1,2,3]
for name, assay_code, prefix, gate, jnk1, jnk2 in assays:
    for r in ['1', '2']:
        mat = []
        for row in range(8):
            for col in range(12):
                well = chr(row+65) + str(col+1).zfill(2)
                if dat_d['VTn_' + assay_code + '_T0'][prefix+'_1_'+well]: # excludes the blanks in FA3 without throwing errors
                    tmp = [name, r, well]
                    for tp in [0, 1, 2, 3]:
                        result = get_ref_counts(dat_d['VTn_' + assay_code + '_T' + str(tp)][prefix+'_'+r+'_'+well], gate) 
                        tmp += list(result) + [(1/np.nanmean(blank_ref_rec[name]))*result[2]]
                    mat.append(tmp)
        tmp_colnames = ['Assay', 'Rep', 'Well'] + [i+t for t in ['0', '1', '2', '3'] for i in ['Ref_Counts_T', 'NonRef_Counts_T', 'Uncorrected_Ref_Freq_T', 'Density_T', 'Ref_Freq_T']]
        td = pd.DataFrame(mat, columns=tmp_colnames)
        td['s'] = td.apply(lambda r: get_fit(r, tps_use), axis=1)
        final_d[name+'_R'+r] = td
        
# The P1 in P3 data uses a different reference than the P3 in P3 data, so we will shift all fitnesses measured there by the fitness difference between the references (about 0.04)
P3_old_raw_data = pd.read_csv('../accessory_files/P3_freq_and_s_data.csv')
# The 2490A clone is in well D04
ref_clone_fits = np.array(P3_old_raw_data[P3_old_raw_data['Well']=='D04'][[i for i in P3_old_raw_data if '_s_R' in i]])[0]
ref_dif = np.median(ref_clone_fits)
for a in ['P1_SC_R1', 'P1_SC_R2']:
    td = final_d[a]
    td['s_raw'] = td['s']
    td['s'] = td['s_raw']+ref_dif
        
#### OUTPUTTING ####
all_data = pd.concat([final_d[d] for d in final_d])
all_data.to_csv('../../output/FACS_assay/all_FACS_counts_and_s.csv', index=False)

comb_d = dict()
assay_names = ['P1_YPD', 'P3_SC', 'P1_SC']
a_plates = {'P1': P1_plate, 'P3': P3_plate}
cols = ['Well', 's', 'Ref_Freq_T0', 'Ref_Freq_T1']
for a in assay_names:
    r1 = all_data[(all_data['Assay']==a) & (all_data['Rep']==1)]
    r2 = all_data[(all_data['Assay']==a) & (all_data['Rep']==2)]
    comb_d[a] = r1[cols].merge(r2[cols], on='Well', how='inner', suffixes=('_R1', '_R2'))
    comb_d[a]['Strain'] = comb_d[a]['Well'].map(a_plates[a[:2]])
    comb_d[a] = comb_d[a][pd.notnull(comb_d[a]['Strain'])]
    comb_d[a]['s'] = np.mean(comb_d[a][['s_R1', 's_R2']], axis=1)
    comb_d[a]['stderr'] = np.std(comb_d[a][['s_R1', 's_R2']], axis=1)/np.sqrt(2)
    comb_d[a]['Ref_Freq_T0'] = np.nanmean(comb_d[a][['Ref_Freq_T0_R1', 'Ref_Freq_T0_R2']], axis=1)
    comb_d[a]['Ref_Freq_T1'] = np.nanmean(comb_d[a][['Ref_Freq_T1_R1', 'Ref_Freq_T1_R2']], axis=1)
    
old_rearrange = dict()
for p in ['P1', 'P3']:
    mat = []
    for jnk, row in pd.read_csv('../accessory_files/'+p+'_freq_and_s_data.csv').iterrows():
        for gen in gens:
            mat.append(['G'+str(gen)+'_'+p+row['Well'], row['Gen'+str(gen)+'_s'], row['Gen'+str(gen)+'_s_scaled']])
        
    old_rearrange[p] = pd.DataFrame(mat, columns=['Gen_pop', 's_VLTE', 's_VLTE_scaled'])
    
compare_results = dict()
for a in comb_d:
    tmd = comb_d[a]
    pops = sorted(set(['_'.join(i.split('_')[:2]) for i in list(tmd['Strain'])]))
    mat = []
    for p in pops:
        tmp = [p]
        td = tmd[tmd['Strain']==p]
        if len(td)==1:
            tmp.append(td.iloc[0]['s'])
            tmp.append(td.iloc[0]['Ref_Freq_T0'])
        else:
            tmp += [np.nan, np.nan]
        td = tmd[tmd['Strain']==p+'_A']
        assert len(td)==1
        tmp.append(td.iloc[0]['s'])
        tmp.append(td.iloc[0]['stderr'])
        tmp.append(td.iloc[0]['Ref_Freq_T0'])
        td = tmd[tmd['Strain']==p+'_B']
        assert len(td)==1
        tmp.append(td.iloc[0]['s'])
        tmp.append(td.iloc[0]['stderr'])
        tmp.append(td.iloc[0]['Ref_Freq_T0'])
        mat.append(tmp)
    
    compare_results[a] = pd.DataFrame(mat, columns=['Gen_pop', 'Pop_s', 'Pop_Ref_Freq_T0', 'Clone_A_s', 'Clone_A_s_stderr', 'Clone_A_Ref_Freq_T0', 'Clone_B_s', 'Clone_B_s_stderr', 'Clone_B_Ref_Freq_T0'])
    compare_results[a]['Clone_s'] = np.nanmean(compare_results[a][['Clone_A_s', 'Clone_B_s']], axis=1)


for a in ['P1_YPD', 'P3_SC']:
    compare_results[a] = compare_results[a].merge(old_rearrange[a.split('_')[0]], on='Gen_pop', how='left')
    compare_results[a]['Gen'] = compare_results[a]['Gen_pop'].str.split('_').str[0].str[1:]
    
env_fix = {'YPD': 'YPD_30C', 'SC': 'SC_37C'}
for a in compare_results:
    compare_results[a]['Sample'] = compare_results[a]['Gen_pop']+'-'+[env_fix[a.split('_')[1]]]*len(compare_results[a])
    compare_results[a]['Fitness'] = np.nanmean(compare_results[a][['Clone_A_s', 'Clone_B_s']], axis=1)
    compare_results[a]['Fitness_std'] = np.nanstd(compare_results[a][['Clone_A_s', 'Clone_B_s']], axis=1, ddof=1)/np.sqrt(2)
    compare_results[a]['Freq_T0'] = 1-np.nanmean(compare_results[a][['Clone_A_Ref_Freq_T0', 'Clone_B_Ref_Freq_T0']], axis=1)
    
    
together = pd.concat([compare_results[a] for a in compare_results])
cols = ['Sample', 'Fitness', 'Fitness_std', 'Freq_T0', 's_VLTE', 's_VLTE_scaled',
        'Clone_A_s', 'Clone_A_s_stderr', 'Clone_A_Ref_Freq_T0', 
        'Clone_B_s', 'Clone_B_s_stderr', 'Clone_B_Ref_Freq_T0', 
        'Pop_Ref_Freq_T0', 'Pop_s']

vtn_x = together[cols]
vtn_x.to_csv('../../output/VTn_x.csv', index=False)

s_cols = {'s': 'Fitness', 's_stderr': 'Fitness_std'}
r1 = vtn_x[['Sample']+['Clone_A_'+s for s in s_cols]].rename(columns={'Clone_A_'+s:s_cols[s] for s in s_cols})
r1['Clone'] = 'A'
r2 = vtn_x[['Sample']+['Clone_B_'+s for s in s_cols]].rename(columns={'Clone_B_'+s:s_cols[s] for s in s_cols})
r2['Clone'] = 'B'
vtn_s_clones = pd.concat([r1, r2])
vtn_s_clones['Sample'] = vtn_s_clones['Sample']+'-'+vtn_s_clones['Clone']
vtn_s_clones[['Sample', 'Fitness', 'Fitness_std']].to_csv('../../output/VTn_x_clones.csv', index=False)