import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
from simple_s_estimation_v2 import s_estimation
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('demult_id_key', help='key to identify which file to parse (which row of the demult file) - for jabba-rays (job arrays)')
parser.add_argument('exp_name', help='experiment name')
args = parser.parse_args()

demult_id_key = int(args.demult_id_key)
exp_name = args.exp_name

ll_cutoff = 40 # log-likelihood ratio cutoff for excluding outliers
bc_target_num = 5 # maximum # of bcs to reduce to (see function in simple_s_estimation)

experiment = 'TP' # this really means which experiment are the plasmid libraries from
input_base = '/n/holyscratch01/desai_lab/mjohnson/' + exp_name + '/output/'
output_base = '/n/holyscratch01/desai_lab/mjohnson/' + exp_name + '/output/s_estimation_V2/'
edge_info_file = '../accessory_files/Tn96_edges_chosen_final.csv'
replicate_info_file = '../accessory_files/'+exp_name+'_replicate_info.csv'

# Reading extra info
edge_info = pd.read_csv(edge_info_file)
# excluding the TDA11 reference since it has effects in some clones / environments:
neut_edges = [i[:15] for i in edge_info.loc[edge_info['new.num.sig']==0]['Edge'] if i!='TATATTGAACTTTAC']
rep_d = pd.read_csv(replicate_info_file)
rep_info = {i[0]: i[1:] for i in rep_d.as_matrix(['Pop_timepoint', 'Clone_A', 'Clone_B'])}
population_timepoints = list(rep_d['Pop_timepoint'])

for i in range(demult_id_key*5, (demult_id_key+1)*5):
    if i == len(population_timepoints):
        break
    plate = population_timepoints[i].split('_')[1][:2]
    if plate == 'P1':
        gens_per_day = 10
    elif plate == 'P3':
        gens_per_day = 8
    if exp_name == 'V2Tn':
        s_estimation([population_timepoints[i]], rep_info, output_base, input_base, experiment, ll_cutoff, bc_target_num, neut_edges, gens_per_day=gens_per_day, exclude_V2Tn_xcontam=True)
    else:
        s_estimation([population_timepoints[i]], rep_info, output_base, input_base, experiment, ll_cutoff, bc_target_num, neut_edges, gens_per_day=gens_per_day)

