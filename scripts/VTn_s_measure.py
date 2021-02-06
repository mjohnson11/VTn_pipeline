import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
from simple_s_estimation import s_estimation
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('demult_id_key', help='key to identify which file to parse (which row of the demult file) - for jabba-rays (job arrays)')
args = parser.parse_args()

demult_id_key = int(args.demult_id_key)

ll_cutoff = 15 # log-likelihood ratio cutoff for excluding outliers
bc_target_num = 7 # maximum # of bcs to reduce to (see function in simple_s_estimation)

experiment = 'VTn'
input_base = '/n/holyscratch01/desai_lab/mjohnson/VTn/output/'
output_base = '/n/holyscratch01/desai_lab/mjohnson/VTn/output/s_estimation/'
edge_info_file = '../accessory_files/Tn96_edges_chosen_final.csv'
replicate_info_file = '../accessory_files/VTn_replicate_info.csv'

# Reading extra info
edge_info = pd.read_csv(edge_info_file)
neut_edges = [i[:15] for i in edge_info.loc[edge_info['new.num.sig']==0]['Edge']]
rep_d = pd.read_csv(replicate_info_file)
rep_info = {i[0]: i[1:] for i in rep_d.as_matrix(['Pop_timepoint', 'Clone_A', 'Clone_B'])}
population_timepoints = list(rep_d['Pop_timepoint'])

s_estimation(population_timepoints[demult_id_key*5:(demult_id_key+1)*5], rep_info, output_base, input_base, experiment, ll_cutoff, bc_target_num, neut_edges)
