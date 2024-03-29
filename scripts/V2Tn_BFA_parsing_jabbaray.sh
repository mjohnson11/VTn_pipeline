#!/bin/bash
#SBATCH -J V2Tn_BFA_parse  #job name for array
#SBATCH -n 1                    # Number of cores
#SBATCH -N 1                    # Ensure that all cores are on one machine
#SBATCH -t 0-08:00              # Runtime in D-HH:MM
#SBATCH -p serial_requeue       # Partition to submit to
#SBATCH --mem=20000               # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o ../../shell_outputs/V2Tn_BFA_%a.out      # File to which STDOUT will be written
#SBATCH -e ../../shell_outputs/V2Tn_BFA_%a.err      # File to which STDERR will be written
#SBATCH --mail-type=ALL              # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=milo.s.johnson.13@gmail.com  # Email to which notifications will be sent

source activate milo_py37

python BT_parse.py ../accessory_files/V2Tn_demult.csv /n/holyscratch01/desai_lab/mjohnson/V2Tn/output/ /n/holyscratch01/desai_lab/mjohnson/V2Tn/data/reads/ "${SLURM_ARRAY_TASK_ID}" V2Tn_BFA
