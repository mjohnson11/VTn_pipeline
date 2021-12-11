#!/bin/bash
#SBATCH -J VTn_fitness  #job name for array
#SBATCH -n 1                    # Number of cores
#SBATCH -N 1                    # Ensure that all cores are on one machine
#SBATCH -t 0-03:00              # Runtime in D-HH:MM
#SBATCH -p serial_requeue       # Partition to submit to
#SBATCH --mem=18000               # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o ../../shell_outputs/VTn_fitness.out      # File to which STDOUT will be written
#SBATCH -e ../../shell_outputs/VTn_fitness.err      # File to which STDERR will be written
#SBATCH --mail-type=ALL              # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=milo.s.johnson.13@gmail.com  # Email to which notifications will be sent

source activate milo_py37
python Fitness_analysis.py
