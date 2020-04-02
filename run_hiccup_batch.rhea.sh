#!/bin/bash
#SBATCH --account=CLI115
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --partition=batch
###SBATCH --partition=gpu
###SBATCH --output=slurm_logs/slurm-%x-%j.out
###SBATCH --mail-user=hannah6@llnl.gov
###SBATCH --mail-type=END,FAIL
#SBATCH --mem=0

# To run this batch script, use the command below to set the output grid from the command line:
# NE=1024 ; sbatch --job-name=hiccup_ne$NE --output=slurm_logs/slurm-%x-%j.out --export=NE=$NE ./run_hiccup_batch.rhea.sh
# NE=1024 ; sbatch --job-name=hiccup_mf_ne$NE --output=slurm_logs/slurm-%x-%j-mf.out --export=NE=$NE ./run_hiccup_batch.rhea.sh
# NE=120 ; sbatch --job-name=hiccup_mf_ne$NE --output=slurm_logs/slurm-%x-%j-mf.out --export=NE=$NE ./run_hiccup_batch.rhea.sh

# Load the python environment
source activate hiccup_env

# Set NE if not set on the command line
if [ -z ${NE+x} ]; then NE=30; fi

python -u ./create_initial_condition_multifile.olcf.py --hgrid=ne${NE}np4

# Notes:
# the "-u" option allows the log file to update in real time, 
# which is useful for montioring the batch job.