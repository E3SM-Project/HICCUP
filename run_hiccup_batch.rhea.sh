#!/bin/bash
#SBATCH --account=CLI115
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --partition=batch
###SBATCH --partition=gpu
###SBATCH --output=logs_slurm/slurm-%x-%j.out
###SBATCH --mail-user=hannah6@llnl.gov
###SBATCH --mail-type=END,FAIL
#SBATCH --mem=0

# To run this batch script, use the command below to set the output grid from the command line:
# NE=1024 ; sbatch --job-name=hiccup_ne$NE --output=logs_slurm/slurm-%x-%j.out --export=NE=$NE ./run_hiccup_batch.rhea.sh
# NE=256 ; sbatch --job-name=hiccup_ne$NE --output=logs_slurm/slurm-%x-%j.out --export=NE=$NE ./run_hiccup_batch.rhea.sh

# Load the python environment
source activate hiccup_env

# Set NE if not set on the command line
if [ -z ${NE+x} ]; then NE=30; fi

# default to L72, but switch to L128 for ne1024 -- not sure why this doesn't work...
# VGRID=L72
# if [$NE -eq 1024]; then VGRID=L128; fi

VGRID=L128

python -u ./create_initial_condition_multifile.olcf.py --hgrid=ne${NE}np4 --vgrid=$VGRID

# Notes:
# the "-u" option allows the log file to update in real time, 
# which is useful for montioring the batch job.