#!/bin/bash
#SBATCH --account=CLI115
#SBATCH --time=20:00:00
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --output=slurm-%x-%j.out
###SBATCH --mail-user=hannah6@llnl.gov
###SBATCH --mail-type=END,FAIL

# To run this batch script, use the command below to set the output grid from the command line:
# NE=1024 ; sbatch --job-name=hiccup_ne$NE --export=NE=$NE ./run_hiccup_batch.rhea.sh

source activate hiccup_env

# Set NE if not set on the command line
if [ -z ${NE+x} ]; then NE=30; fi

time python -u ./create_initial_condition.py --hgrid=ne${NE}np4

# Notes:
# "time" is added to add the execution time to the log file.
# the "-u" option allows the log file to update in real time 
# which is useful for montioring the batch job.