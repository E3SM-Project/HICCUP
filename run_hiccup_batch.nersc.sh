#!/bin/bash
#SBATCH -C knl
#SBATCH -p regular
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --mail-user=hannah6@llnl.gov
#SBATCH --mail-type=END,FAIL
#SBATCH --output=slurm-%x-%j.out
###SBATCH --job-name=hiccup_$(NE)

# To run this batch script, use the command below 
# to set the output grid from the command line:
# NE ~ number of spectral elements on a cube edge 
# (NE=30 roughly corresponds 1 degree grid)
# NE=30 ; sbatch --job-name=hiccup_$NE --export=ALL,NE=$NE ./run_hiccup_batch.nersc.sh

module load python

source activate hiccup_env

# Set NE if not set on the command line
if [ -z ${NE+x} ]; then NE=30; fi

time python -u ./create_initial_condition.py --hgrid=ne${NE}np4