#!/bin/bash
#SBATCH --constraint=cpu
#SBATCH --account=m3312
#SBATCH --qos=regular
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --mail-user=hannah6@llnl.gov
#SBATCH --mail-type=END,FAIL

# To run this batch script, use the command below 
# to set the output grid from the command line:
# (NOte that NE ~ number of spectral elements on a cube edge, and NE=30 roughly corresponds 1 degree grid)
# NE=120 ; sbatch --job-name=hiccup_ne$NE --output=logs_slurm/slurm-%x-%j.out --export=NE=$NE ./run_hiccup_batch.rhea.sh


# sbatch template_scripts/run_hiccup_batch.nersc.sh


# Load the python environment
source activate hiccup_env


# VGRID=L120
# time python -u ./create_initial_condition.py --init_date=2008-10-01 --vgrid=$VGRID

# sbatch template_scripts/run_hiccup_batch.nersc.sh
time python -u custom_scripts/2024.SCREAM_autocalibration.process_nudging_data.py

# Notes:
# "time" is used here to print the execution time to the end of log file.
# the "-u" option allows the log file to update in real time 
# which is useful for montioring the batch job.