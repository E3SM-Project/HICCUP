#!/bin/bash
#SBATCH -C knl
#SBATCH --account=m3312
###SBATCH -p debug
###SBATCH --time=0:30:00
#SBATCH --partition regular
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --mail-user=hannah6@llnl.gov
#SBATCH --mail-type=END,FAIL
###SBATCH --output=logs_slurm/slurm-%x-%j.out
###SBATCH --job-name=hiccup_$(NE)

# To run this batch script, use the command below 
# to set the output grid from the command line:
# NE ~ number of spectral elements on a cube edge 
# (NE=30 roughly corresponds 1 degree grid)
# NE=1024 ; sbatch --job-name=hiccup_ne$NE --output=logs_slurm/slurm-%x-%j.out --export=NE=$NE ./run_hiccup_batch.rhea.sh
# NE=120 ; sbatch --job-name=hiccup_ne$NE --output=logs_slurm/slurm-%x-%j.out --export=NE=$NE ./run_hiccup_batch.rhea.sh
# sbatch --job-name=hiccup_UP --output=logs_slurm/slurm-%x-%j.out ./run_hiccup_batch.rhea.sh

# module load python

# source activate hiccup_env

# Set NE if not set on the command line
# if [ -z ${NE+x} ]; then NE=30; fi

# VGRID=L72
# time python -u ./create_initial_condition.py --hgrid=ne${NE}np4  --vgrid=$VGRID

VGRID=L50
time python -u ./create_initial_condition.py --init_date=2008-10-01 --vgrid=$VGRID
VGRID=L120
time python -u ./create_initial_condition.py --init_date=2008-10-01 --vgrid=$VGRID

# VGRID=L100
# time python -u ./create_initial_condition.py --init_date=2008-10-01 --vgrid=$VGRID
# time python -u ./create_initial_condition.py --init_date=2009-10-01 --vgrid=$VGRID
# time python -u ./create_initial_condition.py --init_date=2010-10-01 --vgrid=$VGRID
# time python -u ./create_initial_condition.py --init_date=2011-10-01 --vgrid=$VGRID
# time python -u ./create_initial_condition.py --init_date=2012-10-01 --vgrid=$VGRID
# time python -u ./create_initial_condition.py --init_date=2013-10-01 --vgrid=$VGRID
# time python -u ./create_initial_condition.py --init_date=2014-10-01 --vgrid=$VGRID
# time python -u ./create_initial_condition.py --init_date=2015-10-01 --vgrid=$VGRID
# time python -u ./create_initial_condition.py --init_date=2016-10-01 --vgrid=$VGRID
# time python -u ./create_initial_condition.py --init_date=2017-10-01 --vgrid=$VGRID
# time python -u ./create_initial_condition.py --init_date=2018-10-01 --vgrid=$VGRID

# Notes:
# "time" is added to add the execution time to the log file.
# the "-u" option allows the log file to update in real time 
# which is useful for montioring the batch job.