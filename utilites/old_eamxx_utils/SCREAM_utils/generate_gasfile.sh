#!/bin/bash
#SBATCH -N 1
###SBATCH -C cpu
#SBATCH -C knl
#SBATCH -t 01:00:00
#SBATCH -q regular
#SBATCH -A e3sm
#SBATCH -o remap_gases-%j.out

# use this command to submit:
# sbatch generate_gasfile.sh

module load python

# conda activate nco-cmem
conda activate hiccup_env

dst_res=256

# output_root=/global/cfs/cdirs/e3sm/bhillma/scream/data/init
output_root=/global/cfs/cdirs/e3sm/whannah/init_data

# input_gas_file=${output_root}/gas_constituents_ne30np4L128_20220711.permuted.nc
input_gas_file=/global/cfs/cdirs/e3sm/bhillma/scream/data/init/gas_constituents_ne30np4L128_20220711.permuted.nc
output_gas_file=${output_root}/gas_constituents_ne${dst_res}np4L128_`date +%Y%m%d`.nc

# map_file=$SCRATCH/grids/ne1024/map_ne30_to_1024_highorder_20220912.nc
map_file=/global/cfs/projectdirs/m3312/whannah/HICCUP/files_map/map_ne30np4_to_ne256np4.20221212.nc
#map_file=${output_root}/ne${dst_res}/map_ne${src_res}_to_${dst_res}_highorder.nc

# Convert data format
srun ncremap -m ${map_file} ${input_gas_file} ${output_gas_file}.permuted
srun ncpdq -a time,ncol,lev ${output_gas_file}.permuted ${output_gas_file}