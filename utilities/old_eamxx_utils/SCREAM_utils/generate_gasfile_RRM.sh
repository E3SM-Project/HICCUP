#!/bin/bash
#SBATCH --nodes=1
#SBATCH --constraint=cpu
#SBATCH --account=m2637
#SBATCH --qos=regular
#SBATCH --time=2:00:00
#SBATCH -o remap_gases-%j.out
#SBATCH --mail-user=hannah6@llnl.gov
#SBATCH --mail-type=END,FAIL

# use this command to submit
# sbatch generate_gasfile_RRM.sh

module load python/3.9-anaconda-2021.11

source activate hiccup_env

# output_root=/global/cfs/cdirs/e3sm/bhillma/scream/data/init
output_root=/global/cfs/cdirs/e3sm/whannah/init_data

# input_gas_file=${output_root}/gas_constituents_ne30np4L128_20220711.permuted.nc
input_gas_file=/global/cfs/cdirs/e3sm/bhillma/scream/data/init/gas_constituents_ne30np4L128_20220711.permuted.nc
output_gas_file=${output_root}/gas_constituents_ne0np4-saomai-128x8-pg2-L128_`date +%Y%m%d`.nc

# commands for map file generation
# SRC_GRID=/global/homes/w/whannah/E3SM/data_grid/ne30.g
# DST_GRID=/global/cfs/cdirs/m2637/jsgoodni/Saomai_2006_ne128x8_lon130E_lat25N.g
# MAP_FILE=/global/cfs/cdirs/m2637/whannah/map_ne30np4_to_ne0np4-saomai-128x8.nc
# ncremap -a tempest --src_grd=${SRC_GRID} --dst_grd=${DST_GRID} --map_file=${MAP_FILE} --wgt_opt='--in_type cgll --in_np 4 --out_type cgll --out_np 4 --out_double'

map_file=/global/cfs/cdirs/m2637/whannah/map_ne30np4_to_ne0np4-saomai-128x8.nc

# Convert data format
CMD=ncremap -m ${map_file} ${input_gas_file} ${output_gas_file}.permuted
echo; echo srun $CMD
echo; srun $CMD

CMD=ncpdq --ovr -a time,ncol,lev ${output_gas_file}.permuted ${output_gas_file}
echo; echo srun $CMD
echo; srun $CMD