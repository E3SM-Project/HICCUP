#!/usr/bin/env python
# ==================================================================================================
# HICCUP - Hindcast Initial Condition Creation Utility/Processor
# This tool automates the creation of atmospheric initial condition files for 
# E3SM using user supplied file for atmospheric and sea surface conditions.
# ==================================================================================================
import os, optparse
import hiccup_data_class as hdc
import hiccup_state_adjustment as hsa
# ------------------------------------------------------------------------------
# Parse the command line options
parser = optparse.OptionParser()
parser.add_option('--hgrid',dest='horz_grid',default=None,help='Sets the output horizontal grid')
parser.add_option('--vgrid',dest='vert_grid',default=None,help='Sets the output vertical grid')
(opts, args) = parser.parse_args()
# ------------------------------------------------------------------------------
# Logical flags for controlling what this script will do
verbose = True            # Global verbosity flag
# ------------------------------------------------------------------------------

# Specify output atmosphere horizontal grid
dst_horz_grid = 'ne4np4'

# Specify output atmosphere vertical grid
# dst_vert_grid = opts.vert_grid if opts.vert_grid is not None else 'L100'
# vert_file_name = os.getenv('HOME')+f'/HICCUP/files_vert/vert_coord_E3SM_{dst_vert_grid}.nc'
# vert_file_name = os.getenv('HOME')+f'/HICCUP/files_vert/UP_L125.nc'
# dst_vert_grid,vert_file_name = 'L100',os.getenv('HOME')+f'/E3SM/vert_grid_files/L100_v1.nc'
dst_vert_grid,vert_file_name = 'L72',os.getenv('HOME')+f'/E3SM/vert_grid_files/L72_v1.nc'


# Specify the output file names
data_root = os.getenv('SCRATCH')+'/HICCUP/data/' # NERSC
# data_root = os.getenv('MEMBERWORK')+'/cli115/HICCUP/data/'  # OLCF
# init_date = '0001-01-01'
# init_year = int(init_date.split('-')[0])
# output_atm_file_name = f'{data_root}HICCUP.atm_era5.{init_date}.{dst_horz_grid}.{dst_vert_grid}.nc'
# output_sst_file_name = f'{data_root}HICCUP.sst_noaa.{init_date}.nc'

# set topo file - replace this with file path if no defalt is set
topo_file_name = hdc.get_default_topo_file_name(dst_horz_grid)
# topo_file_name = '/project/projectdirs/e3sm/inputdata/atm/cam/topo/USGS_conusx4v1-tensor12x_consistentSGH_c150924.nc'

scratch_ic_path = '/global/cfs/cdirs/e3sm/inputdata/atm/cam/inic/homme'
# if dst_horz_grid=='ne4np4' : cami_file = f'{scratch_ic_path}/cami_aquaplanet_ne4np4_L72_c190218.nc'
# if dst_horz_grid=='ne30np4': cami_file = f'{scratch_ic_path}/cami_aquaplanet_ne30np4_L72_c190215.nc'
if dst_horz_grid=='ne4np4' : cami_file = f'{scratch_ic_path}/cami_mam3_Linoz_ne4np4_L72_c160909.nc'
if dst_horz_grid=='ne30np4': cami_file = f'{scratch_ic_path}/cami_mam3_Linoz_ne30np4_L72_c160214.nc'
if dst_horz_grid=='ne45np4': cami_file = f'{scratch_ic_path}/cami_mam3_Linoz_ne45np4_L72_c20200611.nc'

# output_atm_file_name = f'{data_root}HICCUP.AQUA.{dst_horz_grid}.{dst_vert_grid}.nc'
# output_atm_file_name = f'{data_root}HICCUP.{dst_horz_grid}.{dst_vert_grid}.nc'
output_atm_file_name = f'{data_root}HICCUP.cami_mam3_Linoz_{dst_horz_grid}.{dst_vert_grid}_alt.nc'

# Create data class instance, which includes xarray file dataset objects
# and variable name dictionaries for mapping between naming conventions.
# This also checks input files for required variables
hiccup_data = hdc.create_hiccup_data(name='EAM'
                                    ,atm_file=cami_file
                                    ,sfc_file=cami_file
                                    ,topo_file=topo_file_name
                                    ,dst_horz_grid=dst_horz_grid
                                    ,dst_vert_grid=dst_vert_grid
                                    ,output_dir=data_root
                                    ,grid_dir=data_root
                                    ,map_dir=data_root
                                    ,tmp_dir=data_root
                                    ,verbose=verbose)

# ------------------------------------------------------------------------------
# Vertically remap the data
# ------------------------------------------------------------------------------

hiccup_data.remap_vertical(input_file_name=cami_file,
                           output_file_name=output_atm_file_name,
                           vert_file_name=vert_file_name)

# ------------------------------------------------------------------------------
# Print final output file name
# ------------------------------------------------------------------------------

print()
print(f'output_atm_file_name: {output_atm_file_name}')
# print(f'output_sst_file_name: {output_sst_file_name}')
print()

# Print summary of timer info
hdc.print_timer_summary()

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
