#!/usr/bin/env python
# ==================================================================================================
# HICCUP - Hindcast Initial Condition Creation Utility/Processor
# This tool automates the creation of atmospheric initial condition files for 
# E3SM using user supplied file for atmospheric and sea surface conditions.
# ==================================================================================================
import os
from hiccup import hiccup
# ------------------------------------------------------------------------------
# HICCUP uses two independent path roots:
#   hiccup_root - path to your local HICCUP repo, only used to locate the bundled
#                 vertical grid files in files_vert/. Leave this pointed at the
#                 repo - it does NOT need to move to scratch.
#   data_root   - path to your input data AND everything HICCUP generates (grid,
#                 map, tmp, and output files). Point this at scratch space when
#                 working on an HPC system, especially for high resolution grids.
# ------------------------------------------------------------------------------
hiccup_root = os.getenv('HOME')+'/HICCUP'
data_root   = os.getenv('SCRATCH')+'/HICCUP/data/' # NERSC
# data_root = os.getenv('MEMBERWORK')+'/cli115/HICCUP/data/'  # OLCF

# Path for "supported" input data
inputdata_root = '/global/cfs/cdirs/e3sm/inputdata'

# time stamp for output file (= datetime.datetime.utcnow().strftime('%Y%m%d')])
timestamp = '20220707'

# output horizontal grid for atmosphere
dst_horz_grid = 'ne30np4'

# output vertical grid for atmosphere
dst_vert_grid,vert_file_name = 'L80',f'{hiccup_root}/files_vert/L80_for_E3SMv3.nc'

# specify input file name
cami_file = f'{inputdata_root}/atm/cam/inic/homme/cami_mam3_Linoz_ne30np4_L72_c160214.nc'

# specify output file - make sure the destination folder exists
if not os.path.exists(data_root): os.makedirs(data_root)
output_atm_file_name = f'{data_root}HICCUP.eam_i_mam3_Linoz_{dst_horz_grid}_{dst_vert_grid}_c{timestamp}.nc'

# topo file of output grid - replace this with file path if no default is set
topo_file_name = hiccup.get_default_topo_file_name(dst_horz_grid)

# ------------------------------------------------------------------------------
# Create HICCUP data class instance

# this includes xarray file dataset objects and variable 
# name dictionaries for mapping between naming conventions.
# This also checks input files for required variables
hiccup_data = hiccup.create_hiccup_data(src_data_name='EAM',
                                        target_model='EAM', # options: EAM / EAMXX
                                        dst_horz_grid=dst_horz_grid,
                                        dst_vert_grid=dst_vert_grid,
                                        input_file_list=[cami_file],
                                        topo_file=topo_file_name,
                                        verbose=True,)

# Print some informative stuff
print('\n  Input Files')
print(f'    input file:      {hiccup_data.input_file_list[0]}')
print(f'    input topo file: {hiccup_data.topo_file}')
print('\n  Output files')
print(f'    output atm file: {output_atm_file_name}')

# ------------------------------------------------------------------------------
# Vertically remap the data

hiccup_data.remap_vertical(input_file_name=cami_file,
                           output_file_name=output_atm_file_name,
                           vert_file_name=vert_file_name)

# ------------------------------------------------------------------------------
# Print final output file name

print()
print(f'output_atm_file_name: {output_atm_file_name}')
# print(f'output_sst_file_name: {output_sst_file_name}')
print()

# Print summary of timer info
hiccup_data.print_timer_summary()

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
