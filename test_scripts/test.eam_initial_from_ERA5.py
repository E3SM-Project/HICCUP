#!/usr/bin/env python3
# ==============================================================================
# HICCUP - Hindcast Initial Condition Creation Utility/Processor
# This script tests the HICCUP functionality for creating an 
# atmospheric initial condition file from ERA5 reanalysis
# ==============================================================================
import os
from hiccup import hiccup
hiccup.hdc.print_memory_usage = True
# ------------------------------------------------------------------------------
# local path for grid and mapping files (move to scratch space for large grids)
hiccup_root  = os.getenv('HOME')+'/HICCUP'
data_root    = f'{hiccup_root}/test_data'
data_tmp     = f'{hiccup_root}/test_data_tmp'       # local machine
# din_loc_root = os.getenv('HOME')+'/E3SM/inputdata'  # local machine
# ------------------------------------------------------------------------------
# NERSC paths
din_loc_root = '/global/cfs/cdirs/e3sm/inputdata'
data_tmp     = os.getenv('SCRATCH')+'/HICCUP/test_data_tmp'
# ------------------------------------------------------------------------------
# LLNL paths
# din_loc_root = '/p/lustre2/hannah6/inputdata'
# data_tmp     = '/p/lustre1/hannah6/hiccup_scratch/test_data_tmp'
# ------------------------------------------------------------------------------

os.makedirs(data_tmp, exist_ok=True)  # create temporary output data path if it doesn't exist

dst_horz_grid = 'ne30np4' # ne30np4 / ne120np4 / ne512np4 / ne1024np4

dst_vert_grid ='L80'; vert_file_name = f'{hiccup_root}/files_vert/L80_for_E3SMv3.nc' # E3SMv3

# Specify output file names
output_atm_file_name = f'{data_tmp}/HICCUP_TEST_OUTPUT.eam_from_era5.{dst_horz_grid}.{dst_vert_grid}.nc'

if dst_horz_grid=='ne30np4'  : topo_file = f'{data_root}/USGS-gtopo30_ne30np4_16xdel2-PFC-consistentSGH.nc'
if dst_horz_grid=='ne120np4' : topo_file = f'{din_loc_root}/atm/cam/topo/USGS-gtopo30_ne120np4pg2_16xdel2.nc'
if dst_horz_grid=='ne512np4' : topo_file = f'{din_loc_root}/atm/cam/topo/USGS-gtopo30_ne512np4pg2_x6t_20230404.nc'
if dst_horz_grid=='ne1024np4': topo_file = f'{din_loc_root}/atm/cam/topo/USGS-gtopo30_ne1024np4pg2_x6t-SGH.c20210614.nc'

# Create data class instance, which includes xarray file dataset objects
# and variable name dictionaries for mapping between naming conventions.
# This also checks input files for required variables
hiccup_data = hiccup.create_hiccup_data( src_data_name='ERA5',
                                         target_model='EAM',
                                         dst_horz_grid=dst_horz_grid,
                                         dst_vert_grid=dst_vert_grid,
                                         atm_file=f'{data_root}/HICCUP_TEST.ERA5.atm.low-res.nc',
                                         sfc_file=f'{data_root}/HICCUP_TEST.ERA5.sfc.low-res.nc',
                                         topo_file=topo_file,
                                         grid_dir=data_tmp,
                                         map_dir=data_tmp,
                                         tmp_dir=data_tmp,
                                         verbose=True,
                                         check_input_files=True,)

# Print some informative stuff
print('\n  Input Files')
print(f'    input atm files: {hiccup_data.atm_file}')
print(f'    input sfc files: {hiccup_data.sfc_file}')
print(f'    input topo file: {hiccup_data.topo_file}')
print('\n  Output files')
print(f'    output atm file: {output_atm_file_name}')
print()

# ------------------------------------------------------------------------------
# reduced variable set for quicker testing
atm_keys_all = hiccup_data.atm_var_name_dict.copy().keys()
sfc_keys_all = hiccup_data.sfc_var_name_dict.copy().keys()
atm_keys_keep = ['lat','lon','T','Q','U','O3']
sfc_keys_keep = ['PS','TS','PHIS']
for key in atm_keys_all:
  if key not in atm_keys_keep: del hiccup_data.atm_var_name_dict[key]
for key in sfc_keys_all:
  if key not in sfc_keys_keep: del hiccup_data.sfc_var_name_dict[key]
# ------------------------------------------------------------------------------
# Get dict of temporary files for each variable
file_dict = hiccup_data.get_multifile_dict(timestamp=999)
# ------------------------------------------------------------------------------
# Make sure files are "unpacked" (may take awhile, so only do it if you need to)
# hiccup_data.unpack_data_files()
# ------------------------------------------------------------------------------
# Create grid description files needed for the mapping file
hiccup_data.create_src_grid_file()
hiccup_data.create_dst_grid_file()
# Create mapping file
hiccup_data.create_map_file()
# ------------------------------------------------------------------------------
# Horizontally regrid the data
hiccup_data.remap_horizontal_multifile(file_dict=file_dict)
# Rename variables to match what the model expects
hiccup_data.rename_vars_multifile(file_dict=file_dict)
# Add time/date information
hiccup_data.add_time_date_variables_multifile(file_dict=file_dict)
# Do surface adjustments
hiccup_data.surface_adjustment_multifile(file_dict=file_dict)
# Vertically remap the data
hiccup_data.remap_vertical_multifile(file_dict=file_dict,
                                     vert_file_name=vert_file_name)
# Perform final state adjustments on interpolated data and add additional data
hiccup_data.atmos_state_adjustment_multifile(file_dict=file_dict)
# Combine and delete temporary files
hiccup_data.combine_files(file_dict=file_dict,delete_files=True,
                          output_file_name=output_atm_file_name)
# Clean up the global attributes of the file
hiccup_data.clean_global_attributes(file_name=output_atm_file_name)
# ------------------------------------------------------------------------------
# Print final output file name
print(); print(f'output_atm_file_name: {output_atm_file_name}'); print()
# Print summary of timer info
hiccup_data.print_timer_summary()
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
