#!/usr/bin/env python3
# ==================================================================================================
# HICCUP - Hindcast Initial Condition Creation Utility/Processor
# This script tests the HICCUP functionality for creating sea surface 
# temperature and sea ice data files from NOAA data
# ==================================================================================================
import os
from hiccup import hiccup_data_class as hdc
from hiccup import hiccup_state_adjustment as hsa
# ------------------------------------------------------------------------------

# local path for grid and mapping files (move this a scratch space for large grids)
hiccup_root = os.getenv('HOME')+'/HICCUP'
data_root = f'{hiccup_root}/test_data'
data_tmp = f'{hiccup_root}/test_data_tmp'

os.makedirs(data_tmp, exist_ok=True)  # create temporary output data path if it doesn't exist

dst_horz_grid = 'ne30np4'

# Specify output file names
output_sst_file_name = f'{data_tmp}/HICCUP_TEST_OUTPUT.NOAA.sstice.{dst_horz_grid}.nc'

# Create data class instance, which includes xarray file dataset objects
# and variable name dictionaries for mapping between naming conventions.
# This also checks input files for required variables
hiccup_data = hdc.create_hiccup_data(name='NOAA'
                                    ,dst_horz_grid=dst_horz_grid
                                    ,dst_vert_grid='L72'
                                    # ,atm_file=f'{data_root}/HICCUP_TEST.ERA5.atm.low-res.nc'
                                    # ,sfc_file=f'{data_root}/HICCUP_TEST.ERA5.sfc.low-res.nc'
                                    ,sstice_name='NOAA'
                                    ,sst_file=f'{data_root}/HICCUP_TEST.NOAA.sst.nc'
                                    ,ice_file=f'{data_root}/HICCUP_TEST.NOAA.ice.nc'
                                    # ,topo_file=f'{data_root}/USGS-gtopo30_ne30np4_16xdel2-PFC-consistentSGH.nc'
                                    ,grid_dir=data_tmp
                                    ,map_dir=data_tmp
                                    ,tmp_dir=data_tmp
                                    ,verbose=True
                                    ,check_input_files=False)

# Print some informative stuff
print('\n  Input Files')
# print(f'    input atm files: {hiccup_data.atm_file}')
# print(f'    input sfc files: {hiccup_data.sfc_file}')
print(f'    input sst files: {hiccup_data.sst_file}')
print(f'    input ice files: {hiccup_data.ice_file}')
# print(f'    input topo file: {hiccup_data.topo_file}')
print('\n  Output files')
print(f'    output sst file: {output_sst_file_name}')

# Get dict of temporary files for each variable
file_dict = hiccup_data.get_multifile_dict(timestamp=999)

# ------------------------------------------------------------------------------
# Create SST/sea ice files

# create grid and mapping files
overwrite = False
hiccup_data.sstice_create_src_grid_file(force_overwrite=overwrite)
hiccup_data.sstice_create_dst_grid_file(force_overwrite=overwrite)
hiccup_data.sstice_create_map_file(force_overwrite=overwrite)

# Remap the sst/ice data after time slicing and combining (if necessary)
hiccup_data.sstice_slice_and_remap(output_file_name=output_sst_file_name,
                                   time_slice_method='initial')

# Rename the variables and remove unnecessary variables and attributes
hiccup_data.sstice_rename_vars(output_file_name=output_sst_file_name)

# Adjust final SST/ice data to fill in missing values and limit ice fraction
hiccup_data.sstice_adjustments(output_file_name=output_sst_file_name)

# ------------------------------------------------------------------------------
# Print final output file name

print()
print(f'output_sst_file_name: {output_sst_file_name}')
print()

# Print summary of timer info
hdc.print_timer_summary()

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
