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
# local path for grid and mapping files
hiccup_root  = os.getenv('HOME')+'/HICCUP'
data_root    = f'{hiccup_root}/test_data'
data_tmp     = f'{hiccup_root}/test_data_tmp'
# ------------------------------------------------------------------------------

os.makedirs(data_tmp, exist_ok=True)  # create temporary output data path if it doesn't exist

dst_horz_grid = 'ne30np4' # ne30np4 / ne120np4 / ne512np4 / ne1024np4

dst_vert_grid ='L128'; vert_file_name = f'{hiccup_root}/files_vert/vert_coord_E3SM_L128.nc' # SCREAM

# Specify output file names
output_atm_file_name = f'{data_tmp}/HICCUP_TEST_OUTPUT.eamxx_nudging_from_era5.{dst_horz_grid}.{dst_vert_grid}.nc'

topo_file = f'{data_root}/USGS-gtopo30_ne30np4_16xdel2-PFC-consistentSGH.nc'

# Create data class instance, which includes xarray file dataset objects
# and variable name dictionaries for mapping between naming conventions.
# This also checks input files for required variables
hiccup_data = hiccup.create_hiccup_data( src_data_name='ERA5',
                                         target_model='EAMXX-nudging',
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
# Vertically remap the data
hiccup_data.remap_vertical_multifile(file_dict=file_dict,
                                     vert_file_name=vert_file_name)
# Combine and delete temporary files
hiccup_data.combine_files(file_dict=file_dict,
                          delete_files=True,
                          output_file_name=output_atm_file_name)
# Clean up the global attributes of the file
hiccup_data.clean_global_attributes(file_name=output_atm_file_name)
# ------------------------------------------------------------------------------
# Print final output file name
print(); print(f'output_atm_file_name: {output_atm_file_name}'); print()
# Print summary of performance info
hiccup_data.print_memory_summary()
hiccup_data.print_timer_summary()
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
