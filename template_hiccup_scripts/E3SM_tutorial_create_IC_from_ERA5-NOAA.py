#!/usr/bin/env python
# ==================================================================================================
# HICCUP - Hindcast Initial Condition Creation Utility/Processor
# This tool automates the creation of atmospheric initial condition files for 
# E3SM using user supplied file for atmospheric and sea surface conditions.
# ==================================================================================================
import os
from hiccup import hiccup_data_class as hdc
from hiccup import hiccup_state_adjustment as hsa
# ------------------------------------------------------------------------------
# Logical flags for controlling what this script will do (comment out to disable)
verbose = True              # Global verbosity flag
unpack_nc_files   = True    # unpack data files (convert short to float)
create_map_file   = True    # grid and map file creation
remap_data_horz   = True    # horz remap, variable renaming
do_sfc_adjust     = True    # perform surface T and P adjustments
remap_data_vert   = True    # vertical remap
do_state_adjust   = True    # post vertical interpolation adjustments
do_random_perturb = True    # add random perturbations (1% of std dev)
combine_files     = True    # combine temporary data files and delete
create_sst_data   = True    # sst/sea ice file creation
# ------------------------------------------------------------------------------

hiccup_root = os.getenv('HOME')+'/HICCUP'
data_root = '/global/cfs/cdirs/e3sm/www/Tutorials/2024/practicum/day_4/atm_breakout/obs_data'

# Specify output atmosphere horizontal grid
dst_horz_grid = 'ne30np4'

# Specify output atmosphere vertical grid
dst_vert_grid,vert_file_name = 'L80',f'{hiccup_root}/files_vert/L80_for_E3SMv3.nc'

# Specify date of data (and separately specify year for SST/ice files)
init_date = '2023-09-08'
init_year = int(init_date.split('-')[0])

# Specify output file names
output_atm_file_name = f'{data_root}/HICCUP.atm_era5.{init_date}.{dst_horz_grid}.{dst_vert_grid}.nc'
output_sst_file_name = f'{data_root}/HICCUP.sst_noaa.{init_date}.nc'

# Specify topo file
topo_file_name = '/global/cfs/cdirs/e3sm/inputdata/atm/cam/topo/USGS-gtopo30_ne30np4pg2_x6t-SGH.c20210614.nc'

# ------------------------------------------------------------------------------
# Create data class instance, which includes xarray file dataset objects
# and variable name dictionaries for mapping between naming conventions.
# This also checks input files for required variables
hiccup_data = hdc.create_hiccup_data(name='ERA5'
                                    ,dst_horz_grid=dst_horz_grid
                                    ,dst_vert_grid=dst_vert_grid
                                    ,atm_file=f'{data_root}/ERA5.atm.{init_date}.00.nc'
                                    ,sfc_file=f'{data_root}/ERA5.sfc.{init_date}.00.nc'
                                    ,sstice_name='NOAA'
                                    ,sst_file=f'{data_root}/sst.day.mean.{init_year}.nc'
                                    ,ice_file=f'{data_root}/icec.day.mean.{init_year}.nc'
                                    ,topo_file=topo_file_name
                                    ,output_dir=data_root
                                    ,grid_dir=f'{data_root}/files_grid'
                                    ,map_dir=f'{data_root}/files_map'
                                    ,tmp_dir=f'{data_root}/files_tmp'
                                    ,verbose=verbose
                                    ,check_input_files=True)
# ------------------------------------------------------------------------------
# Print some informative stuff
print('\n  Input Files')
print(f'    input atm files: {hiccup_data.atm_file}')
print(f'    input sfc files: {hiccup_data.sfc_file}')
print(f'    input sst files: {hiccup_data.sst_file}')
print(f'    input ice files: {hiccup_data.ice_file}')
print(f'    input topo file: {hiccup_data.topo_file}')
print('\n  Output files')
print(f'    output atm file: {output_atm_file_name}')
print(f'    output sst file: {output_sst_file_name}')

# Get dict of temporary files for each variable
file_dict = hiccup_data.get_multifile_dict(timestamp=999)

# ------------------------------------------------------------------------------
# Make sure files are "unpacked" (may take awhile, so only do it if you need to)
if 'unpack_nc_files' not in locals(): unpack_nc_files = False
if unpack_nc_files:

    hiccup_data.unpack_data_files()

# ------------------------------------------------------------------------------
# Create grid and mapping files
if 'create_map_file' not in locals(): create_map_file = False
if create_map_file :

    # Create grid description files needed for the mapping file
    hiccup_data.create_src_grid_file()
    hiccup_data.create_dst_grid_file()

    # Create mapping file
    hiccup_data.create_map_file()

# ------------------------------------------------------------------------------
# perform multi-file horizontal remap
if 'remap_data_horz' not in locals(): remap_data_horz = False
if remap_data_horz :

    # Horizontally regrid the data
    hiccup_data.remap_horizontal_multifile(file_dict)

    # Rename variables to match what the model expects
    hiccup_data.rename_vars_multifile(file_dict=file_dict)

    # Add time/date information
    hiccup_data.add_time_date_variables_multifile(file_dict=file_dict)

# ------------------------------------------------------------------------------
# Do surface adjustments
if 'do_sfc_adjust' not in locals(): do_sfc_adjust = False
if do_sfc_adjust:

    hiccup_data.surface_adjustment_multifile(file_dict=file_dict)

# ------------------------------------------------------------------------------
# Vertically remap the data
if 'remap_data_vert' not in locals(): remap_data_vert = False
if remap_data_vert :

    hiccup_data.remap_vertical_multifile(file_dict=file_dict
                                        ,vert_file_name=vert_file_name)

# ------------------------------------------------------------------------------
# Perform final state adjustments on interpolated data and add additional data
if 'do_state_adjust' not in locals(): do_state_adjust = False
if do_state_adjust :

    hiccup_data.atmos_state_adjustment_multifile(file_dict=file_dict)

# ------------------------------------------------------------------------------
# Apply random perturbation to the final data
if 'do_random_perturb' not in locals(): do_random_perturb = False
if do_random_perturb :
    
    hiccup_data.atmos_state_apply_perturbations_multifile(file_dict=file_dict)

# ------------------------------------------------------------------------------
# Combine files
if 'combine_files' not in locals(): combine_files = False
if combine_files :

    # Combine and delete temporary files
    hiccup_data.combine_files(file_dict=file_dict
                             ,delete_files=True
                             ,output_file_name=output_atm_file_name)

    # Clean up the global attributes of the file
    hiccup_data.clean_global_attributes(file_name=output_atm_file_name)

# ------------------------------------------------------------------------------
# Create SST/sea ice file
if 'create_sst_data' not in locals(): create_sst_data = False
if create_sst_data :

    # create grid and mapping files
    overwrite = False
    hiccup_data.sstice_create_src_grid_file(force_overwrite=overwrite)
    hiccup_data.sstice_create_dst_grid_file(force_overwrite=overwrite)
    hiccup_data.sstice_create_map_file(force_overwrite=overwrite)

    # Remap the sst/ice data after time slicing and combining (if necessary)
    hiccup_data.sstice_slice_and_remap(output_file_name=output_sst_file_name,
                                       time_slice_method='match_atmos',
                                       atm_file=output_atm_file_name)

    # Rename the variables and remove unnecessary variables and attributes
    hiccup_data.sstice_rename_vars(output_file_name=output_sst_file_name)

    # Adjust final SST/ice data to fill in missing values and limit ice fraction
    hiccup_data.sstice_adjustments(output_file_name=output_sst_file_name)

# ------------------------------------------------------------------------------
# Print final output file names

print()
print(f'output_atm_file_name: {output_atm_file_name}')
print(f'output_sst_file_name: {output_sst_file_name}')
print()

# Print summary of timer info
hdc.print_timer_summary()

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
