#!/usr/bin/env python
# ==================================================================================================
# HICCUP - Hindcast Initial Condition Creation Utility/Processor
# This tool automates the creation of atmospheric initial condition files for 
# E3SM using user supplied file for atmospheric and sea surface conditions.
# ==================================================================================================
'''
# NOTE - the standard "eam.i" file does not include the PHIS variable, which
# is a problem for performing the surface adjustment step. To address this
# I've put example commands below to add the PHIS variable from the original
# topography file into the initial condition file, so that it can be remapped
# to the new grid and used for the surface adjustment.

DIN_LOC_ROOT=/global/cfs/cdirs/e3sm/inputdata
src_init_root=${DIN_LOC_ROOT}/atm/cam/inic/homme
dst_init_root=/global/cfs/cdirs/m4310/whannah/files_init
src_init_file=${src_init_root}/eami_mam4_Linoz_ne30np4_L80_c20231010.nc
dst_init_file=${dst_init_root}/eami_mam4_Linoz_ne30np4_L80_c20231010_w-phis.nc
src_topo_file=${DIN_LOC_ROOT}/atm/cam/topo/USGS-gtopo30_ne30np4pg2_16xdel2.c20200108.nc

cp ${src_init_file}  ${dst_init_file}
ncks -A ${src_topo_file} ${dst_init_file} -v PHIS_d
ncrename -v PHIS_d,PHIS  ${dst_init_file}
echo ; echo ${dst_init_file} ; echo

'''
# ==================================================================================================
import os, optparse, datetime
from hiccup import hiccup
# ------------------------------------------------------------------------------
# Logical flags for controlling what this script will do (comment out to disable)
create_map_file = True    # grid and map file creation
remap_data_horz = True    # horz remap, variable renaming
do_sfc_adjust   = True    # perform surface T and P adjustments
remap_data_vert = True    # vertical remap
do_state_adjust = True    # post vertical interpolation adjustments
combine_files   = True    # combine temporary data files and delete
# ------------------------------------------------------------------------------

hiccup_root   = os.getenv('HOME')   +'/HICCUP' # local HICCUP path
output_root   = os.getenv('SCRATCH')+'/HICCUP' # root path for HICCUP output
dst_horz_grid = 'ne128np4'                     # output horizontal grid for atmosphere
dst_vert_grid = 'L128'                         # output vertical grid for atmosphere
dst_vert_file = f'{hiccup_root}/files_vert/vert_coord_E3SM_L128.nc'
timestamp     = '99999999'                     # time stamp for output file

# Path for supported E3SM input data
inputdata_root = '/global/cfs/cdirs/e3sm/inputdata'

# specify input EAM file name
eami_file = '/global/cfs/cdirs/m4310/whannah/E3SM/init_data/v3.LR.amip_0101/archive/rest/2000-01-01-00000/v3.LR.amip_0101.eam.i.2000-01-01-00000_w-topo-x6t-SGH.nc'

# specify output EAMxx file
output_atm_file_name = f'{output_root}/files_init/v3.LR.amip_0101.eam.i.2000-01-01-00000.EAMxx-format.{dst_horz_grid}.{timestamp}.nc'

# topo file of output grid
topo_file_name = f'{inputdata_root}/atm/cam/topo/USGS-topo_{dst_horz_grid}_smoothedx6t_20250904.nc'

# ------------------------------------------------------------------------------
# Create HICCUP data class instance

# this includes xarray file dataset objects and variable 
# name dictionaries for mapping between naming conventions.
# This also checks input files for required variables
hiccup_data = hiccup.create_hiccup_data(src_data_name='EAM',
                                        target_model='EAMXX', # options: EAM / EAMXX
                                        dst_horz_grid=dst_horz_grid,
                                        dst_vert_grid=dst_vert_grid,
                                        atm_file=eami_file,
                                        sfc_file=eami_file,
                                        topo_file=topo_file_name,
                                        output_dir=f'{output_root}/files_init',
                                        grid_dir=f'{output_root}/files_grid',
                                        map_dir=f'{output_root}/files_map',
                                        tmp_dir=f'{output_root}/files_hiccup_tmp',
                                        verbose=True,)

# Print some informative stuff
print('\n  Input Files')
print(f'    input atm file:  {hiccup_data.atm_file}')
print(f'    input topo file: {hiccup_data.topo_file}')
print('\n  Output files')
print(f'    output atm file: {output_atm_file_name}')

# ------------------------------------------------------------------------------
# Get dict of temporary files for each variable
file_dict = hiccup_data.get_multifile_dict(timestamp=timestamp)
# ------------------------------------------------------------------------------
# Create grid and mapping files
if 'create_map_file' not in locals(): create_map_file = False
if create_map_file :

    # Create grid description files needed for the mapping file
    hiccup_data.create_src_grid_file()
    hiccup_data.create_dst_grid_file()

    # Create mapping files
    hiccup_data.create_map_file()

# ------------------------------------------------------------------------------
# perform multi-file horizontal remap
if 'remap_data_horz' not in locals(): remap_data_horz = False
if remap_data_horz :

    # Horizontally regrid np4 data
    hiccup_data.map_file = hiccup_data.map_file_np
    hiccup_data.remap_horizontal_multifile(file_dict=file_dict)

    # Rename variables to match what the model expects
    hiccup_data.rename_vars_multifile(file_dict=file_dict)

    # Add time/date information
    hiccup_data.add_time_date_variables_multifile_eam(file_dict=file_dict)
# ------------------------------------------------------------------------------
# Do surface adjustments
if 'do_sfc_adjust' not in locals(): do_sfc_adjust = False
if do_sfc_adjust:

    hiccup_data.surface_adjustment_multifile(file_dict=file_dict,
                                             adj_TS=False,
                                             adj_PS=True,
                                             adj_T_eam=True)

# ------------------------------------------------------------------------------
# Vertically remap the data
if 'remap_data_vert' not in locals(): remap_data_vert = False
if remap_data_vert :

  hiccup_data.remap_vertical_multifile(file_dict=file_dict,
                                       vert_file_name=dst_vert_file)

# ------------------------------------------------------------------------------
# Perform final state adjustments on interpolated data and add additional data
if 'do_state_adjust' not in locals(): do_state_adjust = False
if do_state_adjust :

    # only need to adjust data on np4 grid
    hiccup_data.atmos_state_adjustment_multifile(file_dict=file_dict)

# ------------------------------------------------------------------------------
# Combine files
if 'combine_files' not in locals(): combine_files = False
if combine_files :

    # Combine and delete temporary files
    hiccup_data.combine_files(file_dict=file_dict,
                              delete_files=True,
                              output_file_name=output_atm_file_name)

    # Clean up the global attributes of the file
    hiccup_data.clean_global_attributes(file_name=output_atm_file_name)

# ------------------------------------------------------------------------------
# Print final output file name

print()
print(f'output_atm_file_name: {output_atm_file_name}')
print()

# Print summary of timer info
hiccup_data.print_timer_summary()

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
