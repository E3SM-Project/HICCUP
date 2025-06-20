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
# Parse the command line options
parser = optparse.OptionParser()
parser.add_option('--hgrid',dest='horz_grid',default=None,help='Sets the output horizontal grid')
parser.add_option('--vgrid',dest='vert_grid',default=None,help='Sets the output vertical grid')
(opts, args) = parser.parse_args()
# ------------------------------------------------------------------------------
# Logical flags for controlling what this script will do (comment out to disable)
create_map_file = True    # grid and map file creation
remap_data_horz = True    # horz remap, variable renaming
do_sfc_adjust   = True    # perform surface adjustments
do_state_adjust = True    # post vertical interpolation adjustments
combine_files   = True    # combine temporary data files and delete
# ------------------------------------------------------------------------------

output_root   = os.getenv('SCRATCH')+'/HICCUP' # root path for HICCUP output
dst_horz_grid = 'ne16np4'                      # output horizontal grid for atmosphere
dst_vert_grid = 'L72'                          # output vertical grid for atmosphere
timestamp     = '99999999'                     # time stamp for output file

# specify input file name
src_eami_file = f'{output_root}/files_init/eami_mam4_Linoz_ne30np4_L80_c20231010_w-phis.nc'
# src_eami_file = '/global/cfs/cdirs/e3sm/inputdata/atm/cam/inic/homme/cami_mam3_Linoz_ne30np4_L72_c160214.nc'

# specify output file
dst_eami_file = f'{output_root}/files_init/HICCUP.eam_i_mam3_Linoz_{dst_horz_grid}_{dst_vert_grid}_c{timestamp}.nc'

# topo file of output grid - replace this with file path if no default is set
topo_file_name = hdc.get_default_topo_file_name(dst_horz_grid)

# ------------------------------------------------------------------------------
# Create HICCUP data class instance

# this includes xarray file dataset objects and variable 
# name dictionaries for mapping between naming conventions.
# This also checks input files for required variables
hiccup_data = hiccup.create_hiccup_data(src_data_name='EAM',
                                        target_model='EAM', # options: EAM / EAMXX
                                        dst_horz_grid=dst_horz_grid,
                                        dst_vert_grid=dst_vert_grid,
                                        atm_file=src_eami_file,
                                        sfc_file=src_eami_file,
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
print(f'    output atm file: {dst_eami_file}')

# ------------------------------------------------------------------------------
# create file name dictionaries for np4 and pg2 data (both are needed)

dst_horz_grid_pg = dst_horz_grid.replace('np4','pg2')

# use separate timestamp for temporary files to ensure these files 
# are distinct from other instances of HICCUP that might be running.
# The routine to generate the multifile dict will generate it's own
# timestamp, but it wil use minutes and seconds, which can become a 
# problem when iteratively debugging a HICCUP script
tmp_timestamp = datetime.datetime.utcnow().strftime('%Y%m%d')

# Get dict of temporary files for each variable
file_dict_np = hiccup_data.get_multifile_dict_eam(dst_horz_grid=hiccup_data.dst_horz_grid_np,
                                                  var_name_dict=hiccup_data.atm_var_name_dict_np,
                                                  timestamp=tmp_timestamp)
file_dict_pg = hiccup_data.get_multifile_dict_eam(dst_horz_grid=hiccup_data.dst_horz_grid_pg,
                                                  var_name_dict=hiccup_data.atm_var_name_dict_pg,
                                                  timestamp=tmp_timestamp)

# ------------------------------------------------------------------------------
# Create grid and mapping files
if 'create_map_file' not in locals(): create_map_file = False
if create_map_file :

    # Create grid description files needed for the mapping file
    hiccup_data.create_src_grid_file()
    hiccup_data.create_dst_grid_file()

    # Create mapping files
    hiccup_data.create_map_file(src_type='GLL')

# ------------------------------------------------------------------------------
# perform multi-file horizontal remap
if 'remap_data_horz' not in locals(): remap_data_horz = False
if remap_data_horz :

    # Horizontally regrid np4 data
    hiccup_data.map_file = hiccup_data.map_file_np
    hiccup_data.remap_horizontal_multifile_eam(file_dict=file_dict_np)

    # Horizontally regrid pg2 data
    hiccup_data.map_file = hiccup_data.map_file_pg 
    hiccup_data.remap_horizontal_multifile_eam(file_dict=file_dict_pg)

    # Add time/date information
    hiccup_data.add_time_date_variables_multifile_eam(file_dict=file_dict_np)
    hiccup_data.add_time_date_variables_multifile_eam(file_dict=file_dict_pg)

# ------------------------------------------------------------------------------
# Do surface adjustments
if 'do_sfc_adjust' not in locals(): do_sfc_adjust = False
if do_sfc_adjust:

    hiccup_data.surface_adjustment_multifile(file_dict=file_dict_np,
                                             adj_TS=False,adj_PS=True)

# ------------------------------------------------------------------------------
# Perform final state adjustments on interpolated data and add additional data
if 'do_state_adjust' not in locals(): do_state_adjust = False
if do_state_adjust :

    # only need to adjust data on np4 grid
    hiccup_data.atmos_state_adjustment_multifile(file_dict=file_dict_np)

# ------------------------------------------------------------------------------
# Combine files
if 'combine_files' not in locals(): combine_files = False
if combine_files :

    # combine file names into a single dict
    file_dict_all = {}
    file_dict_all.update(file_dict_np)
    file_dict_all.update(file_dict_pg)

    # Combine and delete temporary files
    hiccup_data.combine_files(file_dict=file_dict_all
                             ,delete_files=True
                             ,output_file_name=dst_eami_file)

    # Clean up the global attributes of the file
    hiccup_data.clean_global_attributes(file_name=dst_eami_file)

# ------------------------------------------------------------------------------
# Print final output file name

print()
print(f'dst_eami_file: {dst_eami_file}')
print()

# Print summary of timer info
hiccup_data.print_timer_summary()

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
