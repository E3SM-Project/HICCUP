#!/usr/bin/env python
# ==================================================================================================
# HICCUP - Hindcast Initial Condition Creation Utility/Processor
# This tool automates the creation of atmospheric initial condition files for 
# E3SM using user supplied file for atmospheric and sea surface conditions.
# ==================================================================================================
import os, optparse
from hiccup import hiccup_data_class as hdc
# ------------------------------------------------------------------------------
# Parse the command line options
parser = optparse.OptionParser()
parser.add_option('--hgrid',dest='horz_grid',default=None,help='Sets the output horizontal grid')
parser.add_option('--vgrid',dest='vert_grid',default=None,help='Sets the output vertical grid')
(opts, args) = parser.parse_args()
# ------------------------------------------------------------------------------
# Logical flags for controlling what this script will do
create_map_file = True    # grid and map file creation
remap_data_horz = True    # horz remap, variable renaming
remap_data_vert = True    # vertical remap
do_state_adjust = True    # post vertical interpolation adjustments
combine_files   = True    # combine temporary data files and delete
verbose = True            # Global verbosity flag
# ------------------------------------------------------------------------------

# Path for data output
data_root = os.getenv('SCRATCH')+'/HICCUP/data/' # NERSC
# data_root = os.getenv('MEMBERWORK')+'/cli115/HICCUP/data/'  # OLCF

# Path for supported E3SM input data
inputdata_path = '/global/cfs/cdirs/e3sm/inputdata'

# time stamp for output file (= datetime.datetime.utcnow().strftime('%Y%m%d')])
timestamp = '20220707'

# output horizontal grid for atmosphere
dst_horz_grid = 'ne16np4'

# output vertical grid for atmosphere
dst_vert_grid,vert_file_name = 'L60',os.getenv('HOME')+f'/E3SM/vert_grid_files/L60.nc'

# specify input file name
cami_file = f'{inputdata_path}/atm/cam/inic/homme/cami_mam3_Linoz_ne30np4_L72_c160214.nc'

# specify output file
output_atm_file_name = f'{data_root}HICCUP.eam_i_mam3_Linoz_{dst_horz_grid}_{dst_vert_grid}_c{timestamp}.nc'

# topo file of output grid - replace this with file path if no default is set
topo_file_name = hdc.get_default_topo_file_name(dst_horz_grid)

# ------------------------------------------------------------------------------
# Create HICCUP data class instance
# ------------------------------------------------------------------------------
# this includes xarray file dataset objects and variable 
# name dictionaries for mapping between naming conventions.
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

# Print some informative stuff
print('\n  Input Files')
print(f'    input atm file:  {hiccup_data.atm_file}')
print(f'    input topo file: {hiccup_data.topo_file}')
print('\n  Output files')
print(f'    output atm file: {output_atm_file_name}')

# ------------------------------------------------------------------------------
# create file name dictionaries for np4 and pg2 data (both are needed)
# ------------------------------------------------------------------------------

dst_horz_grid_pg = dst_horz_grid.replace('np4','pg2')

# Get dict of temporary files for each variable
file_dict_np4 = hiccup_data.get_multifile_dict_eam(dst_horz_grid=dst_horz_grid   ,var_name_dict=hiccup_data.atm_var_name_dict)
file_dict_pg2 = hiccup_data.get_multifile_dict_eam(dst_horz_grid=dst_horz_grid_pg,var_name_dict=hiccup_data.atm_var_name_dict_pg2)

file_dict_all = {}
file_dict_all.update(file_dict_np4)
file_dict_all.update(file_dict_pg2)

# ------------------------------------------------------------------------------
# Create grid and mapping files
# ------------------------------------------------------------------------------
if 'create_map_file' not in locals(): create_map_file = False
# if create_map_file :

    # Create grid description files needed for the mapping file
    # hiccup_data.create_src_grid_file()
    # hiccup_data.create_dst_grid_file()
    
    # hiccup_data.src_grid_file = '/global/homes/w/whannah/E3SM/data_grid/ne30.g'
    # hiccup_data.dst_grid_file = '/global/homes/w/whannah/E3SM/data_grid/ne16.g'

    # # Create mapping file
    # hiccup_data.create_map_file()

# ncremap --alg_typ=tempest --src_grd=${HOME}/E3SM/data_grid/ne45.g           --dst_grd=${HOME}/E3SM/data_grid/ne16.g           --map_file=${HOME}/maps/map_ne45np4_to_ne16np4_mono_20220707.nc --wgt_opt=' --in_type cgll --in_np 4 --out_type cgll --out_np 4  --out_double '
# ncremap --alg_typ=tempest --src_grd=${HOME}/E3SM/data_grid/ne45pg2_scrip.nc --dst_grd=${HOME}/E3SM/data_grid/ne16pg2_scrip.nc --map_file=${HOME}/maps/map_ne45pg2_to_ne16pg2_aave_20220707.nc --wgt_opt=' --in_type fv   --in_np 2 --out_type fv   --out_np 2  --out_double '

map_file_np = os.getenv('HOME')+f'/maps/RRM/map_ne45np4_to_{dst_horz_grid   }_mono_20220707.nc'
map_file_pg = os.getenv('HOME')+f'/maps/RRM/map_ne45pg2_to_{dst_horz_grid_pg}_aave_20220707.nc'

# ------------------------------------------------------------------------------
# perform multi-file horizontal remap
# ------------------------------------------------------------------------------
if 'remap_data_horz' not in locals(): remap_data_horz = False
if remap_data_horz :

    # Horizontally regrid np4 data
    hiccup_data.map_file = map_file_np
    hiccup_data.remap_horizontal_multifile_eam(file_dict=file_dict_np4)

    # Horizontally regrid pg2 data
    hiccup_data.map_file = map_file_pg 
    hiccup_data.remap_horizontal_multifile_eam(file_dict=file_dict_pg2)

    # Add time/date information
    hiccup_data.add_time_date_variables_multifile_eam(file_dict=file_dict_np4)
    hiccup_data.add_time_date_variables_multifile_eam(file_dict=file_dict_pg2)

# ------------------------------------------------------------------------------
# Perform final state adjustments on interpolated data and add additional data
# ------------------------------------------------------------------------------
if 'do_state_adjust' not in locals(): do_state_adjust = False
if do_state_adjust :

    # only adjust data on np4 grid
    hiccup_data.atmos_state_adjustment_multifile(file_dict=file_dict_np4)

# ------------------------------------------------------------------------------
# Vertically remap the data
# ------------------------------------------------------------------------------
if 'remap_data_vert' not in locals(): remap_data_vert = False
if remap_data_vert :

  hiccup_data.remap_vertical_multifile(file_dict=file_dict_np4,
                                       vert_file_name=vert_file_name):

  hiccup_data.remap_vertical_multifile(file_dict=file_dict_pg2,
                                       vert_file_name=vert_file_name):

# ------------------------------------------------------------------------------
# Combine files
# ------------------------------------------------------------------------------
if 'combine_files' not in locals(): combine_files = False
if combine_files :

    # Combine and delete temporary files
    hiccup_data.combine_files(file_dict=file_dict_all
                             ,delete_files=True
                             ,output_file_name=output_atm_file_name)

    # Clean up the global attributes of the file
    hiccup_data.clean_global_attributes(file_name=output_atm_file_name)

# ------------------------------------------------------------------------------
# Print final output file name
# ------------------------------------------------------------------------------

print()
print(f'output_atm_file_name: {output_atm_file_name}')
print()

# Print summary of timer info
hdc.print_timer_summary()

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
