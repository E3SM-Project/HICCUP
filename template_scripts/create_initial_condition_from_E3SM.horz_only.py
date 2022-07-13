#!/usr/bin/env python
# ==================================================================================================
# HICCUP - Hindcast Initial Condition Creation Utility/Processor
# This tool automates the creation of atmospheric initial condition files for 
# E3SM using user supplied file for atmospheric and sea surface conditions.
# ==================================================================================================
import os, optparse, datetime
from hiccup import hiccup_data_class as hdc
# from hiccup import hiccup_state_adjustment as hsa # this gets imported as part of hdc
# ------------------------------------------------------------------------------
# Parse the command line options
parser = optparse.OptionParser()
parser.add_option('--hgrid',dest='horz_grid',default=None,help='Sets the output horizontal grid')
parser.add_option('--vgrid',dest='vert_grid',default=None,help='Sets the output vertical grid')
(opts, args) = parser.parse_args()
# ------------------------------------------------------------------------------
# Logical flags for controlling what this script will do
# create_map_file = True    # grid and map file creation
remap_data_horz = True    # horz remap, variable renaming
# remap_data_vert = True    # vertical remap
# do_state_adjust = True    # post vertical interpolation adjustments
combine_files   = True    # combine temporary data files and delete
verbose = True            # Global verbosity flag
# ------------------------------------------------------------------------------

### Path for data output
data_root = os.getenv('SCRATCH')+'/HICCUP/data/' # NERSC
# data_root = os.getenv('MEMBERWORK')+'/cli115/HICCUP/data/'  # OLCF

### Path for "supported" input data
scratch_ic_path = '/global/cfs/cdirs/e3sm/inputdata/atm/cam/inic/homme'

### time stamp for output file (= datetime.datetime.utcnow().strftime('%Y%m%d')])
timestamp = '20210907'

### output horizontal grid for atmosphere
# dst_horz_grid = 'ne4np4'
# dst_horz_grid = 'ne45np4'
# dst_horz_grid = 'ne120np4'
dst_horz_grid = 'RRM-cubeface-grad'

### output vertical grid for atmosphere
dst_vert_grid,vert_file_name = 'L60',os.getenv('HOME')+f'/E3SM/vert_grid_files/L60.nc'


# L60 spinup
# data_root = '/global/cscratch1/sd/whannah/e3sm_scratch/init_files/L60_initial_conditions'
# cami_file = f'{data_root}/eam_i_mam3_Linoz_ne30np4_{dst_vert_grid}_c20210901.nc'

scratch_ic_path = '/global/homes/w/whannah/E3SM/init_scratch/L60_initial_conditions'
cami_file = f'{scratch_ic_path}/eam_i_aquaplanet_ne45np4_L60_c20210823.nc'


if '_mam3_' in cami_file: output_atm_file_name = f'{data_root}/eam_i_mam3_Linoz_{dst_horz_grid}_{dst_vert_grid}_c{timestamp}.nc'
if 'aqua'   in cami_file: output_atm_file_name = f'{data_root}/eam_i_aqua_{dst_horz_grid}_{dst_vert_grid}_c{timestamp}.nc'
if 'rcemip' in cami_file: output_atm_file_name = f'{data_root}/eam_i_rcemip_{dst_horz_grid}_{dst_vert_grid}_c{timestamp}.nc'

### topo file of output grid - replace this with file path if no defalt is set
topo_file_name = hdc.get_default_topo_file_name('ne45np4')
# topo_file_name = hdc.get_default_topo_file_name(dst_horz_grid)
# topo_file_name = '/project/projectdirs/e3sm/inputdata/atm/cam/topo/USGS_conusx4v1-tensor12x_consistentSGH_c150924.nc'

### Create data class instance, which includes xarray file dataset objects
### and variable name dictionaries for mapping between naming conventions.
### This also checks input files for required variables
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

### Print some informative stuff
print('\n  Input Files')
print(f'    input atm file:  {hiccup_data.atm_file}')
print(f'    input topo file: {hiccup_data.topo_file}')
# print(f'    input sfc file:  {hiccup_data.sfc_file}')
# print(f'    input sst file:  {hiccup_data.sst_file}')
# print(f'    input ice file:  {hiccup_data.ice_file}')
print('\n  Output files')
print(f'    output atm file: {output_atm_file_name}')
# print(f'    output sst file: {output_sst_file_name}')

# exit()

dst_horz_grid_pg = dst_horz_grid.replace('np4','pg2')

### Get dict of temporary files for each variable
# file_dict = hiccup_data.get_multifile_dict(eam=True)
file_dict_np4 = hiccup_data.get_multifile_dict_eam(dst_horz_grid=dst_horz_grid   ,var_name_dict=hiccup_data.atm_var_name_dict)
file_dict_pg2 = hiccup_data.get_multifile_dict_eam(dst_horz_grid=dst_horz_grid_pg,var_name_dict=hiccup_data.atm_var_name_dict_pg2)

file_dict_all = {}
file_dict_all.update(file_dict_np4)
file_dict_all.update(file_dict_pg2)


# ------------------------------------------------------------------------------
# Create grid and mapping files
# ------------------------------------------------------------------------------
# if 'create_map_file' not in locals(): create_map_file = False
# if create_map_file :

#     # Create grid description files needed for the mapping file
#     hiccup_data.create_src_grid_file()
#     hiccup_data.create_dst_grid_file()
#     # hiccup_data.dst_grid_file = '/global/homes/w/whannah/E3SM/data_grid/conusx4v1.g'

#     # Create mapping file
#     hiccup_data.create_map_file()

map_file_np = os.getenv('HOME')+f'/maps/RRM/map_ne45np4_to_{dst_horz_grid}_mono_20210901.nc'
map_file_pg = os.getenv('HOME')+f'/maps/RRM/map_ne45pg2_to_{dst_horz_grid_pg}_aave_20210901.nc'

# ------------------------------------------------------------------------------
# perform multi-file horizontal remap
# ------------------------------------------------------------------------------
if 'remap_data_horz' not in locals(): remap_data_horz = False
if remap_data_horz :

    ### Horizontally regrid np4 data
    hiccup_data.map_file = map_file_np
    hiccup_data.remap_horizontal_multifile_eam(file_dict=file_dict_np4)

    ### Horizontally regrid pg2 data
    hiccup_data.map_file = map_file_pg 
    hiccup_data.remap_horizontal_multifile_eam(file_dict=file_dict_pg2)

    ### Add time/date information
    hiccup_data.add_time_date_variables_multifile_eam(file_dict=file_dict_np4)
    hiccup_data.add_time_date_variables_multifile_eam(file_dict=file_dict_pg2)

# ------------------------------------------------------------------------------
# Perform final state adjustments on interpolated data and add additional data
# ------------------------------------------------------------------------------
if 'do_state_adjust' not in locals(): do_state_adjust = False
if do_state_adjust :

    ### only adjust data on np4 grid
    hiccup_data.atmos_state_adjustment_multifile(file_dict=file_dict_np4)

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
# print(f'output_sst_file_name: {output_sst_file_name}')
print()

# Print summary of timer info
hdc.print_timer_summary()

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
