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

### Path for data output
data_root = os.getenv('SCRATCH')+'/HICCUP/data/' # NERSC
# data_root = os.getenv('MEMBERWORK')+'/cli115/HICCUP/data/'  # OLCF

### Path for "supported" input data
inputdata_path = '/global/cfs/cdirs/e3sm/inputdata'

### time stamp for output file (= datetime.datetime.utcnow().strftime('%Y%m%d')])
timestamp = '20220707'

### output horizontal grid for atmosphere
dst_horz_grid = 'ne30np4'

### output vertical grid for atmosphere
dst_vert_grid,vert_file_name = 'L60',os.getenv('HOME')+f'/E3SM/vert_grid_files/L60.nc'

### specify input file name
cami_file = f'{inputdata_path}/atm/cam/inic/homme/cami_mam3_Linoz_ne30np4_L72_c160214.nc'

### specify output file
output_atm_file_name = f'{data_root}HICCUP.eam_i_mam3_Linoz_{dst_horz_grid}_{dst_vert_grid}_c{timestamp}.nc'

### topo file of output grid - replace this with file path if no default is set
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
                                    ,verbose=True)

### Print some informative stuff
print('\n  Input Files')
print(f'    input atm file:  {hiccup_data.atm_file}')
print(f'    input topo file: {hiccup_data.topo_file}')
print('\n  Output files')
print(f'    output atm file: {output_atm_file_name}')

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
