#!/usr/bin/env python
#===================================================================================================
# This script automates the creation of an E3SM atmospheric initial condition 
# file from an ERA5 reanalysis file. NCO must be installed locally and xarray 
# must be installed in the python environment. 
#===================================================================================================
import os
import subprocess as sp
import glob
import datetime
import xarray as xr
from . import state_adjustments
#===============================================================================
# Specify file names
#===============================================================================

input_file_name = 'input.nc'

output_file_name = 'output.nc'

tmp_file_name = 'tmp.nc'

# This dict provides a map between variable names in the input and output files
var_rename_dict = {}

# DYAMOND - hybrid coords renamed because they use a different dimension name
var_rename_dict.update({'hyam':'ohyam','hybm':'ohybm','hyai':'ohyai','hybi':'ohybi'})
var_rename_dict.update({'t':'T','q':'Q','u':'U','v':'V'})
var_rename_dict.update({'clwc':'CLDLIQ','ciwc':'CLDICE','z_2':'PHIS','o3':'O3'})
var_rename_dict.update({'stl1':'TS1','stl2':'TS2','stl3':'TS3','stl4':'TS4'})

#===============================================================================
# Check that input file has required variables
#===============================================================================


#===============================================================================
# Make a copy of the input file and rename/subset variables
#===============================================================================
os.system(f'cp {input_file_name} {tmp_file_name}'

# Rename the variables to match the output names
# for key in var_rename_dict :
#   os.system(f'ncrename -v {key},{var_rename_dict[key]}  {tmp_file_name} ')


# Insert new PS variable
# os.system(f'ncap2 -s "PS=PHIS" {tmp_file_name} ')

# Insert new P0 variable
# os.system(f'ncap2 -s "P0=1.0D0" {tmp_file_name} ')

#===============================================================================
# Horizontally regrid the data
#===============================================================================
# Create mapping file
src_grid = '????'
dst_grid = 'ne30np4'
map_file = f'map_{src_grid}_to_{dst_grid}_aave.nc'

# os.system(f' ncks --map {map_file}  {file_in}  {file_out} ')

#===============================================================================
# Prepare the vertical grid file for vertical regridding
#===============================================================================

#===============================================================================
# Vertically regrid the data
#===============================================================================

#===============================================================================
# Adjust surface pressure and temperature to match new surface height
#===============================================================================

# state_adjustments.adjust_surface_state(plev, ncol, temperature,     \
#                                        pressure_mid, pressure_int,  \
#                                        phis_old, ps_old, ts_old,    \
#                                        phis_new, ps_new, ts_new))

#===============================================================================
# Perform additional state variable adjustments
#===============================================================================


#===============================================================================
# Clean up
#===============================================================================

#===============================================================================
#===============================================================================
