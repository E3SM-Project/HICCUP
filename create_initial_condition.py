#!/usr/bin/env python
#===================================================================================================
# HICCUP - Hindcast Initial Condition Creation Utility/Processor
# This tool automates the creation of atmospheric initial condition files
# for E3SM using a user supplied reanalysis file, such as EAR5 data. 
# Requires NCO and xarray.
#===================================================================================================
import os
import subprocess as sp
import glob
import datetime
import xarray as xr
import hiccup_data_class as hdc
import hiccup_state_adjustment
#-------------------------------------------------------------------------------
# Specify file names
#-------------------------------------------------------------------------------

input_file_atm = 'ERA5.HICCUP_TEST.atm.remap.nc'
input_file_sfc = 'ERA5.HICCUP_TEST.sfc.remap.nc'

output_file_name = 'HICCUP.output.nc'
output_grid_name = 'ne30np4'

tmp_file_name = 'tmp.nc'

# Options for output state adjustment
adjust_ts   = False   # Adjust surface temperature to match new surface height
adjust_ps   = False   # Adjust surface pressure to match new surface height
adjust_mass = False   # adjust surface pressure to retain dry mass of atmosphere
adjust_qv   = False   # adjust qv to eliminate supersaturation
adjust_cw   = False   # adjust cloud water to remove negative values
adjust_cf   = False   # adjust cloud fraction to remove values outside of [0,1]

#-------------------------------------------------------------------------------
# Load input data 
#-------------------------------------------------------------------------------

# Create data class instance, which includes xarray file dataset objects 
# and variable name dictionaries for mapping between naming conventions
hiccup_data = hdc.create_hiccup_data(name='ERA5',atm_file=input_file_atm,sfc_file=input_file_sfc)


# print(hiccup_data.ds_atm)
# exit()

#-------------------------------------------------------------------------------
# Check input files for for required variables
#-------------------------------------------------------------------------------

# Create list of variables in the files
atm_file_vars = []
sfc_file_vars = []
for key in hiccup_data.ds_atm.variables.keys(): atm_file_vars.append(key)
for key in hiccup_data.ds_sfc.variables.keys(): sfc_file_vars.append(key)

# Check that all required data exists in the atm file
for key in hiccup_data.atm_var_name_dict : 
    if hiccup_data.atm_var_name_dict[key] not in atm_file_vars: 
        raise ValueError(f'{hiccup_data.atm_var_name_dict[key]} is not in ATM dataset: ({input_file_atm})')

# Check that all required data exists in the sfc file
for key in hiccup_data.sfc_var_name_dict : 
    if hiccup_data.sfc_var_name_dict[key] not in sfc_file_vars: 
        raise ValueError(f'{hiccup_data.sfc_var_name_dict[key]} is not in SFC dataset: ({input_file_sfc})')

#-------------------------------------------------------------------------------
# Make a copy of the input file and rename/subset variables
#-------------------------------------------------------------------------------

# os.system(f'cp {input_file_name} {tmp_file_name}')

# Rename the variables to match the output names
# for key in var_rename_dict :
#   os.system(f'ncrename -v {key},{var_rename_dict[key]}  {tmp_file_name} ')

# Insert new PS variable
# os.system(f'ncap2 -s "PS=PHIS" {tmp_file_name} ')

# Insert new P0 variable
# os.system(f'ncap2 -s "P0=1.0D0" {tmp_file_name} ')

#-------------------------------------------------------------------------------
# Horizontally regrid the data
#-------------------------------------------------------------------------------

hiccup_data.create_src_grid_file()
hiccup_data.create_dst_grid_file(grid=output_grid_name)

print(hiccup_data.grid_file)

# Create mapping file
# src_grid = '????'
# dst_grid = 'ne30np4'
# map_file = f'map_{src_grid}_to_{dst_grid}_aave.nc'

# os.system(f' ncks --map {map_file}  {file_in}  {file_out} ')

#-------------------------------------------------------------------------------
# Prepare the vertical grid file for vertical regridding
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Vertically regrid the data
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Perform state adjustments on interpolated data
#-------------------------------------------------------------------------------

exit()

# Load the file into an xarray dataset
ds = xr.open_dataset(output_file_name)

# Adjust surface temperature to match new surface height
if adjust_ts : 
  state_adjustment.adjust_surface_temperature( ncol, phis_old, ts_old, phis_new, ts_new )

# Adjust surface pressure to match new surface height
if adjust_ps : 
  state_adjustment.adjust_surface_pressure( plev, ncol, temperature_mid,  \
                                            pressure_mid, pressure_int,   \
                                            phis_old, ps_old, phis_new, ps_new )

# adjust surface pressure to retain dry mass of atmosphere?
# if adjust_mass : 
#   state_adjustment.dry_mass_fixer( ncol, plev, hyai, hybi, wgt, qv, mass_ref, ps_in, ps_out )

# adjust qv to eliminate supersaturation?
if adjust_qv :
  state_adjustment.remove_supersaturation()

# adjust cloud water to remove negative values?
if adjust_cw :
  cld_liq = cld_liq.where(cld_liq>=0,other=0.)
  cld_ice = cld_ice.where(cld_ice>=0,other=0.)

# adjust cloud fraction to remove values outside of [0,1]
if adjust_cf :
  state_adjustment.adjust_cloud_fraction()

# Write the dataset back to the file
ds.to_netcdf(output_file_name)

#-------------------------------------------------------------------------------
# Clean up
#-------------------------------------------------------------------------------


