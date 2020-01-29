#!/usr/bin/env python
#===================================================================================================
# HICCUP - Hindcast Initial Condition Creation Utility/Processor
# This tool automates the creation of atmospheric initial condition files
# for E3SM using a user supplied reanalysis file, such as EAR5 data. 
# Requires NCO, TempestRemap, and xarray.
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

output_file_name = 'HICCUP_TEST.output.nc'
output_grid_name = 'ne30pg2'

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
hiccup_data = hdc.create_hiccup_data(name='ERA5'                \
                                    ,atm_file=input_file_atm    \
                                    ,sfc_file=input_file_sfc    \
                                    ,dst_grid_name=output_grid_name )

# Check input files for for required variables
hiccup_data.check_file_vars()

#-------------------------------------------------------------------------------
# Horizontally regrid the data
#-------------------------------------------------------------------------------

# Create grid description files needed for the mapping file
hiccup_data.create_src_grid_file()
hiccup_data.create_dst_grid_file()

# Create mapping file
hiccup_data.create_map_file()

# regrid the atm and sfc data to temporary files
# hiccup_data.remap_input_data()
atm_tmp_file_name = 'tmp_atm_data.nc'
sfc_tmp_file_name = 'tmp_sfc_data.nc'
sp.call(f'ncremap --map_file={hiccup_data.map_file} --in_file={hiccup_data.atm_file} --out_file={atm_tmp_file_name} ', shell=True)
sp.call(f'ncremap --map_file={hiccup_data.map_file} --in_file={hiccup_data.sfc_file} --out_file={sfc_tmp_file_name} ', shell=True)

# Combine the temporary files into the final output file
sp.call(f'ncks -A {atm_tmp_file_name} {output_file_name} ', shell=True)
sp.call(f'ncks -A {sfc_tmp_file_name} {output_file_name} ', shell=True)

# delete the temporary files
sp.call(f'rm {sfc_tmp_file_name} {atm_tmp_file_name} ', shell=True)

# Rename variables to match what the model expects
hiccup_data.rename_vars(output_file_name)

print(f'\noutput_file_name: {output_file_name}\n')

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


