#!/usr/bin/env python
# ==================================================================================================
# HICCUP - Hindcast Initial Condition Creation Utility/Processor
# This tool automates the creation of atmospheric initial condition files
# for E3SM using a user supplied reanalysis file, such as EAR5 data.
# Requires NCO, TempestRemap, numpy, and xarray.
# 
# NOTE ABOUT VERTICAL GRID FILES: 
# The current E3SM vertical grid was created through an iterative process 
# involving numerous, undocumented, subjective decisions mainly by Phil Rasch 
# and Po-Lun Ma who did not document the process, so there is no recipe to 
# recreate the grid from scratch. To create the vertical coordinate file it is 
# easiest to extract it from a pre-existing model data file as follows:
#   1. Dump the vertical grid data into a text file using ncdump:
#      ncdump -v P0,hyam,hybm,hyai,hybi,lev,ilev <history_file> > vert_coord.txt
#   2. manually edit the file to remove extra header info,
#      but keep the general CDL format created by ncdump
#   3. Generate a new netcdf file from the edited text file using ncgen:
#      ncgen vert_coord.txt -o vert_coord.nc
# 
# ==================================================================================================
import os
import glob
import datetime
import numpy as np
import xarray as xr
import subprocess as sp
import hiccup_data_class as hdc
import hiccup_state_adjustment as hsa

verbose = True

# Logical flags for debugging
unpack_nc_files = False
create_map_file = False    # flag for grid and map file creation
remap_data_horz = True    # toggle horizontal remap, variable renaming, and reference pressure
remap_data_vert = True    # toggle vertical remap
do_state_adjust = True    # toggle for all adjustment calculations
create_sst_data = False    # sst/sea ice file creation

output_atm_file_name = 'data/HICCUP_TEST.output.atm.nc'
output_sst_file_name = 'data/HICCUP_TEST.output.sst.nc'

vert_file_name = 'vert_coord_L72.nc'

# topo_file_path = '/project/projectdirs/acme/inputdata/atm/cam/topo/'            # path for NERSC 
# topo_file_name = 'data/USGS-gtopo30_ne30np4pg2_16xdel2.c20200108.nc'
topo_file_name = 'data/USGS-gtopo30_ne30np4_16xdel2-PFC-consistentSGH.nc'

# Create data class instance, which includes xarray file dataset objects
# and variable name dictionaries for mapping between naming conventions.
# This also checks input files for required variables
hiccup_data = hdc.create_hiccup_data(name='ERA5'
                                    ,atm_file='data/HICCUP_TEST.ERA5.atm.low-res.nc'
                                    ,sfc_file='data/HICCUP_TEST.ERA5.sfc.low-res.nc'
                                    # ,atm_file='data/HICCUP_TEST.ERA5.atm.upack.nc'
                                    # ,sfc_file='data/HICCUP_TEST.ERA5.sfc.upack.nc'
                                    ,sstice_name='NOAA'
                                    ,sst_file='data/sst.day.mean.2018.nc'
                                    ,ice_file='data/icec.day.mean.2018.nc'
                                    ,dst_horz_grid='ne30np4'
                                    ,dst_vert_grid='L72'
                                    ,verbose=verbose)

# override the xarray default netcdf format of 
# NETCDF4 to avoid file permission issue
nc_format = 'NETCDF3_64BIT'

# ------------------------------------------------------------------------------
# Make sure files are "unpacked" (may take awhile, so only do it if you need to)
# ------------------------------------------------------------------------------
if unpack_nc_files:

    hiccup_data.unpack_data_files()

# ------------------------------------------------------------------------------
# Create grid and mapping files
# ------------------------------------------------------------------------------
if create_map_file :

    # Create grid description files needed for the mapping file
    hiccup_data.create_src_grid_file()
    hiccup_data.create_dst_grid_file()

    # Create mapping file
    hiccup_data.create_map_file()

# ------------------------------------------------------------------------------
# Horizontally remap the data
# ------------------------------------------------------------------------------
if remap_data_horz :

    # Horizontally regrid the data
    hiccup_data.remap_horizontal(output_file_name=output_atm_file_name)

    # Rename variables to match what the model expects
    hiccup_data.rename_vars(file_name=output_atm_file_name)

    # add P0 variable
    hiccup_data.add_reference_pressure(file_name=output_atm_file_name)

    # Clean up the global attributes of the file
    hiccup_data.clean_global_attributes(file_name=output_atm_file_name)

    # Add time/date information
    ds_data = xr.open_dataset(output_atm_file_name).load()
    hiccup_data.add_time_date_variables( ds_data )
    ds_data.to_netcdf(output_atm_file_name,format=nc_format)
    ds_data.close()

# ------------------------------------------------------------------------------
# Adjust sfc temperature and pressure before vertical interpolation
# ------------------------------------------------------------------------------
if do_state_adjust :

    # Load the file into an xarray dataset
    ds_data = xr.open_dataset(output_atm_file_name).load()
    ds_topo = xr.open_dataset(topo_file_name)

    # Adjust surface temperature to match new surface height
    hsa.adjust_surface_temperature( ds_data, ds_topo )

    # Adjust surface pressure to match new surface height
    hsa.adjust_surface_pressure( ds_data, ds_topo \
                                ,lev_coord_name='plev' \
                                ,pressure_var_name='plev' )

    # Write the adjusted dataset back to the file
    ds_data.to_netcdf(output_atm_file_name,format=nc_format)
    ds_data.close()

# ------------------------------------------------------------------------------
# Vertically remap the data
# ------------------------------------------------------------------------------
if remap_data_vert :

    # Do the vertical interpolation
    hiccup_data.remap_vertical(input_file_name=output_atm_file_name
                              ,output_file_name=output_atm_file_name
                              ,vert_file_name=vert_file_name)

# ------------------------------------------------------------------------------
# Perform final state adjustments on interpolated data and add additional data
# ------------------------------------------------------------------------------
if do_state_adjust :

    # Load the file into an xarray dataset
    ds_data = xr.open_dataset(output_atm_file_name).load()

    # adjust water vapor to eliminate supersaturation
    hsa.remove_supersaturation( ds_data, hybrid_lev=True )

    # adjust cloud water to remove negative values?
    hsa.adjust_cld_wtr( ds_data )

    # adjust cloud fraction to remove values outside of [0,1]
    # hsa.adjust_cloud_fraction( ds_data )

    # adjust surface pressure to retain dry mass of atmosphere
    # hsa.dry_mass_fixer( ds_data )

    # Add extra variable that weren't included in input data - DO WE NEED THIS?
    # hiccup_data.add_extra_data_variables( ds_data )

    # Write the final dataset back to the file
    ds_data.to_netcdf(output_atm_file_name,format=nc_format)
    ds_data.close()

# ------------------------------------------------------------------------------
# Create SST/sea ice file
# ------------------------------------------------------------------------------
if create_sst_data :

    hiccup_data.create_sstice(output_file_name=output_sst_file_name)

# ------------------------------------------------------------------------------
# Print final output file name
# ------------------------------------------------------------------------------

print()
print(f'output_atm_file_name: {output_atm_file_name}')
print(f'output_sst_file_name: {output_sst_file_name}')
print()

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
