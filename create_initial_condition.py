#!/usr/bin/env python
# ==================================================================================================
# HICCUP - Hindcast Initial Condition Creation Utility/Processor
# This tool automates the creation of atmospheric initial condition files for 
# E3SM using user supplied file for atmospheric and sea surface conditions.
# ==================================================================================================
import os
import glob
import datetime
import numpy as np
import xarray as xr
import subprocess as sp
import hiccup_data_class as hdc
import hiccup_state_adjustment as hsa
from optparse import OptionParser
# ------------------------------------------------------------------------------
# Logical flags for controlling what this script will do
verbose = True            # Global verbosity flag
unpack_nc_files = False    # unpack data files (convert short to float)
create_map_file = False    # grid and map file creation
remap_data_horz = True    # horizontal remap and variable renaming
do_state_adjst1 = True    # post horizontal interpolation adjustments
remap_data_vert = True    # vertical remap
do_state_adjst2 = True    # post vertical interpolation adjustments
create_sst_data = False    # sst/sea ice file creation
# ------------------------------------------------------------------------------
# Parse the command line options
parser = OptionParser()
parser.add_option('--hgrid',dest='horz_grid',default=None,help='Sets the output horizontal grid')
parser.add_option('--vgrid',dest='vert_grid',default=None,help='Sets the output vertical grid')
(opts, args) = parser.parse_args()
# ------------------------------------------------------------------------------

# Specify output atmosphere horizontal grid ( ne30np4 / ne120np4 / ne1024np4 )
dst_horz_grid = opts.horz_grid if opts.horz_grid is not None else 'ne30np4'

# Specify output atmosphere vertical grid
dst_vert_grid = opts.vert_grid if opts.vert_grid is not None else 'L72'
vert_file_name = f'vert_coord_{dst_vert_grid}.nc'

# Specify the output file names
# data_root = '/global/cscratch1/sd/whannah/HICCUP/data/'
data_root = '/gpfs/alpine/scratch/hannah6/cli115/HICCUP/data/'
output_atm_file_name = f'{data_root}HICCUP_TEST.output.atm.{dst_horz_grid}.{dst_vert_grid}.nc'
output_sst_file_name = f'{data_root}HICCUP_TEST.output.sst.{dst_horz_grid}.{dst_vert_grid}.nc'

# set topo file
# topo_file_path = data_root
# topo_file_path = '/project/projectdirs/acme/inputdata/atm/cam/topo/'  # NERSC
topo_file_path = '/gpfs/alpine/world-shared/csc190/e3sm/cesm/inputdata/atm/cam/topo/' # OLCF
if dst_horz_grid=='ne1024np4': topo_file_name = f'{topo_file_path}USGS-gtopo30_ne1024np4_16xconsistentSGH_20190528.nc'
if dst_horz_grid=='ne120np4' : topo_file_name = f'{topo_file_path}USGS-gtopo30_ne120np4_16xdel2-PFC-consistentSGH.nc'
if dst_horz_grid=='ne30np4'  : topo_file_name = f'{topo_file_path}USGS-gtopo30_ne30np4_16xdel2-PFC-consistentSGH.nc'

# Create data class instance, which includes xarray file dataset objects
# and variable name dictionaries for mapping between naming conventions.
# This also checks input files for required variables
hiccup_data = hdc.create_hiccup_data(name='ERA5'
                                    # ,atm_file='data/HICCUP_TEST.ERA5.atm.low-res.nc'
                                    # ,sfc_file='data/HICCUP_TEST.ERA5.sfc.low-res.nc'
                                    ,atm_file=f'{data_root}ERA5.atm.2018-01-01.nc'
                                    ,sfc_file=f'{data_root}ERA5.sfc.2018-01-01.nc'
                                    # ,sstice_name='NOAA'
                                    # ,sst_file=f'{data_root}sst.day.mean.2018.nc'
                                    # ,ice_file=f'{data_root}icec.day.mean.2018.nc'
                                    # ,sstice_name='ERA5'
                                    # ,sstice_combined_file='data_scratch/HICCUP_TEST.ERA5.sfc.upack.nc'
                                    # ,dst_horz_grid='ne30np4'
                                    ,dst_horz_grid=dst_horz_grid
                                    ,dst_vert_grid='L72'
                                    ,output_dir=data_root
                                    ,grid_dir=data_root
                                    ,map_dir=data_root
                                    ,tmp_dir=data_root
                                    ,verbose=verbose)

# override the xarray default netcdf format of 
# NETCDF4 to avoid file permission issue
nc_format = 'NETCDF3_64BIT'
# nc_format = 'NETCDF4_CLASSIC'
# nc_format = 'NETCDF4'

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
    ds_data = xr.open_dataset(output_atm_file_name)#.load()
    hiccup_data.add_time_date_variables( ds_data )
    ds_data.to_netcdf(output_atm_file_name,format=nc_format,mode='a')
    # ds_data.to_netcdf(output_atm_file_name,format=nc_format,mode='w')
    ds_data.close()

# ------------------------------------------------------------------------------
# Adjust sfc temperature and pressure before vertical interpolation
# ------------------------------------------------------------------------------
if do_state_adjst1 :

    # Load the file into an xarray dataset
    ds_data = xr.open_dataset(output_atm_file_name)
    ds_topo = xr.open_dataset(topo_file_name)

    # Adjust surface temperature to match new surface height
    hsa.adjust_surface_temperature( ds_data, ds_topo, verbose=verbose )

    # Adjust surface pressure to match new surface height
    hsa.adjust_surface_pressure( ds_data, ds_topo \
                                ,lev_coord_name='plev' \
                                ,pressure_var_name='plev'
                                ,verbose=verbose )

    # Write the adjusted dataset back to the file
    ds_data.to_netcdf(output_atm_file_name,format=nc_format,mode='a')
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
if do_state_adjst2 :

    # Load the file into an xarray dataset
    ds_data = xr.open_dataset(output_atm_file_name)

    # adjust water vapor to eliminate supersaturation
    hsa.remove_supersaturation( ds_data, hybrid_lev=True, verbose=verbose )

    # adjust cloud water to remove negative values?
    hsa.adjust_cld_wtr( ds_data, verbose=verbose )

    # adjust cloud fraction to remove values outside of [0,1] - DO WE NEED THIS?
    # hsa.adjust_cloud_fraction( ds_data, verbose=verbose )

    # adjust surface pressure to retain dry mass of atmosphere - NOT TESTED
    # hsa.dry_mass_fixer( ds_data )

    # Add extra variable that weren't included in input data - DO WE NEED THIS?
    # hiccup_data.add_extra_data_variables( ds_data )

    # Write the final dataset back to the file
    ds_data.to_netcdf(output_atm_file_name,format=nc_format,mode='a')
    ds_data.close()

# ------------------------------------------------------------------------------
# Create SST/sea ice file
# ------------------------------------------------------------------------------
if create_sst_data :

    # create grid and mapping files
    overwrite = False
    hiccup_data.sstice_create_src_grid_file(force_overwrite=overwrite)
    hiccup_data.sstice_create_dst_grid_file(force_overwrite=overwrite)
    hiccup_data.sstice_create_map_file(force_overwrite=overwrite)

    # Remap the sst/ice data after time slicing and combining (if necessary)
    hiccup_data.sstice_slice_and_remap(output_file_name=output_sst_file_name,
                                       time_slice_method='initial',
                                       atm_file=output_atm_file_name)

    # Rename the variables and remove unnecessary variables and attributes
    hiccup_data.sstice_rename_vars(output_file_name=output_sst_file_name)

    # Adjust final SST/ice data to fill in missing values and limit ice fraction
    hiccup_data.sstice_adjustments(output_file_name=output_sst_file_name)

# ------------------------------------------------------------------------------
# Print final output file name
# ------------------------------------------------------------------------------

print()
print(f'output_atm_file_name: {output_atm_file_name}')
print(f'output_sst_file_name: {output_sst_file_name}')
print()

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
