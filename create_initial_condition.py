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
import xarray as xr
import subprocess as sp
import hiccup_data_class as hdc
import hiccup_state_adjustment as hsa

verbose = True

recreate_map_file = False
remap_data_horz   = False

# Adjustment options
adjust_sfc_temp = True     # Adjust surface temperature to match new surface height
adjust_sfc_pres = True     # Adjust surface pressure to match new surface height
adjust_supersat = False     # adjust qv to eliminate supersaturation
adjust_glb_mass = False     # adjust surface pressure to retain dry mass of atmosphere
adjust_cld_wtr  = False     # adjust cloud water to remove negative values
adjust_cld_frac = False     # adjust cloud fraction to remove values outside of [0,1]

output_file_name = 'HICCUP_TEST.output.nc'

# topo_file_path = '/project/projectdirs/acme/inputdata/atm/cam/topo/'
topo_file_name = 'USGS-gtopo30_ne30np4pg2_16xdel2.c20200108.nc'

# Create data class instance, which includes xarray file dataset objects
# and variable name dictionaries for mapping between naming conventions
hiccup_data = hdc.create_hiccup_data(name='ERA5'
                                    ,atm_file='HICCUP_TEST.ERA5.atm.remap.nc'
                                    ,sfc_file='HICCUP_TEST.ERA5.sfc.remap.nc'
                                    ,dst_horz_grid='ne30pg2'
                                    ,dst_vert_grid='L72'
                                    ,verbose=verbose)

# Check input files for for required variables
hiccup_data.check_file_vars()   # cjones note: let's fold this into the create_hiccup_data() call

# ------------------------------------------------------------------------------
# Create grid and mapping files
# ------------------------------------------------------------------------------
if recreate_map_file :

    # Create grid description files needed for the mapping file
    hiccup_data.create_src_grid_file()
    hiccup_data.create_dst_grid_file()

    # Create mapping file
    hiccup_data.create_map_file()

# ------------------------------------------------------------------------------
# Horizontally remap the data
# ------------------------------------------------------------------------------

if False :

  # Horizontally regrid the data
  hiccup_data.remap_horizontal(output_file_name=output_file_name)

  # Rename variables to match what the model expects
  hiccup_data.rename_vars(file_name=output_file_name)

  # add P0 variable
  hiccup_data.add_reference_pressure(file_name=output_file_name)

  # Clean up the global attributes of the file
  hiccup_data.clean_global_attributes(file_name=output_file_name)

# exit(f'\n{output_file_name}\n')

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if any([adjust_sfc_temp, adjust_sfc_pres]):

    # Load the file into an xarray dataset
    ds_data = xr.open_dataset(output_file_name)
    ds_topo = xr.open_dataset(topo_file_name)

    # Adjust surface temperature to match new surface height
    if adjust_sfc_temp : hsa.adjust_surface_temperature( ds_data, ds_topo )

    # Adjust surface pressure to match new surface height
    if adjust_sfc_pres : hsa.adjust_surface_pressure( ds_data, ds_topo, lev_coord_name='plev', debug=True )

    # Write the adjusted dataset back to the file
    # ds.to_netcdf(output_file_name)

exit(f'\n{output_file_name}\n')

# ------------------------------------------------------------------------------
# Vertically remap the data
# ------------------------------------------------------------------------------

# Define list of variables that will be vertical remapped
# vert_remap_var_list = 'T,Q,U,V,P0,PS,PHIS,CLDLIQ,CLDICE,O3,date,datesec,hyam,hybm,hyai,hybi,lev,ilev'
# vert_remap_var_list = 'T,Q,U,V,P0,PS,CLDLIQ,CLDICE,O3'
# vert_remap_var_list = 'T'

vert_tmp_file_name = output_file_name.replace('.nc',f'.{hiccup_data.dst_vert_grid}.nc')

hiccup_data.remap_vertical(input_file_name=output_file_name
                          ,output_file_name=vert_tmp_file_name
                          ,vert_file_name=f'vert_coord_{hiccup_data.dst_vert_grid}.nc')

sp.call(f'mv {vert_tmp_file_name} {output_file_name} '.split())

# ------------------------------------------------------------------------------
# Perform state adjustments on interpolated data
# ------------------------------------------------------------------------------

if any([adjust_glb_mass, adjust_supersat, adjust_cld_wtr, adjust_cld_frac]):

    # Load the file into an xarray dataset
    ds = xr.open_dataset(output_file_name)

    # adjust surface pressure to retain dry mass of atmosphere
    # if adjust_glb_mass :
    #   hsa.dry_mass_fixer( ncol, plev, hyai, hybi, wgt, qv, mass_ref, ps_in, ps_out )

    # adjust water vapor to eliminate supersaturation
    if adjust_supersat:
        hsa.remove_supersaturation()

    # adjust cloud water to remove negative values?
    if adjust_cld_wtr:
        cld_liq = cld_liq.where(cld_liq >= 0, other=0.)
        cld_ice = cld_ice.where(cld_ice >= 0, other=0.)

    # adjust cloud fraction to remove values outside of [0,1]
    if adjust_cld_frac:
        hsa.adjust_cloud_fraction()

    # Write the adjusted dataset back to the file
    ds.to_netcdf(output_file_name)

# ------------------------------------------------------------------------------
# Print final output file name
# ------------------------------------------------------------------------------

print(f'\noutput_file_name: {output_file_name}\n')

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
