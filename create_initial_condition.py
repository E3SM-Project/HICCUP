#!/usr/bin/env python
# ===================================================================================================
# HICCUP - Hindcast Initial Condition Creation Utility/Processor
# This tool automates the creation of atmospheric initial condition files
# for E3SM using a user supplied reanalysis file, such as EAR5 data.
# Requires NCO, TempestRemap, and xarray.
# ===================================================================================================
import os
import subprocess as sp
import glob
import datetime
import xarray as xr
import hiccup_data_class as hdc
import hiccup_state_adjustment

# Options for output state adjustment
adjust_ts = False    # Adjust surface temperature to match new surface height
adjust_ps = False    # Adjust surface pressure to match new surface height
adjust_mass = False  # adjust surface pressure to retain dry mass of atmosphere
adjust_qv = False    # adjust qv to eliminate supersaturation
adjust_cw = False    # adjust cloud water to remove negative values
adjust_cf = False    # adjust cloud fraction to remove values outside of [0,1]

verbose = True

output_file_name = 'HICCUP_TEST.output.nc'

# Create data class instance, which includes xarray file dataset objects
# and variable name dictionaries for mapping between naming conventions
hiccup_data = hdc.create_hiccup_data(name='ERA5', atm_file='HICCUP_TEST.ERA5.atm.remap.nc',
                                     sfc_file='HICCUP_TEST.ERA5.sfc.remap.nc',
                                     dst_horz_grid='ne30pg2', dst_vert_grid='L72')

# Check input files for for required variables
hiccup_data.check_file_vars()   # cjones note: let's fold this into the create_hiccup_data() call

# -------------------------------------------------------------------------------
# Create grid and mapping files
# -------------------------------------------------------------------------------

# Create grid description files needed for the mapping file
hiccup_data.create_src_grid_file(verbose=verbose)
hiccup_data.create_dst_grid_file(verbose=verbose)

# Create mapping file
hiccup_data.create_map_file(verbose=verbose)

# -------------------------------------------------------------------------------
# Remap the data
# -------------------------------------------------------------------------------

# Horizontally regrid the data
hiccup_data.remap_horizontal(output_file_name=output_file_name, verbose=verbose)

# Rename variables to match what the model expects
hiccup_data.rename_vars(output_file_name, verbose=verbose)

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# add P0 variable
sp.call(f"ncap2 -O -s 'P0=100000.' {output_file_name}".split())
sp.call(f'ncatted -O -a long_name,P0,a,c,"reference pressure" {output_file_name}'.split())
sp.call(f'ncatted -O -a units,P0,a,c,"Pa" {output_file_name}'.split())

# Rename pressure variable and change type to double (needed for vertical remap)
new_lev_name = 'plev'
sp.call(f'ncrename -d level,{new_lev_name} -v level,{new_lev_name} {output_file_name}'.split())
sp.call(f"ncap2 -O -s '{new_lev_name}={new_lev_name}.convert(NC_DOUBLE)' {output_file_name} {output_file_name}".split())

# Remove lat/lon vertices variables
sp.call(f'ncks -C -O  -x -v lat_vertices,lon_vertices {output_file_name} {output_file_name}'.split())

# Clean up up the global file attributes
if verbose:
    print('\nCleaning up excessive global attributes...')
global_att_list = ['history_of_appended_files', 'nco_openmp_thread_number', 'input_file',
                   'map_file', 'remap_version', 'remap_hostname', 'remap_command', 'remap_script',
                   'NCO', 'history']
for att in global_att_list:
    cmd = f'ncatted -O -a {att},global,d,, {output_file_name} {output_file_name}'
    if verbose:
        print(f'  {cmd}')
    sp.call(cmd.split())

# cmd = f'ncatted -O -a history,global,m,, {output_file_name} {output_file_name}'
# print(f'\n  {cmd}\n')
# sp.call(cmd, shell=True)


# -------------------------------------------------------------------------------
# Vertically remap the data
# -------------------------------------------------------------------------------

# To create the vertical coordinate file it is easiest to
# extract from a pre-existing model data file as follows:
# 1. ncdump -v P0,hyam,hybm,hyai,hybi,lev,ilev <history_file> > vert_coord.txt
# 2. < edit the file to remove extra header info - but keep the general CDL format >
# 3. ncgen vert_coord.txt -o vert_coord.nc

tmp_vert_file_name = f'vert_coord_{hiccup_data.dst_vert_grid}.nc'

# Define list of variables that will be vertical remapped
# vert_remap_var_list = 'T,Q,U,V,P0,PS,PHIS,CLDLIQ,CLDICE,O3,date,datesec,hyam,hybm,hyai,hybi,lev,ilev'
# vert_remap_var_list = 'T,Q,U,V,P0,PS,CLDLIQ,CLDICE,O3'
vert_remap_var_list = 'T'

vert_output_file = output_file_name.replace('.nc', f'.{hiccup_data.dst_vert_grid}.nc')

# Perform the vertical remapping
cmd = f'ncremap --vrt_fl={tmp_vert_file_name} -v {vert_remap_var_list} {output_file_name} {vert_output_file}'
if verbose:
    print(f'\n  {cmd}\n')
sp.call(cmd.split())

# Delete the temporary files
# sp.call(f'rm {tmp_vert_file_text} {tmp_vert_file_name}', shell=True)

print(f'\n{vert_output_file}\n')
exit()

# -------------------------------------------------------------------------------
# Perform state adjustments on interpolated data
# -------------------------------------------------------------------------------

if any([adjust_ts, adjust_ps, adjust_mass, adjust_qv, adjust_cw, adjust_cf]):

    # Load the file into an xarray dataset
    ds = xr.open_dataset(output_file_name)

    # Adjust surface temperature to match new surface height
    if adjust_ts:
        state_adjustment.adjust_surface_temperature(
            ncol, phis_old, ts_old, phis_new, ts_new)

    # Adjust surface pressure to match new surface height
    if adjust_ps:
        state_adjustment.adjust_surface_pressure(plev, ncol, temperature_mid,
                                                 pressure_mid, pressure_int,
                                                 phis_old, ps_old, phis_new, ps_new)

    # adjust surface pressure to retain dry mass of atmosphere?
    # if adjust_mass :
    #   state_adjustment.dry_mass_fixer( ncol, plev, hyai, hybi, wgt, qv, mass_ref, ps_in, ps_out )

    # adjust qv to eliminate supersaturation?
    if adjust_qv:
        state_adjustment.remove_supersaturation()

    # adjust cloud water to remove negative values?
    if adjust_cw:
        cld_liq = cld_liq.where(cld_liq >= 0, other=0.)
        cld_ice = cld_ice.where(cld_ice >= 0, other=0.)

    # adjust cloud fraction to remove values outside of [0,1]
    if adjust_cf:
        state_adjustment.adjust_cloud_fraction()

    # Write the adjusted dataset back to the file
    ds.to_netcdf(output_file_name)

# -------------------------------------------------------------------------------
# Print final output file name
# -------------------------------------------------------------------------------

print(f'\noutput_file_name: {output_file_name}\n')

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------
