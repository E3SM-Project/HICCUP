#!/usr/bin/env python
# ==================================================================================================
# HICCUP - Hindcast Initial Condition Creation Utility/Processor
# ==================================================================================================
import os, datetime, pandas as pd
import xarray as xr
from hiccup import hiccup_data_class as hdc
from hiccup import hiccup_state_adjustment as hsa
# ------------------------------------------------------------------------------
# Logical flags for controlling what this script will do (comment out to disable)
verbose = True            # Global verbosity flag
# unpack_nc_files = True    # unpack data files (convert short to float)
# create_map_file = True    # grid and map file creation
# remap_data_horz = True    # horz remap, variable renaming
# remap_data_vert = True    # vertical remap
# combine_files   = True    # combine temporary data files and delete
# add_pressure    = True
combine_daily   = True
# ------------------------------------------------------------------------------

# Create list of date and time 
beg_date = datetime.datetime.strptime('2006-08-09 00', '%Y-%m-%d %H')
end_date = datetime.datetime.strptime('2006-08-13 21', '%Y-%m-%d %H')
datetime_list = pd.date_range(beg_date, end_date, freq='3H')

# local path for grid, map, and data files
hiccup_root = os.getenv('HOME')+'/HICCUP'
src_data_root = '/global/cfs/cdirs/m2637/whannah/nudge_data'
dst_data_root = '/global/cfs/cdirs/m2637/whannah/nudge_data'

# Specify output atmosphere horizontal grid
dst_horz_grid = 'ne0np4-saomai-128x8'

# Specify output atmosphere vertical grid
dst_vert_grid,vert_file_name = 'L128',f'{hiccup_root}/files_vert/vert_coord_E3SM_L128.nc'

dst_grid_file = '/global/cfs/cdirs/m2637/jsgoodni/Saomai_2006_ne128x8_lon130E_lat25Npg2.scrip.nc'
# map_file = '/global/cfs/cdirs/m2637/whannah/'

hdc.target_model = 'EAMXX-nudging'
hdc.verbose_indent = '  '

# ------------------------------------------------------------------------------
output_file_list = []
for t in datetime_list:

    nudge_datetime = t.strftime("%Y-%m-%d-%H")

    src_atm_file = f'{src_data_root}/ERA5.atm.{nudge_datetime}:00.nc'
    src_sfc_file = f'{src_data_root}/ERA5.sfc.{nudge_datetime}:00.nc'

    output_atm_file_name = f'{dst_data_root}/HICCUP.nudging_uv_era5.{nudge_datetime}.{dst_horz_grid}.{dst_vert_grid}.nc'
    output_file_list.append(output_atm_file_name)

    # Create data class instance
    hiccup_data = hdc.create_hiccup_data(name='ERA5'
                                        ,dst_horz_grid=dst_horz_grid
                                        ,dst_vert_grid=dst_vert_grid
                                        ,atm_file=src_atm_file
                                        ,sfc_file=src_sfc_file
                                        ,topo_file=None
                                        ,output_dir=dst_data_root
                                        ,grid_dir=dst_data_root
                                        ,map_dir=dst_data_root
                                        ,tmp_dir=dst_data_root
                                        ,verbose=verbose
                                        ,check_input_files=False) # only need U/V
    # Print some informative stuff
    print('\n  Input Files')
    print(f'    input atm files: {hiccup_data.atm_file}')
    print(f'    input sfc files: {hiccup_data.sfc_file}')
    print('\n  Output files')
    print(f'    output atm file: {output_atm_file_name}')
    # Get dict of temporary files for each variable
    file_dict = hiccup_data.get_multifile_dict()
    # # Specify mapping file
    # hiccup_data.map_file = map_file
    # --------------------------------------------------------------------------
    # Make sure files are "unpacked" (may take awhile, so only do it if you need to)
    if 'unpack_nc_files' not in locals(): unpack_nc_files = False
    if unpack_nc_files:
        hiccup_data.unpack_data_files()
    # ------------------------------------------------------------------------------
    # Create grid and mapping files
    if 'create_map_file' not in locals(): create_map_file = False
    if create_map_file :
        # Create destination grid description file
        hiccup_data.create_src_grid_file()
        # specify destination grid file
        hiccup_data.dst_grid_file = dst_grid_file
        # Create mapping file
        hiccup_data.create_map_file(src_type='FV',dst_type='FV')
    # --------------------------------------------------------------------------
    # perform multi-file horizontal remap
    if 'remap_data_horz' not in locals(): remap_data_horz = False
    if remap_data_horz :
        # Horizontally regrid the data
        hiccup_data.remap_horizontal_multifile(file_dict)
        # Rename variables to match what the model expects
        hiccup_data.rename_vars_multifile(file_dict=file_dict)
        # Add time/date information
        hiccup_data.add_time_date_variables_multifile(file_dict=file_dict)
    # --------------------------------------------------------------------------
    # Vertically remap the data
    if 'remap_data_vert' not in locals(): remap_data_vert = False
    if remap_data_vert :
        hiccup_data.remap_vertical_multifile(file_dict=file_dict
                                            ,vert_file_name=vert_file_name)
    # --------------------------------------------------------------------------
    # Combine files
    if 'combine_files' not in locals(): combine_files = False
    if combine_files :
        # Combine and delete temporary files
        hiccup_data.combine_files(file_dict=file_dict
                                 ,delete_files=True
                                 ,output_file_name=output_atm_file_name)
        # Clean up the global attributes of the file
        hiccup_data.clean_global_attributes(file_name=output_atm_file_name)
    # --------------------------------------------------------------------------
    # Use NCO to add p_mid field and case_t0 attribute
    if 'add_pressure' not in locals(): add_pressure = False
    if add_pressure:

        print(hdc.verbose_indent+'\nAdding pressure field to nudging data...')

        cmd = f'ncap2 -O -s \'p_mid[time,lev,ncol]=100000.0*hyam+PS*hybm\' {output_atm_file_name} {output_atm_file_name}'
        hdc.run_cmd(cmd,verbose=True,prepend_line=False,shell=True)

        datetime_att = t.strftime("%Y-%m-%d-%H000")
        cmd = f'ncatted -O -a case_t0,global,c,c,\'{datetime_att}\' {output_atm_file_name}'
        hdc.run_cmd(cmd,verbose=True,prepend_line=False,shell=True)
    # --------------------------------------------------------------------------
    # # Use xarray to add number concentration fields
    # if 'add_number' not in locals(): add_number = False
    # if add_number:
    #     with xr.open_dataset(ofile) as ds:
    #         ds.load()
    #         ds['nc'] = xr.zeros_like( ds['qv'] )
    #         ds['ni'] = xr.zeros_like( ds['qv'] )
    #         ds['nc'].attrs['long_name'] = 'Specific ice number concentration'
    #         ds['ni'].attrs['long_name'] = 'Specific liq number concentration'
    #         ds['nc'].attrs['units'] = 'kg**-1'
    #         ds['ni'].attrs['units'] = 'kg**-1'
    ### NOTE - we don't need this as long as the defaults are set properly for 
    ### each grid in => components/eamxx/cime_config/namelist_defaults_scream.xml
# ------------------------------------------------------------------------------
# Print final output file names
print('\noutput files:')
for ofile in output_file_list: print(f'  {ofile}')
print()
# ------------------------------------------------------------------------------
if 'combine_daily' not in locals(): combine_daily = False
if combine_daily:
    print('\nCombining sub-daily files into daily files...')

    for t in datetime_list:

        hr_datetime = t.strftime("%Y-%m-%d-%H")
        dy_datetime = t.strftime("%Y-%m-%d")

        hr_file_name = f'{dst_data_root}/HICCUP.nudging_uv_era5.{hr_datetime}.{dst_horz_grid}.{dst_vert_grid}.nc'
        dy_file_name = f'{dst_data_root}/HICCUP.nudging_uv_era5.{dy_datetime}.{dst_horz_grid}.{dst_vert_grid}.nc'

        cmd = f'ncks -A --no_tmp_fl {hr_file_name} {dy_file_name} '
        hdc.run_cmd(cmd,verbose=True,prepend_line=False,shell=True)

# ------------------------------------------------------------------------------
