#!/usr/bin/env python
# ==================================================================================================
# HICCUP - Hindcast Initial Condition Creation Utility/Processor
# ==================================================================================================
import os, datetime, pandas as pd
import xarray as xr
from hiccup import hiccup_data_class as hdc
from hiccup import hiccup_state_adjustment as hsa
from time import perf_counter
# ------------------------------------------------------------------------------
''' commands to generate model grid file
NE=4
HICCUP_ROOT=/lustre/orion/cli115/proj-shared/hannah6/HICCUP
GenerateCSMesh --alt --res ${NE} --file ${HICCUP_ROOT}/files_grid/exodus_ne${NE}.g
GenerateVolumetricMesh --in ${HICCUP_ROOT}/files_grid/exodus_ne${NE}.g --out ${HICCUP_ROOT}/files_grid/exodus_ne${NE}pg2.g --np 2 --uniform
ConvertMeshToSCRIP --in ${HICCUP_ROOT}/files_grid/exodus_ne${NE}pg2.g --out ${HICCUP_ROOT}/files_grid/scrip_ne${NE}pg2.nc
'''
# ------------------------------------------------------------------------------
''' commands to generate TR bilinear mapping file for online gridding of nudging data
NE_SRC=4 ; NE_DST=30
HICCUP_ROOT=/lustre/orion/cli115/proj-shared/hannah6/HICCUP
SRC_GRID_FILE=${HICCUP_ROOT}/files_grid/scrip_ne${NE_SRC}pg2.nc
DST_GRID_FILE=${HICCUP_ROOT}/files_grid/scrip_ne${NE_DST}pg2.nc
MAP_FILE=${HICCUP_ROOT}/files_map/map_ne${NE_SRC}pg2_to_ne${NE_DST}pg2_trbilin.20240201.nc
ncremap -5 -a trbilin --a2o --src_grd=${SRC_GRID_FILE} --dst_grd=${DST_GRID_FILE} --map_file=${MAP_FILE}
'''
# ------------------------------------------------------------------------------
# Logical flags for controlling what this script will do (comment out to disable)
verbose = True            # Global verbosity flag
unpack_nc_files = True    # unpack data files (convert short to float)
create_map_file = True    # grid and map file creation
remap_data_horz = True    # horz remap, variable renaming
remap_data_vert = True    # vertical remap
combine_files   = True    # combine temporary data files and delete
add_pressure    = True
adjust_fillval  = True
add_case_t0     = True
transpose_dims  = True
remove_ilev     = True
combine_daily   = True
# ------------------------------------------------------------------------------

# bdate_str,edate_str = '2020-01-20 00','2020-01-26 21' # 2024 SCREAM - autocalibration
bdate_str,edate_str = '2020-01-20 00','2020-01-20 21'
# bdate_str,edate_str = '2020-01-21 00','2020-01-26 21'

# Create list of date and time 
beg_date = datetime.datetime.strptime(bdate_str, '%Y-%m-%d %H')
end_date = datetime.datetime.strptime(edate_str, '%Y-%m-%d %H')
datetime_list = pd.date_range(beg_date, end_date, freq='3H')

# datetime_list = [ datetime_list[0] ] # use one file for debugging

# local path for grid, map, and data files
hiccup_root = os.getenv('HOME')+'/HICCUP'
src_data_root = '/lustre/orion/cli115/proj-shared/hannah6/scream_scratch/nudge_data'
dst_data_root = '/lustre/orion/cli115/proj-shared/hannah6/scream_scratch/nudge_data'

# Specify output atmosphere horizontal grid
# dst_horz_grid = 'ne4pg2'
dst_horz_grid = 'ne128pg2'
# dst_horz_grid = 'ne256pg2'
# dst_horz_grid = 'ne512pg2'

# Specify output atmosphere vertical grid
dst_vert_grid,vert_file_name = 'L128',f'{hiccup_root}/files_vert/vert_coord_E3SM_L128.nc'

hdc.target_model = 'EAMXX-nudging'
hdc.verbose_indent = '  '

ref_date = beg_date.strftime("%Y-%m-%d")

# ------------------------------------------------------------------------------
def get_nudge_file(datetime):
    return f'{dst_data_root}/HICCUP.nudging_uv_era5.{datetime}.{dst_horz_grid}.{dst_vert_grid}.nc'
# ------------------------------------------------------------------------------
timer_start = perf_counter() # start timer for entire script
# ------------------------------------------------------------------------------
map_file_created = False
output_file_list = []
for t in datetime_list:

    nudge_datetime = t.strftime("%Y-%m-%d-%H")

    src_atm_file = f'{src_data_root}/ERA5.atm.{nudge_datetime}:00.nc'
    src_sfc_file = f'{src_data_root}/ERA5.sfc.{nudge_datetime}:00.nc'

    output_atm_file_name = get_nudge_file(nudge_datetime)
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
    if verbose:
        print('\n  Input Files')
        print(f'    input atm files: {hiccup_data.atm_file}')
        print(f'    input sfc files: {hiccup_data.sfc_file}')
        print('\n  Output files')
        print(f'    output atm file: {output_atm_file_name}')

    # Get dict of temporary files for each variable
    # dict_timestamp = datetime.datetime.utcnow().strftime('%Y%m%d')
    dict_timestamp = nudge_datetime
    file_dict = hiccup_data.get_multifile_dict(verbose=verbose
                                              ,timestamp=dict_timestamp)
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
    if create_map_file and not map_file_created:
        # Create grid description files
        hiccup_data.create_src_grid_file()
        # hiccup_data.src_grid_file = f'{hiccup_data.grid_dir}/nudge_data/scrip_ERA5_721x1440.nc'
        hiccup_data.create_dst_grid_file()
        # make sure we're using the pg2 grid for output
        hiccup_data.dst_grid_file = hiccup_data.dst_grid_file_pg
        # Create mapping file
        hiccup_data.create_map_file(src_type='FV',dst_type='FV')
        # save the paths and skip this step for further iterations
        map_file_created = True
    else:
        hiccup_data.src_grid_file = f'{hiccup_data.grid_dir}/scrip_ERA5_721x1440.nc'
        hiccup_data.dst_grid_file = f'{hiccup_data.grid_dir}/scrip_{dst_horz_grid}.nc'
        hiccup_data.map_file      = f'{hiccup_data.map_dir}/map_721x1440_to_{dst_horz_grid}.nc'
    # --------------------------------------------------------------------------
    # perform multi-file horizontal remap
    if 'remap_data_horz' not in locals(): remap_data_horz = False
    if remap_data_horz :
        # Horizontally regrid the data
        hiccup_data.remap_horizontal_multifile(file_dict)
        # Rename variables to match what the model expects
        hiccup_data.rename_vars_multifile(file_dict=file_dict)
        # Add time/date information
        hiccup_data.add_time_date_variables_multifile(file_dict=file_dict
                                                     ,ref_date=ref_date)
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
        # cmd = f'ncap2 -O -s \'p_mid[time,lev,ncol]=100000.0*hyam+PS*hybm\' {output_atm_file_name} {output_atm_file_name}'
        cmd = f'ncap2 -O -s \'p_mid[time,ncol,lev]=100000.0*hyam+PS*hybm\' {output_atm_file_name} {output_atm_file_name}'
        hdc.run_cmd(cmd,verbose=True,prepend_line=False,shell=True)
        # add _FillValue attribute
        cmd = f'ncatted -O -a _FillValue,p_mid,o,d,1e30 {output_atm_file_name}  {output_atm_file_name}'
        hdc.run_cmd(cmd,verbose=True,prepend_line=False,shell=True)
    # --------------------------------------------------------------------------
    # Use NCO to add p_mid field and case_t0 attribute
    if 'adjust_fillval' not in locals(): adjust_fillval = False
    if adjust_fillval:
        ds = xr.open_dataset(output_atm_file_name)
        print()
        fval_var_list = []
        mval_var_list = []
        xval_var_list = []
        for v in (list(ds.coords)+list(ds.data_vars)):
            if '_FillValue'    in ds[v].encoding.keys(): fval_var_list.append(v)
            if 'missing_value' in ds[v].encoding.keys(): mval_var_list.append(v)
            if 'eulaVlliF_'    in ds[v].attrs.keys():    xval_var_list.append(v) # not sure why this attribute gets added
            if 'eulaVlliF_'    in ds[v].encoding.keys(): xval_var_list.append(v) # not sure why this attribute gets added
        # update _FillValue to ensure consistency and a large positive value required by SCREAM
        if fval_var_list != []:
            for v in fval_var_list:
                cmd = f'ncatted -h -O -a _FillValue,{v},o,d,1e30 {output_atm_file_name}  {output_atm_file_name}'
                hdc.run_cmd(cmd,verbose=True,prepend_line=False,shell=True)
        # delete the redundant missing_value attribute from all variables
        if mval_var_list != []:
            for v in mval_var_list:
                cmd = f'ncatted -h -O -a missing_value,{v},d,, {output_atm_file_name}  {output_atm_file_name}'
                hdc.run_cmd(cmd,verbose=True,prepend_line=False,shell=True)
        # fix mangled fill value attribute
        if xval_var_list != []:
            for v in xval_var_list:
                cmd = f'ncatted -h -O -a eulaVlliF_,{v},d,, {output_atm_file_name}  {output_atm_file_name}'
                hdc.run_cmd(cmd,verbose=True,prepend_line=False,shell=True)
    # --------------------------------------------------------------------------
    # Use NCO to add global attribute case_t0
    if 'add_case_t0' not in locals(): add_case_t0 = False
    if add_case_t0:
        # This is temporary - won't be needed after PR to overhaul time interpolation in SCREAM
        case_t0 = f'{ref_date}-00000'
        cmd = f'ncatted -O -a case_t0,global,c,c,\'{case_t0}\' {output_atm_file_name}'
        hdc.run_cmd(cmd,verbose=True,prepend_line=False,shell=True)
    # --------------------------------------------------------------------------
    if 'transpose_dims' not in locals(): transpose_dims = False
    if transpose_dims:
        print(hdc.verbose_indent+'\nTransposing data dimensions...')
        cmd = f'ncpdq -O --rdr=time,ncol,lev {output_atm_file_name} {output_atm_file_name}'
        hdc.run_cmd(cmd,verbose=True,prepend_line=False,shell=True)
    # --------------------------------------------------------------------------
    if 'remove_ilev' not in locals(): remove_ilev = False
    if remove_ilev:
        print(hdc.verbose_indent+'\nRemoving variables with ilev dimension...')
        cmd = f'ncks -C -O -x -v ilev,hyai,hybi {output_atm_file_name} {output_atm_file_name}'
        hdc.run_cmd(cmd,verbose=True,prepend_line=False,shell=True)
# ------------------------------------------------------------------------------
# Print final output file names
print('\nhourly output files:')
for ofile in output_file_list: print(f'  {ofile}')
print()
# ------------------------------------------------------------------------------
if 'combine_daily' not in locals(): combine_daily = False
if combine_daily:
    print('\nCombining sub-daily files into daily files...')

    # build list of daily files and identify old daily files to delete
    dy_file_list = []
    dy_datestr_list = []
    dy_datetime_list = []
    for t in datetime_list:
        dy_datetime = t.strftime("%Y-%m-%d")
        # add to list of daily dates (no subdaily dates)
        if dy_datetime not in dy_datestr_list: 
            dy_datestr_list.append(dy_datetime)
            dy_datetime_list.append(t)
        # check if the daily file already exists
        dy_file_name = get_nudge_file(dy_datetime)
        if dy_file_name not in dy_file_list:
            if os.path.isfile(dy_file_name):
                dy_file_list.append(dy_file_name)

    # delete old daily files
    if dy_file_list!=[]:
        print(hdc.verbose_indent+'removing old daily files...')
        for dy_file_name in dy_file_list:
            hdc.run_cmd(f'rm {dy_file_name} ',verbose=True,prefix=hdc.verbose_indent)
        print()

    # combine hourly files for each day
    dy_file_list = []
    for d in dy_datetime_list:
        dy_date = d.strftime("%Y-%m-%d")
        dy_file_name = get_nudge_file(dy_date)
        dy_file_list.append(dy_file_name)
        # build list of sub-daily files for this date
        hr_file_list = []
        for h in datetime_list:
            hr_date = h.strftime("%Y-%m-%d")
            if hr_date==dy_date:
                hr_datetime = h.strftime("%Y-%m-%d-%H")
                hr_file_list.append(get_nudge_file(hr_datetime))

        hr_file_list_str = ' '.join(hr_file_list)
        cmd = f'ncrcat -5 -O -h {hr_file_list_str} {dy_file_name} '
        hdc.run_cmd(cmd,verbose=True,shell=True,prefix=hdc.verbose_indent)

    print('\noutput files:')
    for dfile in dy_file_list: print(f'  {dfile}')
    print()
# ------------------------------------------------------------------------------

hdc.print_timer(timer_start,caller=f'Total',print_msg=False); print()

# ------------------------------------------------------------------------------
