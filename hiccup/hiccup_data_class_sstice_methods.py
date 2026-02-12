import os
import numpy as np
import xarray as xr
import pandas as pd
from time import perf_counter
from hiccup.hiccup_data_class import xarray_sst_nc_format
from hiccup.hiccup_data_class import hdr_pad
from hiccup.hiccup_data_class import ncremap_file_fmt
from hiccup.hiccup_utilities import check_dependency
from hiccup.hiccup_utilities import run_cmd
# ------------------------------------------------------------------------------
# HICCUP SST and sea ice methods
# ------------------------------------------------------------------------------
def get_sst_file(self):
    """
    Return file name that contains input SST data
    """
    sst_file = self.sstice_combined_file
    if sst_file is None : sst_file = self.sst_file
    if sst_file is None : raise ValueError('No valid sst data file found!')
    return sst_file
# ------------------------------------------------------------------------------
def sstice_create_src_grid_file(self,diagnose_grid=True,nlat=None,nlon=None,
                                force_overwrite=False,verbose=None):
    """
    Create a source grid file to use for remapping the SST and sea ice data.
    The SST and ice data are assumed to exist on the same grid.
    """
    if self.do_timers: timer_start = perf_counter()
    if verbose is None : verbose = self.verbose

    if diagnose_grid:
        # Load the SST file as xarray datasets to read the grid dimensions
        sst_file = self.sstice_combined_file
        if sst_file is None : sst_file = self.sst_file
        if sst_file is None : raise ValueError('No valid sst data file found!')

        ds_sst = xr.open_dataset(sst_file)

        # Determine input grid 
        if self.sstice_name=='NOAA': lat_name, lon_name = 'lat', 'lon'
        if self.sstice_name=='ERA5': lat_name, lon_name = 'latitude', 'longitude'
        self.sstice_nlat_src = len( ds_sst[lat_name].values )
        self.sstice_nlon_src = len( ds_sst[lon_name].values )

        # Close the dataset
        ds_sst.close()
    else:
        # Use supplied grid dimensions
        if nlat is None : raise ValueError('nlat can not be None if diagnose_grid is False!')
        if nlon is None : raise ValueError('nlon can not be None if diagnose_grid is False!')
        self.sstice_nlat_src = nlat
        self.sstice_nlon_src = nlon

    src_grid = f'{self.sstice_nlat_src}x{self.sstice_nlon_src}'

    # Define the source grid file to be created
    self.sstice_src_grid_file = f'{self.grid_dir}/scrip_{src_grid}.nc'

    # Create the source grid file
    if force_overwrite or not os.path.isfile(self.sstice_src_grid_file) :
        if verbose : print(self.verbose_indent+f'\nCreating source grid file for SST and sea ice data...')
        cmd  = f'ncremap'
        cmd += f' --tmp_dir={self.tmp_dir}'
        cmd += f' -G ttl=\'Equi-Angular grid {src_grid}\'' 
        cmd += f'#latlon={self.sstice_nlat_src},{self.sstice_nlon_src}'
        cmd +=  '#lat_typ=uni'
        cmd +=  '#lon_typ=grn_ctr'
        if self.sstice_name=='NOAA': cmd += '#lat_drc=s2n'
        if self.sstice_name=='ERA5': cmd += '#lat_drc=n2s'
        cmd += f' -g {self.sstice_src_grid_file} '
        run_cmd(cmd,verbose,shell=True,prepend_line=False)

    if self.do_timers: self.print_timer(timer_start)
    return
# ------------------------------------------------------------------------------
def open_combined_sstice_dataset(self):
    """
    Return xarray dataset of combined SST and sea ice dataset
    """
    if self.do_timers: timer_start = perf_counter()
    if self.sstice_combined_file is None :
        # if sstice_combined_file is not defined then we need to 
        # combine the individual SST and sea ice data files
        if self.sst_file is None and self.ice_file is None: 
            raise ValueError('sst_file and ice_file must be set if sstice_combined_file is not set!')
        # Combine individual SST and sea ice data files
        ds_sst = xr.open_dataset(self.sst_file).load()
        ds_ice = xr.open_dataset(self.ice_file).load()
        ds_out = xr.merge([ds_sst,ds_ice])
        ds_sst.close(); ds_ice.close()
    else:
        # Open combined sst/ice data file
        ds_out = xr.open_dataset(self.sstice_combined_file)

    if self.do_timers: self.print_timer(timer_start)
    return ds_out
# ------------------------------------------------------------------------------
def sstice_create_dst_grid_file(self,output_grid_spacing=1,force_overwrite=False,
                                verbose=None):
    """
    Create a target grid file to use for remapping the SST and sea ice data.
    The SST and ice data are assumed to exist on the same grid.
    """
    if self.do_timers: timer_start = perf_counter()
    if verbose is None : verbose = self.verbose

    # Define output grid dimensions
    self.sstice_nlat_dst = int( 180/output_grid_spacing )
    self.sstice_nlon_dst = int( 360/output_grid_spacing )

    dst_grid = f'{self.sstice_nlat_dst}x{self.sstice_nlon_dst}'

    # Define the destination grid file to be created
    # (add 's2n' in case input data is same grid with opposite orientation)
    self.sstice_dst_grid_file = f'{self.grid_dir}/scrip_{dst_grid}_s2n.nc'

    # Create the destination grid file
    if force_overwrite or not os.path.isfile(self.sstice_dst_grid_file) :
        if verbose : print(self.verbose_indent+f'\nCreating target grid file for SST and sea ice data...')
        cmd  = f'ncremap'
        cmd += f' --tmp_dir={self.tmp_dir}'
        cmd += f' -G ttl=\'Equi-Angular grid {dst_grid}\'' 
        cmd += f'#latlon={self.sstice_nlat_dst},{self.sstice_nlon_dst}'
        cmd +=  '#lat_typ=uni'
        cmd +=  '#lon_typ=grn_ctr'
        cmd +=  '#lat_drc=s2n'
        cmd += f' -g {self.sstice_dst_grid_file} '
        run_cmd(cmd,verbose,shell=True)

    if self.do_timers: self.print_timer(timer_start)
    return
# ------------------------------------------------------------------------------
def sstice_create_map_file(self,force_overwrite=False,verbose=None):
    """
    Create a mapping file to be used for SST and sea ice data
    """
    if self.do_timers: timer_start = perf_counter()
    if verbose is None : verbose = self.verbose

    src_grid = f'{self.sstice_nlat_src}x{self.sstice_nlon_src}'
    dst_grid = f'{self.sstice_nlat_dst}x{self.sstice_nlon_dst}'

    self.sstice_map_file = f'{self.map_dir}/map_{src_grid}_to_{dst_grid}_s2n.nc'

    # Generate mapping file
    if force_overwrite or not os.path.isfile(self.sstice_map_file) :
        if verbose : print(self.verbose_indent+f'\nCreating mapping file for SST and sea ice data...')
        cmd  = f'ncremap -a fv2fv'
        cmd += f' --src_grd={self.sstice_src_grid_file}'
        cmd += f' --dst_grd={self.sstice_dst_grid_file}'
        cmd += f' --map_file={self.sstice_map_file}'
        run_cmd(cmd,verbose,shell=True,prepend_line=False)

    if self.do_timers: self.print_timer(timer_start)
    return
# ------------------------------------------------------------------------------
def sstice_slice_and_remap(self,output_file_name,
                           time_slice_method='match_atmos',
                           atm_file=None,verbose=None):
    """
    Horizontally remap the SST and sea ice data after time slicing and 
    combining SST and sea ice files (if they are separate)
    Extract a temporal subset of the SST and sea ice data
    and combine into new temporary file for later regridding
    supported time_slice_method options:
    - initial           use initial time
    - match_atmos       use corresponding time in atmos initial condition
    - use_all           use all times available - simplest way to have transient SST
    """
    if self.do_timers: timer_start = perf_counter()
    if verbose is None : verbose = self.verbose
    if verbose : print(self.verbose_indent+f'\nTime slicing {self.sstice_name} SST and sea ice data...')

    check_dependency('ncatted')
    check_dependency('ncremap')

    # Define temporary file to hold the time sliced data for regridding
    sstice_tmp_file_name = f'{self.tmp_dir}/tmp_sstice_timeslice_data.nc'

    # Check that the time_slice_method is supported
    if time_slice_method not in ['initial','match_atmos','use_all'] :
        raise ValueError(f'time_slice_method: {time_slice_method} is not a supported method')

    if time_slice_method == 'initial' :
        time_slice = 0
    
    if time_slice_method == 'match_atmos' :
        if atm_file is None : 
            raise ValueError('atm_file can not be None for time_slice_method= \'match_atmos\'')
        sst_file = self.get_sst_file()
        ds_sst = xr.open_dataset(sst_file)
        ds_atm = xr.open_dataset(atm_file)
        time_chk = ds_atm['time'].values[0] == ds_sst['time'].values
        ds_sst.close(); ds_atm.close()
        time_slice = np.flatnonzero(time_chk)
        if time_slice.size == 0 :
            raise ValueError('No matching time slice found between SST and atmos data!')
        # If there are multiple matching times, just use the first instance
        if time_slice.size > 1 : time_slice = time_slice[0]
    
    if time_slice_method == 'use_all' :
        sst_file = self.get_sst_file()
        ds_sst = xr.open_dataset(sst_file)
        time_slice = slice(0,len(ds_sst.time)+1)
        ds_sst.close()

    # Load time sliced data
    ds_out = self.open_combined_sstice_dataset()
    if 'time_slice' in locals(): ds_out = ds_out.isel(time=time_slice).load()
    # Drop any extra variables
    for var in ds_out.variables :
        if var not in [self.sst_name,self.ice_name] and var not in ds_out.coords : 
            ds_out = ds_out.drop_vars(var)
    # write out to temporary file
    ds_out.to_netcdf(sstice_tmp_file_name,format=xarray_sst_nc_format)
    ds_out.close()

    # Replace nan values with missing_value to avoid remapping issues
    if self.sstice_name=='NOAA':
        missing_value = 1e36
        # ds = xr.open_dataset(sstice_tmp_file_name).load()
        # ds[self.sst_name] = ds[self.sst_name].where( ds[self.sst_name].values != np.nan, missing_value )
        # ds.to_netcdf(sstice_tmp_file_name,format=xarray_sst_nc_format)
        # ds.close()
        cmd = f'ncatted -h -O -a   xxxx,o,f,{missing_value} {sstice_tmp_file_name}'
        run_cmd(cmd.replace('xxxx',   f'_FillValue,{self.sst_name}'),verbose,shell=True,prepend_line=False)
        run_cmd(cmd.replace('xxxx',   f'_FillValue,{self.ice_name}'),verbose,shell=True,prepend_line=False)
        run_cmd(cmd.replace('xxxx',f'missing_value,{self.sst_name}'),verbose,shell=True,prepend_line=False)
        run_cmd(cmd.replace('xxxx',f'missing_value,{self.ice_name}'),verbose,shell=True,prepend_line=False)
    # if self.sstice_name=='ERA5':
    #     missing_value = 1e36
    #     # ds = xr.open_dataset(sstice_tmp_file_name).load()
    #     # ds.fillna(value=missing_value)
    #     # ds.to_netcdf(sstice_tmp_file_name,format=xarray_sst_nc_format)
    #     # ds.close()
    #     cmd = f'ncatted -h -O -a   xxxx,o,f,{missing_value} {sstice_tmp_file_name}'
    #     run_cmd(cmd.replace('xxxx',   f'_FillValue,{self.sst_name}'),verbose,shell=True,prepend_line=False)
    #     run_cmd(cmd.replace('xxxx',   f'_FillValue,{self.ice_name}'),verbose,shell=True,prepend_line=False)
    #     run_cmd(cmd.replace('xxxx',f'missing_value,{self.sst_name}'),verbose,shell=True,prepend_line=False)
    #     run_cmd(cmd.replace('xxxx',f'missing_value,{self.ice_name}'),verbose,shell=True,prepend_line=False)

    if verbose : print(self.verbose_indent+f'\nRemapping {self.sstice_name} SST and sea ice data...')

    # make sure output file is deleted (overwrite flag not working?)
    if os.path.isfile(output_file_name): 
        run_cmd(f'rm {output_file_name}',verbose)

    # remap the SST data onto the target grid for the model
    cmd =  f'ncremap'
    cmd += f" --nco_opt='-O --no_tmp_fl --hdr_pad={hdr_pad}' "
    cmd += f' --vars={self.sst_name},{self.ice_name} '
    cmd += f' --map_file={self.sstice_map_file} '
    cmd += f' --in_file={sstice_tmp_file_name} '
    cmd += f' --out_file={output_file_name} '
    cmd += f' --fl_fmt={ncremap_file_fmt} '
    run_cmd(cmd,verbose,shell=True)

    # delete the temporary file
    run_cmd(f'rm {sstice_tmp_file_name}',verbose)

    if self.do_timers: self.print_timer(timer_start)
    return
# ------------------------------------------------------------------------------
def sstice_rename_vars(self, output_file_name, new_sst_name='SST_cpl',
                       new_ice_name='ice_cov', verbose=None):
    """
    Rename sst and sea icea variables and remove unnecessary variables
    """
    if self.do_timers: timer_start = perf_counter()
    if verbose is None : verbose = self.verbose
    if verbose : print(self.verbose_indent+'\nRenaming SST and sea ice variables...')

    check_dependency('ncrename')
    check_dependency('ncks')
    check_dependency('ncatted')

    # rename variables
    cmd = f'ncrename --hst'
    cmd+= f' --variable {self.sst_name},{new_sst_name}'
    cmd+= f' --variable {self.ice_name},{new_ice_name}'
    cmd+= f' {output_file_name}'
    run_cmd(cmd,verbose,shell=True)

    # Make sure dimensions names are correct
    cmd = f'ncrename --hst'
    cmd+= f' --dimension .longitude,lon'
    cmd+= f' --dimension .latitude,lat'
    cmd+= f' {output_file_name}'
    run_cmd(cmd,verbose,shell=True)

    # Drop unnecessary variables
    run_cmd(f'ncks -C -O -x -v area,gw,lat_bnds,lon_bnds {output_file_name} {output_file_name}',
            verbose,prepend_line=False,shell=True)

    # Drop unnecessary attributes
    global_att_list = ['history_of_appended_files']
    for att in global_att_list:
        run_cmd(f'ncatted -O -a {att},global,d,, {output_file_name} {output_file_name}',
                verbose,prepend_line=False)

    # Reset the history attribute
    run_cmd(f'ncatted -h -O -a history,global,o,c, {output_file_name} {output_file_name}',
            verbose,prepend_line=False)

    if self.do_timers: self.print_timer(timer_start)
    return
# ------------------------------------------------------------------------------
def sstice_adjustments(self, output_file_name, verbose=None):
    """
    Perform miscellaneous adjustments to the final SST/ice file
    - make sure SST units are Celsius
    - limit sea ice fraction  
    - interpolate to fill in missing SST data, extrapolate where necessary (i.e. Antarctica)
    - make sure time is a coordinate
    - add date and datesec variables
    """
    if self.do_timers: timer_start = perf_counter()
    if verbose is None : verbose = self.verbose
    if verbose : print(self.verbose_indent+'\nAdjusting SST and sea ice data values...')

    # Open remapped and combined data file
    ds = xr.open_dataset(output_file_name).load()

    # Convert units to Celsius
    if 'units' not in ds['SST_cpl'].attrs:
        ds['SST_cpl'].attrs['units'] = 'deg_C'
    if ds['SST_cpl'].attrs['units'] in ['K','degrees_K','deg_K','Kelvin'] :
        ds['SST_cpl'] = ds['SST_cpl']-hc.tk_zero
        ds['SST_cpl'].attrs['units'] = 'deg_C'

    # Set invalid values to np.nan before using interpolate_na()
    ds['SST_cpl'] = xr.where( np.fabs(ds['SST_cpl'].values) < 999, ds['SST_cpl'], np.nan)

    # fill in missing SST values over continents by linearly interpolating
    ds['SST_cpl'] = ds['SST_cpl'].interpolate_na(dim='lon',period=360)

    # extrapolate in latitude direction to deal with poles (requires scipy)
    ds['SST_cpl'] = ds['SST_cpl'].interpolate_na(dim='lat',method='nearest'
                                                ,fill_value='extrapolate')

    # Fill missing ice values with zero
    ds['ice_cov'] = ds['ice_cov'].where( ds['ice_cov']>=0, 0)
    ds['ice_cov'] = ds['ice_cov'].where( ds['ice_cov']<=2, 0)

    # Also limit ice values slightly above 1 that might result from remap
    ds['ice_cov'] = ds['ice_cov'].where( ds['ice_cov']<=1, 1)

    # Make sure time is a coordinate and dimension of the dataset
    if 'time' not in ds.dims : ds = ds.expand_dims('time',axis=0)
    if 'time' not in ds.coords : ds = ds.assign_coords(coords={'time':ds['time']})

    # # For single time value, append extra dummy time for temporal interpolation
    # # STILL NOT SURE IF WE NEED THIS OR NOT
    # if ds['time'].size == 1:
    #     dt = np.timedelta64(10,'D')
    #     ds_dummy = ds.copy(deep=True)
    #     ds_dummy['time'] = ds_dummy['time']-dt
    #     ds = xr.concat([ds_dummy,ds],dim='time')
    #     ds_dummy['time'] = ds_dummy['time']+dt*2
    #     ds = xr.concat([ds,ds_dummy],dim='time')
    #     ds_dummy.close()

    # Create date and datesec variables
    time_index = pd.DatetimeIndex( ds['time'].values )
    date = np.array( time_index.year*1e4+time_index.month*1e2+time_index.day, dtype=int )
    datesec = np.array( time_index.second, dtype=int )

    # # Create alternate time index - doesn't work
    # cftime_index = pd.to_datetime(time_index).astype('cftime.DatetimeNoLeap')
    # cftime_index = pd.to_datetime(time_index)

    # # Create cftime index - doesn't work
    # cftime_index = cftime.datetime(year=time_index.year[0]
    #                               ,month=time_index.month[0]
    #                               ,day=time_index.day[0]
    #                               ,hour=time_index.hour[0]
    #                               ,minute=time_index.minute[0]
    #                               ,second=time_index.second [0] )


    # Add date and datesec variables to dataset
    ds['date'] = xr.DataArray( date, coords={'time':ds['time'].values}, dims=['time'] )
    ds['date'].attrs['long_name'] = 'current date (YYYYMMDD)'
    ds['datesec'] = xr.DataArray( datesec, coords={'time':ds['time'].values}, dims=['time'] )
    ds['datesec'].attrs['long_name'] = 'current seconds of current date'

    # --------------------------------------------------------------------------
    # Experimental - Open climatological data and use it to replace time coordinate
    # --------------------------------------------------------------------------
    # ds2 = xr.open_dataset('/global/cfs/cdirs/acme/inputdata/atm/cam/sst/sst_HadOIBl_bc_1x1_clim_c101029.nc')
    # # ds.drop('time')
    # ds = ds.assign_coords(coords={'time':ds2['time'][0:2]})
    # # ds['time'] = ds2['time'][0:1]
    # if 'time' not in ds.dims : ds = ds.expand_dims('time',axis=0)
    # if 'time' not in ds.coords : ds = ds.assign_coords(coords={'time':ds['time']})
    # ds['date'] = ds2['date'][0:2]
    # ds['datesec'] = ds2['datesec'][0:2]
    # ds2.close()
    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------

    # reset attributes to avoid error reading file during model run
    ds['SST_cpl'].attrs = {'long_name':'Sea Surface Temperature','units':'deg_C'}
    ds['ice_cov'].attrs = {'long_name':'Sea Ice Fraction','units':'fraction'}
    ds['lat'].attrs = {'long_name':'latitude','units':'degrees_north'}
    ds['lon'].attrs = {'long_name':'longitude','units':'degrees_east'}

    # Write back to final file
    ds.to_netcdf(output_file_name
                ,unlimited_dims=['time'] 
                ,encoding={'time':{'dtype':'float64'}}
                ,format=xarray_sst_nc_format
                )
    ds.close()

    run_cmd(f'ncatted --hst -a calendar,time,m,c,\'365_day\' {output_file_name}',
            verbose,prepend_line=False,shell=True)

    run_cmd(f'ncatted -O -a _FillValue,time,d,, {output_file_name}',
            verbose,prepend_line=False)

    # A confusing issue occurs with the dimension ordering 
    # the commands below provide a clunky workaround
    tmp_file = output_file_name.replace('.nc','.ncpdq_tmp.nc')
    cmd1 = f'ncpdq --permute=lon,lat,time --overwrite {output_file_name} {tmp_file}'
    cmd2 = f'ncpdq --permute=time,lat,lon --overwrite {tmp_file} {output_file_name}'
    run_cmd(cmd1,verbose)
    run_cmd(cmd2,verbose,prepend_line=False)
    run_cmd(f'rm {tmp_file} ',verbose,prepend_line=False)

    # Reset the history attribute
    run_cmd(f'ncatted -h -O -a history,global,o,c, {output_file_name} {output_file_name}',
            verbose,prepend_line=False)

    if self.do_timers: self.print_timer(timer_start)
    return
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
