
import os, datetime, xarray as xr, pandas as pd

bdate_str,edate_str = '2020-01-21 00','2020-01-26 21' # 2024 SCREAM - autocalibration

# Create list of date and time 
beg_date = datetime.datetime.strptime(bdate_str, '%Y-%m-%d %H')
end_date = datetime.datetime.strptime(edate_str, '%Y-%m-%d %H')
datetime_list = pd.date_range(beg_date, end_date, freq='1D')

data_root = '/lustre/orion/cli115/proj-shared/hannah6/scream_scratch/nudge_data'
dst_horz_grid = 'ne128pg2'
dst_vert_grid = 'L128'

ref_date = '2000-01-01'

# ------------------------------------------------------------------------------
for t in datetime_list:

    nudge_datetime = t.strftime("%Y-%m-%d")

    data_file_name = f'{data_root}/HICCUP.nudging_uv_era5.{nudge_datetime}.{dst_horz_grid}.{dst_vert_grid}.nc'

    print();print(' '*4+f'{data_file_name}')

    print();print(' '*4+'Transposing data dimensions...')

    cmd = f'ncpdq -5 -a time,ncol,lev {data_file_name}'
    hdc.run_cmd(cmd,verbose=True,shell=True,prefix=hdc.verbose_indent)

    # ds = xr.open_dataset(data_file_name)
    # transpose_flag = False
    # for v in ds.variables:
    #   if ds[v].dims==('time','lev','ncol'):
    #     print(' '*6+f'var: {v}')
    #     transpose_flag = True
    #     ds[v] = ds[v].transpose('time','ncol','lev')
    # if transpose_flag:
    #   time_encoding_dict = {'time':{'units': f'hours since {ref_date} 00:00:00'}}
    #   ds.to_netcdf(data_file_name,format='NETCDF3_64BIT_DATA',mode='w',encoding=time_encoding_dict)
    #   # ds.to_netcdf(data_file_name,format='NETCDF4',mode='w',encoding=time_encoding_dict)
    # ds.close()
# ------------------------------------------------------------------------------