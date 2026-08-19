import xarray as xr, numpy as np, os
from optparse import OptionParser
#-------------------------------------------------------------------------------
help_header = 'usage: ./%prog --file=[file],[file] --var=[var],[var]\n'
help_header += '\nThis script will check the data values files for inf and nan values'
parser = OptionParser(usage=help_header)
parser.add_option('--file',dest='file',default='all',help='comma separated list of files to check')
parser.add_option('--var',dest='var',default='all',help='comma separated list of variables to check')
(opts, args) = parser.parse_args()
#-------------------------------------------------------------------------------
# ifile = '/pscratch/sd/e/ebercosh/SCREAM/IC_20050827/HICCUP.atm_era5.2005-08-27.ne256np4.L128.nc'
# ifile = '/pscratch/sd/w/whannah/HICCUP/test_ebercosh/HICCUP.atm_era5.2005-08-27.ne1024np4.L128.nc'

# var = 'ps'
# var = 'T_mid'

#-------------------------------------------------------------------------------
ds = xr.open_dataset(opts.file).isel(time=0)

data = ds[var]
if 'lev' in data.dims:
  data = data.isel(lev=-1)
max_ind = data.argmax(dim='ncol')
#-------------------------------------------------------------------------------
print()
print(f'  {"ncol":6} : {max_ind.values}')
print(f'  {var   :6} : {data     .isel(ncol=max_ind).values}')
print(f'  {"lat" :6} : {ds["lat"].isel(ncol=max_ind).values}')
print(f'  {"lon" :6} : {ds["lon"].isel(ncol=max_ind).values}')
print()
#-------------------------------------------------------------------------------
