import xarray as xr
import hapy_common as hc

ifile_path = '/global/cfs/cdirs/e3sm/inputdata/atm/cam/topo/USGS-gtopo30_ne30np4pg2_16xdel2.c20200108.nc'

ofile_path = 'USGS-gtopo30_ne30np4pg2_16xdel2.DUMMY.nc'

ds = xr.open_dataset(ifile_path)

print(ds.LANDM_COSLAT)

for key in ds.variables.keys():
  print('  '+key)
  if key in ['PHIS','PHIS_d','LANDFRAC','LANDM_COSLAT','SGH','SGH30']:
    # print()
    # hc.print_stat(ds[key],compact=True)
    ds[key] = ds[key]*0.
    # hc.print_stat(ds[key],compact=True)

print()
print(f'  writing to file: {ofile_path}')
ds.to_netcdf(path=ofile_path,mode='w')
