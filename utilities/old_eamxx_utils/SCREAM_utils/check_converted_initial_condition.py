import numpy as np, xarray as xr

# var_name,file_path = 'o3','/p/lustre1/hannah6/hiccup_scratch/gas_constituents_ne30np4L128_20220711.permuted.nc'
# var_name,file_path = 'o3_volume_mix_ratio','/p/lustre1/hannah6/hiccup_scratch/HICCUP.atm_era5.2006-08-09.ne30np4.L128.converted.nc'
# var_name,file_path = 'O3','/p/lustre1/hannah6/hiccup_scratch/HICCUP.atm_era5.2006-08-09.ne30np4.L128.nc'

# var_name,file_path = 'o3','/p/lustre1/hannah6/2024-saomai-data/gas_constituents_ne30np4L128_20220711.permuted.nc'
# var_name,file_path = 'o3','/p/lustre1/hannah6/2024-saomai-data/gas_constituents_ne0np4-saomai-128x8-pg2-L128_20240110.nc'
# var_name,file_path = 'o3','/p/lustre1/hannah6/2024-saomai-data/gas_constituents_ne0np4-saomai-128x8-pg2-L128_20240110.permuted.nc'
var_name,file_path = 'o3_volume_mix_ratio','/p/lustre1/hannah6/2024-saomai-data/HICCUP.atm_era5.2006-08-09.Saomai_2006_ne128x8_lon130E_lat25N.L128.converted.nc'

print()
print(f'file: {file_path}')
print(f'var:  {var_name}')
print()

ds = xr.open_dataset(file_path)

data = ds[var_name]

any_isnan = np.any(np.isnan(data.values))
any_isinf = np.any(np.isinf(data.values))
any_isneg = np.any( data.values<0 )

print()
print(f'  any isnan: {any_isnan}')
print(f'  any isinf: {any_isnan}')
print(f'  any isneg: {any_isneg}')
print()

# if True:
if any_isneg:
  for e in reversed(range(10)):
    p = -1*e
    chk_val = -1* np.power(10.,p)
    chk = np.any( data.values < chk_val )
    print(f'  any less than -10^{p:2d} ( {chk_val:+1.10f} ) : {chk}')

print()