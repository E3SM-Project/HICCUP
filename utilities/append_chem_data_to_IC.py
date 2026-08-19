#!/usr/bin/env python
import os, datetime, subprocess as sp
import xarray as xr, numpy as np
#---------------------------------------------------------------------------------------------------
class tclr: END,RED,GREEN,MAGENTA,CYAN = '\033[0m','\033[31m','\033[32m','\033[35m','\033[36m'
#---------------------------------------------------------------------------------------------------
def run_cmd(cmd,suppress_output=False):
  if suppress_output : cmd = cmd + ' > /dev/null'
  msg = tclr.GREEN + cmd + tclr.END ; print(f'\n{msg}')
  os.system(cmd); return
#---------------------------------------------------------------------------------------------------
''' commands to remap source file for chemsitry

GRID_ROOT=/global/cfs/projectdirs/m3312/whannah/HICCUP/files_grid
MAPS_ROOT=/global/cfs/projectdirs/m3312/whannah/HICCUP/files_map
SRC_GRID=${GRID_ROOT}/exodus_ne30.g
DST_GRID=${GRID_ROOT}/exodus_ne256.g
MAP_FILE=${MAPS_ROOT}/map_ne30np4_to_ne256np4_se2se_20250702.nc

mv /global/cfs/projectdirs/m3312/whannah/HICCUP/files_map/map_ne30np4_to_ne256_traave_20241001.nc /global/cfs/projectdirs/m3312/whannah/HICCUP/files_map/map_ne30np4_to_ne256_se2se_20250702.nc
mv /global/cfs/projectdirs/m3312/whannah/HICCUP/files_map/map_ne30np4_to_ne256_se2se_20250702.nc /global/cfs/projectdirs/m3312/whannah/HICCUP/files_map/map_ne30np4_to_ne256np4_se2se_20250702.nc

ncremap -a se2se --src_grd=${SRC_GRID} --dst_grd=${DST_GRID} --map_file=${MAP_FILE}

chm_root=/global/cfs/cdirs/m3312/whannah/e3smv3_amip/v3.LR.amip_0201/archive/rest/2023-01-01-00000
chm_src_file=$chm_root/v3.LR.amip_0201.eam.i.2023-01-01-00000.nc
chm_dst_file=$chm_root/v3.LR.amip_0201.eam.i.2023-01-01-00000.remap_ne256.nc

ncremap -m $MAP_FILE -i $chm_src_file -o $chm_dst_file
'''
#---------------------------------------------------------------------------------------------------

src_root = '/global/cfs/projectdirs/e3sm/whannah/HICCUP'
chm_root = '/global/cfs/cdirs/m3312/whannah/e3smv3_amip/v3.LR.amip_0201/archive/rest/2023-01-01-00000/'

# src_file = f'{src_root}/HICCUP.atm_era5.2012-10-01.ne256np4.L80.nc'
# dst_file = f'{src_root}/HICCUP.atm_era5.2012-10-01.ne256np4.L80.w-chem.nc'
# chm_file = f'{chm_root}/v3.LR.amip_0201.eam.i.2023-01-01-00000.remap_ne256.nc'

src_file = f'{src_root}/HICCUP.atm_era5.2012-10-01.ne30np4.L80.nc'
dst_file = f'{src_root}/HICCUP.atm_era5.2012-10-01.ne30np4.L80.w-chem.nc'
chm_file = f'{chm_root}/v3.LR.amip_0201.eam.i.2023-01-01-00000.nc'

print()
print(f'  src_file: {tclr.GREEN}{src_file}{tclr.END}')
print(f'  dst_file: {tclr.GREEN}{dst_file}{tclr.END}')
print(f'  chm_file: {tclr.GREEN}{chm_file}{tclr.END}')
print()

# exit()

# ------------------------------------------------------------------------------
# # use this to check basic things about the datasets
# print(); print( xr.open_dataset(src_file) )
# print(); print( xr.open_dataset(dst_file) )
# print(); print( xr.open_dataset(chm_file) )
# print()
# exit()
# ------------------------------------------------------------------------------
# # use this to compare the new file after it has been created
# var = 'O3' # PS / O3
# da_src = xr.open_dataset(src_file)[var].values
# da_dst = xr.open_dataset(dst_file)[var].values
# da_chm = xr.open_dataset(chm_file)[var].values
# print()
# print(f'src {var} min / max : {np.nanmin(da_src):10.2g}  /  {np.nanmax(da_src):10.2g}')
# print(f'dst {var} min / max : {np.nanmin(da_dst):10.2g}  /  {np.nanmax(da_dst):10.2g}')
# print(f'chm {var} min / max : {np.nanmin(da_chm):10.2g}  /  {np.nanmax(da_chm):10.2g}')
# print()
# exit()
# ------------------------------------------------------------------------------
# # use this to check for what variables might be missing from the new file (dst_file)
# ds_dst = xr.open_dataset(dst_file)
# ds_chm = xr.open_dataset(chm_file)
# var_list = []
# for v in ds_chm.data_vars:
#   if v not in ds_dst.data_vars:
#     if 'ncol_d' in ds_chm[v].dims:
#       var_list.append(v)
# if var_list!=[]:
#   print('variables in source chem file not found in destination file: ')
#   for v in var_list:
#     print(v)
#     # print(f'{v:12}  {ds_chm[v].dims}')
# exit()
# ------------------------------------------------------------------------------
def check_and_modify_fill_value(da):
  fill_value = 1e30 # this will replace _FillValue=NaN
  # check attributes
  if hasattr(da, '_FillValue'):
    if np.isnan(da._FillValue):
      da._FillValue = fill_value
  # check encoding
  if '_FillValue' in da.encoding:
    da.encoding['_FillValue'] = fill_value
    # if np.isnan(da.encoding['_FillValue']): da.encoding['_FillValue'] = fill_value
  return da
# ------------------------------------------------------------------------------

chm_var_list = []
chm_var_list.append('CO')
chm_var_list.append('C2H6')
chm_var_list.append('C3H8')
chm_var_list.append('CH3COCH3')
chm_var_list.append('E90')
chm_var_list.extend(['N2OLNZ','NOYLNZ','CH4LNZ','H2OLNZ'])
chm_var_list.append('DMS')
chm_var_list.append('SO2')
chm_var_list.append('H2SO4')
chm_var_list.extend(['so4_a1','so4_a2','so4_a3','so4_a5'])
chm_var_list.extend(['pom_a1','pom_a3','pom_a4'])
chm_var_list.extend(['bc_a1','bc_a3','bc_a4'])
chm_var_list.extend(['dst_a1','dst_a3'])
chm_var_list.extend(['ncl_a1','ncl_a2','ncl_a3'])
chm_var_list.extend(['mom_a1','mom_a2','mom_a3'])
chm_var_list.extend(['mom_a4','num_a1','num_a2','num_a3','num_a4','num_a5'])
# chm_var_list.append('O3') # we can keep O3 from ERA5
chm_var_list.append('OH')
chm_var_list.append('HO2')
chm_var_list.append('H2O2')
chm_var_list.append('CH2O')
chm_var_list.append('CH3O2')
chm_var_list.append('CH3OOH')
chm_var_list.append('NO')
chm_var_list.append('NO2')
chm_var_list.append('NO3')
chm_var_list.append('N2O5')
chm_var_list.append('HNO3')
chm_var_list.append('HO2NO2')
chm_var_list.append('PAN')
chm_var_list.append('C2H5O2')
chm_var_list.append('C2H5OOH')
chm_var_list.append('CH3CHO')
chm_var_list.append('CH3CO3')
chm_var_list.append('C2H4')
chm_var_list.append('ROHO2')
chm_var_list.append('ISOP')
chm_var_list.append('ISOP_VBS')
chm_var_list.append('ISOPO2')
chm_var_list.append('C10H16')
chm_var_list.append('MVKMACR')
chm_var_list.append('MVKO2')
chm_var_list.extend(['soa_a1','soa_a2','soa_a3'])
chm_var_list.extend(['SOAG0','SOAG15','SOAG24','SOAG31','SOAG32','SOAG33','SOAG34','SOAG35'])

# these were identified later as possibly needed
chm_var_list.extend(['BVRIM','CLDRIM','NUMICE','NUMLIQ','NUMRAI','RAINQM'])
# chm_var_list.append('SICTHK')
# chm_var_list.append('TSICE')
# chm_var_list.append('cosp_ht_bnds')
# chm_var_list.append('cosp_htmisr_bnds')
# chm_var_list.append('cosp_prs_bnds')
# chm_var_list.append('cosp_reffice_bnds')
# chm_var_list.append('cosp_reffliq_bnds')
# chm_var_list.append('cosp_sr_bnds')
# chm_var_list.append('cosp_tau_bnds')
# chm_var_list.append('cosp_tau_modis_bnds')
# chm_var_list.append('cosp_temp_bnds')


# ------------------------------------------------------------------------------
# attempt #1 - the PS variable is erroneously overwritten

# chm_var_list_str = ','.join(chm_var_list)

# # make a copy of the main eam.i file
# run_cmd(f'cp {src_file} {dst_file}')

# # append the chem vars
# run_cmd(f'ncks -A -v {chm_var_list_str} --exclude PS {chm_file} {dst_file}')

# ------------------------------------------------------------------------------
# attempt #2 - errors due to type mismatch and ncol_d/ncol issues

# alt_var_list = []
# alt_var_list.append('CLDICE')
# alt_var_list.append('CLDLIQ')
# # alt_var_list.append('ICEFRAC')
# alt_var_list.append('O3')
# alt_var_list.append('PHIS')
# alt_var_list.append('PS')
# alt_var_list.append('Q')
# alt_var_list.append('SNOWHICE')
# alt_var_list.append('T')
# alt_var_list.append('TS')
# alt_var_list.append('TS1')
# alt_var_list.append('TS2')
# alt_var_list.append('TS3')
# alt_var_list.append('TS4')
# alt_var_list.append('U')
# alt_var_list.append('V')
# alt_var_list.append('date')
# alt_var_list.append('datesec')
# alt_var_list.append('nbdate')
# alt_var_list.append('nbsec')
# alt_var_list.append('ndbase')
# alt_var_list.append('ndcur')
# alt_var_list.append('nsbase')
# alt_var_list.append('nscur')
# alt_var_list.append('nsteph')
# alt_var_list.append('time')
# alt_var_list_str = ','.join(alt_var_list)

# # make a copy of the "chem" eam.i file
# run_cmd(f'cp {chm_file} {dst_file}')

# # change format of variables in destination file
# # ???

# # append the non-chem vars
# run_cmd(f'ncks -A -v {alt_var_list_str} {src_file} {dst_file}')

# ------------------------------------------------------------------------------

# make a copy of the main eam.i file
run_cmd(f'cp {src_file} {dst_file}')


# ds_src = xr.open_dataset(src_file)
ds_dst = xr.open_dataset(dst_file)
ds_chm = xr.open_dataset(chm_file)

print()
print('copying data from chem file to new eam.i file...')

# for var in chm_var_list:
for var in ['soa_a2']:
  if var in ds_chm:
    print(f'  copying variable: {tclr.CYAN}{var}{tclr.END}')
    da_chm = ds_chm[var]
    # rename ncol coordinate
    da_chm = da_chm.rename({'ncol_d':'ncol'})
    # replace time cooridate to avoid overwriting destination time coord
    da_chm['time'] = ds_dst['time']
    ds_dst[var] = da_chm
  else:
    print(f'  skipping variable: {tclr.RED}{var}{tclr.END} (not found in source data file)')


# update _FillValue attribute/encoding
# note - this doesn't seem to matter until we add the chemistry tracers (not sure why)
for var in ds_dst.variables:
  ds_dst[var] = check_and_modify_fill_value(ds_dst[var])

ds_dst.load()
# ds_dst.to_netcdf(dst_file)
ds_dst.to_netcdf(dst_file,encoding={'_FillValue':1e30})
ds_dst.close()

print('done.')
print()

run_cmd(f'ncks -O  -5 {dst_file} {dst_file}')

# ------------------------------------------------------------------------------
# check results
var = 'PS'
da_src = xr.open_dataset(src_file,decode_times=False)[var].values
da_dst = xr.open_dataset(dst_file,decode_times=False)[var].values
da_chm = xr.open_dataset(chm_file,decode_times=False)[var].values
print()
print(f'src {var} min / max : {np.nanmin(da_src):10.2f}  /  {np.nanmax(da_src):10.2f}')
print(f'dst {var} min / max : {np.nanmin(da_dst):10.2f}  /  {np.nanmax(da_dst):10.2f}')
print(f'chm {var} min / max : {np.nanmin(da_chm):10.2f}  /  {np.nanmax(da_chm):10.2f}')
print()
exit()

# ------------------------------------------------------------------------------
