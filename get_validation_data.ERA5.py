#!/usr/bin/env python
# Script for downloading ERA5 data for hindcast validation
# links to CDS web interface:
#   https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-pressure-levels
#   https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels
# A list of available variables can also be found in the ERA5 documentation:
#   https://confluence.ecmwf.int/display/CKB/ERA5%3A+data+documentation

# Typically this data and the hindcast data should be remapped to a common grid
# So here are some helpful commands for putting everything on a 1 deg grid:

# Unpack the data: FILE=data_scratch/ERA5_validation.Q.2016-08-01.nc ; ncpdq --ovr --unpack $FILE $FILE

# Generate 1-deg grid file: ncremap -a tempest -G ttl='Equi-Angular grid 180x360'#latlon=180,360#lat_typ=uni#lat_drc=s2n#lon_typ=grn_ctr -g $HOME/HICCUP/files_grid/scrip_180x360_s2n.nc
# Generate 2-deg grid file: ncremap -a tempest -G ttl='Equi-Angular grid 90x180'#latlon=90,180#lat_typ=uni#lat_drc=s2n#lon_typ=grn_ctr -g $HOME/HICCUP/files_grid/scrip_90x180_s2n.nc

# Generate map file for ERA5: ncremap --alg_typ=tempest -a fv2fv --src_grd=./files_grid/scrip_ERA5_721x1440.nc --dst_grd=./files_grid/scrip_90x180_s2n.nc --map_file=./files_mapping/map_721x1440_n2s_to_90x180_s2n.nc

# Generate map file for E3SM ne30np4: ncremap --alg_typ=tempest --src_grd=./files_grid/exodus_ne30.g --dst_grd=./files_grid/scrip_90x180_s2n.nc --map_file=/global/homes/w/whannah/maps/map_ne30np4_to_90x180.nc --wgt_opt='--in_type cgll --in_np 4 --out_type fv --out_np 2 --out_double'

# Generate map file for E3SM ne30pg2: ncremap --alg_typ=tempest --src_grd=./files_grid/exodus_ne30pg2.nc --dst_grd=./files_grid/scrip_90x180_s2n.nc --map_file=/global/homes/w/whannah/maps/map_ne30pg2_to_90x180.nc --wgt_opt='--in_type fv --in_np 2 --out_type fv --out_np 2 --out_double'

# Remap the data: ncremap -m <map file> -i <input file> -o <output file>
# obs example: FILE=data_scratch/ERA5_validation.Z.2016-08-01 ; ncremap -m ./files_mapping/map_721x1440_n2s_to_90x180_s2n.nc -i $FILE.nc -o $FILE.remap_90x180.nc
# hindcast example:
# export MSCRATCH=/global/cscratch1/sd/whannah/e3sm_scratch/cori-knl/
# CASE=E3SM_HINDCAST-TEST_2016-08-01_ne30_FC5AV1C-L_00 ; FILE=$MSCRATCH/$CASE/run/$CASE.cam.h1.2016-08-01-00000 ; ncremap -m $HOME/maps/map_ne30np4_to_90x180.nc -i $FILE.nc -o $FILE.remap_90x180.nc
# CASE=E3SM_HINDCAST-TEST_2016-08-01_ne30pg2_FC5AV1C-L_00 ; FILE=$MSCRATCH/$CASE/run/$CASE.cam.h1.2016-08-01-00000 ; ncremap -m $HOME/maps/map_ne30pg2_to_90x180.nc -i $FILE.nc -o $FILE.remap_90x180.nc


import os
import cdsapi
import datetime
server = cdsapi.Client()

get_atm = False
get_sfc = True


# Build a list of year,month,day values
ndays = 5
sdate = datetime.date(2016, 8, 1)


yr_list,mn_list,dy_list = [],[],[]
# for i in range( (edate-sdate).days + 1):
for i in range( (datetime.timedelta(days=ndays)).days ):
  tdate = sdate + datetime.timedelta(days=i)
  yr_list.append(str(tdate.year).zfill(4))
  mn_list.append(str(tdate.month).zfill(2))
  dy_list.append(str(tdate.day).zfill(2))

time_list = ['00:00','03:00','06:00','09:00','12:00','15:00','18:00','21:00']

# lev = [ '50','100','150','200','300','400','500','600','700','750','800','850','900','950','1000']


atm_var_dict = {}
atm_var_dict.update({'Z':'geopotential'})
atm_var_dict.update({'T':'temperature'})
atm_var_dict.update({'Q':'specific_humidity'})
atm_var_dict.update({'U':'u_component_of_wind'})
atm_var_dict.update({'V':'v_component_of_wind'})

sfc_var_dict = {}
sfc_var_dict.update({'TS':'skin_temperature'})
sfc_var_dict.update({'PS':'surface_pressure'})
# sfc_var_dict.update({'':'sea_surface_temperature'})
# sfc_var_dict.update({'':'sea_ice_cover'})
# sfc_var_dict.update({'':'snow_depth'})

# output_path = os.getenv('PWD')+'/validation_data/'
output_path = '/global/cscratch1/sd/whannah/HICCUP/data/'

#-------------------------------------------------------------------------------
# atmossphere pressure level data
if get_atm:
  for key in atm_var_dict.keys():
    if key=='Z': lev = ['100','500','700']
    if key=='T': lev = ['500','850']
    if key=='U': lev = ['200','850']
    if key=='V': lev = ['200','850']
    if key=='Q': lev = ['850']
    output_file = output_path+f'ERA5_validation.{key}.{yr_list[0]}-{mn_list[0]}-{dy_list[0]}.nc'
    server.retrieve('reanalysis-era5-pressure-levels',{
        'product_type'  : 'reanalysis',
        'pressure_level': lev,
        'time'          : time_list,
        'day'           : dy_list,
        'month'         : mn_list,
        'year'          : yr_list,
        'format'        : 'netcdf',
        'variable'      : [var_dict[key]],
    }, output_file)
#-------------------------------------------------------------------------------
# surface data
if get_sfc:
  for key in sfc_var_dict.keys():
    output_file = output_path+f'ERA5_validation.{key}.{yr_list[0]}-{mn_list[0]}-{dy_list[0]}.nc'
    server.retrieve('reanalysis-era5-single-levels',{
        'product_type'  : 'reanalysis',
        'time'          : time_list,
        'day'           : dy_list,
        'month'         : mn_list,
        'year'          : yr_list,
        'format'        : 'netcdf',
        'variable'      : [sfc_var_dict[key]],
    }, output_file)
#-------------------------------------------------------------------------------
