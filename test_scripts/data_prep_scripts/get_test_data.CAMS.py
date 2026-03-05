#!/usr/bin/env python
# Script for downloading CAMS reanalysis and ERA5 surface data for testing
# See the CAMS data documentation for more information:
#   https://confluence.ecmwf.int/display/CKB/CAMS%3A+Reanalysis+data+documentation
# NOTE: requires cdsapi >= 0.7.0
# ~/.adsapirc should contain:
#   url: https://ads.atmosphere.copernicus.eu/api
#   key: <your-personal-access-token>
#-----------------------------------------------------------------------------
import os, yaml, cdsapi, datetime, pandas as pd
#-----------------------------------------------------------------------------
get_atm = True
get_sfc = True
#-------------------------------------------------------------------------------
# obtain single date for testing
yr,mn,dy,time = '2008','10','01','00:00'

# Small level set for testing
lev = [ '50','100','200','300','400','500','600','700','800','900','1000']

#-----------------------------------------------------------------------------
# Specify output file names
output_path = os.getenv('PWD')+'/test_data'
output_file_atm = output_path+f'/HICCUP_TEST.CAMS.atm.nc'
output_file_sfc = output_path+f'/HICCUP_TEST.CAMS.sfc.nc'
#-------------------------------------------------------------------------------
# atmossphere pressure level data
if get_atm:
    with open(os.getenv('HOME')+'/.adsapirc', 'r') as f: credentials = yaml.safe_load(f)
    ads_server = cdsapi.Client(url=credentials['url'], key=credentials['key'])
    ads_server.retrieve('cams-global-reanalysis-eac4',{
        'type'          : 'an',
        'stream'        : 'oper',
        'data_format'   : 'netcdf',
        'grid'          : [0.75, 0.75],
        'pressure_level': lev,
        'time'          : [time],
        'date': f'{yr}-{mn}-{dy}/{yr}-{mn}-{dy}',
        'variable'      : ['dust_aerosol_0.03-0.55um_mixing_ratio'
                          ,'dust_aerosol_0.55-0.9um_mixing_ratio'
                          ,'dust_aerosol_0.9-20um_mixing_ratio'
                          ],
    }, output_file_atm)
    del ads_server
#-------------------------------------------------------------------------------
# surface data
if get_sfc:
    cds_server = cdsapi.Client()
    cds_server.retrieve('reanalysis-era5-single-levels',{
        'product_type'  : ['reanalysis'],
        'grid'          : [0.75, 0.75],
        'time'          : [time],
        'day'           : [dy],
        'month'         : [mn],
        'year'          : [yr],
        'data_format'   : 'netcdf',
        'variable'      : ['surface_pressure'],
    }, output_file_sfc)
    del cds_server
#-------------------------------------------------------------------------------
