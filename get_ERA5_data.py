#!/usr/bin/env python
import os
import cdsapi
server = cdsapi.Client()

yr1,mn1,dy1 = 2018,1,1
yr2,mn2,dy2 = 2018,1,2
date_range = f'{yr1}{mn1}{dy1}/to/{yr2}{mn2}{dy2}'

var = ['Ps','Ts','U','T','Q']

lev = ['50','70','100','125','150','175','200'
      ,'225','250','300','350','400','450','500'
      ,'550','600','650','700','750','775','800'
      ,'825','850','875','900','925','950','975','1000'],

output_path = os.getenv('HOME')+'/Data/Obs/ERA5/raw_netcdf/'
output_file = output_path+'ERA5.HICCUP_TEST.nc'

server.retrieve('reanalysis-era5-complete',{
    'product_type'  : 'reanalysis',
    'levtype'       : 'pl', 
    'levelist'      : lev,
    'time'          : ['00:00','03:00','06:00','09:00','12:00','15:00','18:00','21:00'],
    'date'          : date_range,
    'format'        : 'netcdf',
    'variable'      : ['temperature','specific_humidity','geopotential'
                      ,'u_component_of_wind','v_component_of_wind'
                      ,'surface_pressure','skin_temperature'
                      ,'sea_surface_temperature'
                      ],
}, output_file)


# ,'mean_sea_level_pressure'  
