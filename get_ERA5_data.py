#!/usr/bin/env python
import os
import cdsapi
server = cdsapi.Client()

get_atm = True
get_sfc = True
get_lnd = True

yr_list = ['2018']
mn_list = ['01']
dy_list = ['01']

# time_list = ['00:00','03:00','06:00','09:00','12:00','15:00','18:00','21:00']
time_list = ['00:00']

lev = [ '50', '70','100','125','150','175','200','225','250','300','350'
      ,'400','450','500','550','600','650','700','750','775','800','825'
      ,'850','875','900','925','950','975','1000']

output_path = os.getenv('PWD')

output_file_plv = output_path+'ERA5.HICCUP_TEST.atm.nc'
output_file_sfc = output_file_plv.replace('atm.nc','sfc.nc')
output_file_lnd = output_file_plv.replace('atm.nc','lnd.nc')

#-------------------------------------------------------------------------------
# atmossphere pressure level data
if get_atm:
    server.retrieve('reanalysis-era5-pressure-levels',{
        'product_type'  : 'reanalysis',
        'pressure_level': lev,
        'time'          : time_list,
        'day'           : dy_list,
        'month'         : mn_list,
        'year'          : yr_list,
        'format'        : 'netcdf',
        'variable'      : ['temperature'
                          ,'specific_humidity'
                          ,'geopotential'
                          ,'u_component_of_wind'
                          ,'v_component_of_wind'
                          ],
    }, output_file_plv)
#-------------------------------------------------------------------------------
# atmosphere surface data
if get_sfc:
    server.retrieve('reanalysis-era5-single-levels',{
        'product_type'  : 'reanalysis',
        'time'          : time_list,
        'day'           : dy_list,
        'month'         : mn_list,
        'year'          : yr_list,
        'format'        : 'netcdf',
        'variable'      : ['surface_pressure'
                          ,'skin_temperature'
                          ,'sea_surface_temperature'
                          ],
    }, output_file_sfc)
#-------------------------------------------------------------------------------
# land model data
if get_lnd:
    server.retrieve('reanalysis-era5-land',{
            'format': 'netcdf',
            'variable': [
                'leaf_area_index_high_vegetation', 'leaf_area_index_low_vegetation', 'skin_reservoir_content',
                'skin_temperature', 'snow_albedo', 'snow_cover',
                'snow_density', 'snow_depth', 'snow_depth_water_equivalent',
                'soil_temperature_level_1', 'soil_temperature_level_2', 'soil_temperature_level_3',
                'soil_temperature_level_4', 'temperature_of_snow_layer', 'volumetric_soil_water_layer_1',
                'volumetric_soil_water_layer_2', 'volumetric_soil_water_layer_3', 'volumetric_soil_water_layer_4',
            ],
            'time'          : time_list,
            'day'           : dy_list,
            'month'         : mn_list,
            'year'          : yr_list,
    },output_file_lnd)
#-------------------------------------------------------------------------------
