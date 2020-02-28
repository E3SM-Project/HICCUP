#!/usr/bin/env python
# Script for downloading ERA5 pressure level and surface data
# links to CDS web interface:
#   https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-pressure-levels
#   https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels
# A list of available variables can also be found in the ERA5 documentation:
#   https://confluence.ecmwf.int/display/CKB/ERA5%3A+data+documentation

import os
import cdsapi
server = cdsapi.Client()

get_atm = False
get_sfc = True
get_lnd = False

yr_list = ['2018']
mn_list = ['01']
dy_list = ['01']

# time_list = ['00:00','03:00','06:00','09:00','12:00','15:00','18:00','21:00']

# use single time for testing
time_list = ['00:00']

# lev = [ '50', '70','100','125','150','175','200','225','250','300','350'
#       ,'400','450','500','550','600','650','700','750','775','800','825'
#       ,'850','875','900','925','950','975','1000']

# Small level set for testing
lev = [ '50','100','150','200','300','400','500','600'
      ,'700','750','800','850','900','950','1000']

output_path = os.getenv('PWD')+'/data/'

output_file_plv = output_path+'HICCUP_TEST.ERA5.atm.nc'
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
                          ,'ozone_mass_mixing_ratio'
                          ,'specific_cloud_ice_water_content'
                          ,'specific_cloud_liquid_water_content'
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
                          ,'soil_temperature_level_1'
                          ,'soil_temperature_level_2'
                          ,'soil_temperature_level_3'
                          ,'soil_temperature_level_4'
                          # ,'leaf_area_index_high_vegetation'
                          # ,'leaf_area_index_low_vegetation'
                          # ,'skin_reservoir_content'
                          # ,'snow_albedo'
                          # ,'snow_density'
                          ,'snow_depth'
                          ,'temperature_of_snow_layer'
                          # ,'volumetric_soil_water_layer_1'
                          # ,'volumetric_soil_water_layer_2'
                          # ,'volumetric_soil_water_layer_3'
                          # ,'volumetric_soil_water_layer_4'
                          ,'orography'
                          ,'sea_ice_cover'
                          ],
    }, output_file_sfc)
#-------------------------------------------------------------------------------
# land model data
if get_lnd:
    server.retrieve('reanalysis-era5-land',{
            'format': 'netcdf',
            'variable': [
                'leaf_area_index_high_vegetation'
                ,'leaf_area_index_low_vegetation'
                ,'skin_reservoir_content'
                ,'snow_albedo'
                ,'snow_cover'
                ,'snow_density'
                ,'snow_depth'
                ,'snow_depth_water_equivalent'
                # ,'soil_temperature_level_1'
                # ,'soil_temperature_level_2'
                # ,'soil_temperature_level_3'
                # ,'soil_temperature_level_4'
                ,'temperature_of_snow_layer'
                ,'volumetric_soil_water_layer_1'
                ,'volumetric_soil_water_layer_2'
                ,'volumetric_soil_water_layer_3'
                ,'volumetric_soil_water_layer_4'
            ],
            'time'          : time_list,
            'day'           : dy_list,
            'month'         : mn_list,
            'year'          : yr_list,
    },output_file_lnd)
#-------------------------------------------------------------------------------
