#!/usr/bin/env python
# Script for downloading ERA5 pressure level and surface data
# links to CDS web interface:
#   https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-pressure-levels
#   https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels
# A list of available variables can also be found in the ERA5 documentation:
#   https://confluence.ecmwf.int/display/CKB/ERA5%3A+data+documentation
#-----------------------------------------------------------------------------
import os
import cdsapi
server = cdsapi.Client()
#-----------------------------------------------------------------------------
get_plv = True
get_mlv = True
get_sfc = True
get_lnd = False
#-----------------------------------------------------------------------------
# obtain single date for testing
yr,mn,dy,time = '2008','10','01','00:00'

# Small level set for testing
plev = [ '50','100','150','200','300','400','500','600'
       ,'700','750','800','850','900','950','1000']

#-----------------------------------------------------------------------------
# Specify output file names
output_path = os.getenv('PWD')+'/test_data'
output_file_plv = output_path+f'/HICCUP_TEST.ERA5.atm.nc'
output_file_mlv = output_path+f'/HICCUP_TEST.ERA5.mlv.nc'
output_file_sfc = output_path+f'/HICCUP_TEST.ERA5.sfc.nc'
output_file_lnd = output_path+f'/HICCUP_TEST.ERA5.lnd.nc'
#-------------------------------------------------------------------------------
# atmossphere pressure level data
if get_plv:
  server.retrieve('reanalysis-era5-pressure-levels',{
      'product_type'  : 'reanalysis',
      'pressure_level': plev,
      'time'          : time,
      'day'           : dy,
      'month'         : mn,
      'year'          : yr,
      'format'        : 'netcdf',
      'variable'      : [ 'temperature',
                          'specific_humidity',
                          'geopotential',
                          'u_component_of_wind',
                          'v_component_of_wind',
                          'ozone_mass_mixing_ratio',
                          'specific_cloud_ice_water_content',
                          'specific_cloud_liquid_water_content',
                        ],
  }, output_file_plv)
#-------------------------------------------------------------------------------
# model level data
if get_mlv:
    server.retrieve('reanalysis-era5-complete',{
        'product_type'  : 'reanalysis',
        # 'pressure_level': mlev,
        'levtype': 'ml', 
        'time'          : time,
        'day'           : dy,
        'month'         : mn,
        'year'          : yr,
        'format'        : 'netcdf',
        'variable'      : [ 'temperature',
                            'specific_humidity',
                            'geopotential',
                            'u_component_of_wind',
                            'v_component_of_wind',
                            'ozone_mass_mixing_ratio',
                            'specific_cloud_ice_water_content',
                            'specific_cloud_liquid_water_content',
                          ],
    }, output_file_mlv)
#-------------------------------------------------------------------------------
# surface data
if get_sfc:
  server.retrieve('reanalysis-era5-single-levels',{
      'product_type'  : 'reanalysis',
      # 'time'          : time_list,
      # 'day'           : dy_list,
      # 'month'         : mn_list,
      # 'year'          : yr_list,
      'time'          : time,
      'day'           : dy,
      'month'         : mn,
      'year'          : yr,
      'format'        : 'netcdf',
      'variable'      : [ 'surface_pressure',
                          'skin_temperature',
                          'sea_surface_temperature',
                          '2m_temperature',
                          'soil_temperature_level_1',
                          'soil_temperature_level_2',
                          'soil_temperature_level_3',
                          'soil_temperature_level_4',
                          'snow_depth',
                          'temperature_of_snow_layer',
                          'geopotential',
                          'sea_ice_cover',
                          # 'leaf_area_index_high_vegetation',
                          # 'leaf_area_index_low_vegetation',
                          # 'skin_reservoir_content',
                          # 'snow_albedo',
                          # 'snow_density',
                          # 'volumetric_soil_water_layer_1',
                          # 'volumetric_soil_water_layer_2',
                          # 'volumetric_soil_water_layer_3',
                          # 'volumetric_soil_water_layer_4',
                          # 'orography',
                        ],
  }, output_file_sfc)
#-------------------------------------------------------------------------------
# land model data
if get_lnd:
  server.retrieve('reanalysis-era5-land',{
          'format': 'netcdf',
          'time'          : time,
          'day'           : dy,
          'month'         : mn,
          'year'          : yr,
          'variable': [
                        'leaf_area_index_high_vegetation',
                        'leaf_area_index_low_vegetation',
                        'skin_reservoir_content',
                        'snow_albedo',
                        'snow_cover',
                        'snow_density',
                        'snow_depth',
                        'snow_depth_water_equivalent',
                        'soil_temperature_level_1',
                        'soil_temperature_level_2',
                        'soil_temperature_level_3',
                        'soil_temperature_level_4',
                        'temperature_of_snow_layer',
                        'volumetric_soil_water_layer_1',
                        'volumetric_soil_water_layer_2',
                        'volumetric_soil_water_layer_3',
                        'volumetric_soil_water_layer_4',
                      ],
          
          
  },output_file_lnd)
#-------------------------------------------------------------------------------
