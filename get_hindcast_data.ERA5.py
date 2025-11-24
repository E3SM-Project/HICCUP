#!/usr/bin/env python
#---------------------------------------------------------------------------------------------------
# links to CDS web interface:
#   https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-pressure-levels
#   https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels
# A list of available variables can also be found in the ERA5 documentation:
#   https://confluence.ecmwf.int/display/CKB/ERA5%3A+data+documentation
#---------------------------------------------------------------------------------------------------
import os, cdsapi, datetime, pandas as pd
server = cdsapi.Client()
#---------------------------------------------------------------------------------------------------
class clr:END,RED,GREEN,MAGENTA,CYAN = '\033[0m','\033[31m','\033[32m','\033[35m','\033[36m'
#---------------------------------------------------------------------------------------------------
usage = f'''
This script is intended to streamline the aquisition of ERA5 pressure level 
and surface data for generating E3SM atmospheric initial conditions with HICCUP

Minimum example of requesting single data file for single initialization date/time:
  {clr.GREEN}python get_hindcast_data.ERA5.py --start-date=<yyyymmdd>{clr.END}

Example of requesting multiple data files spanning a range of dates/times every 3-hours:
  {clr.GREEN}python get_hindcast_data.ERA5.py --start-date=<yyyymmdd> --final-date=<yyyymmdd> --start-hour=<hh> --final-hour=<hh> --data-freq=3h --output-root=<path>{clr.END}
'''
from optparse import OptionParser
parser = OptionParser(usage=usage)
parser.add_option('--start-date',  dest='start_date',  default=None,  help='date of first file [yyyymmdd]')
parser.add_option('--final-date',  dest='final_date',  default=None,  help='date of last file [yyyymmdd] (optional)')
parser.add_option('--start-hour',  dest='start_hour',  default='00',  help='UTC hour of first file (default=00Z)')
parser.add_option('--final-hour',  dest='final_hour',  default=None,  help='UTC hour of last file (default=00Z)')
parser.add_option('--data-freq',   dest='data_freq',   default='24h', help='frequency of data files (default=24h)')
parser.add_option('--output-root', dest='output_root', default='./',  help='Output path for data files (default is PWD)')
(opts, args) = parser.parse_args()
#---------------------------------------------------------------------------------------------------
# check that input arguments are valid
if opts.start_date is None: raise ValueError(f'{clr.RED}initialization date was not specified{clr.END}')
if opts.final_date is None: opts.final_date = opts.start_date
if opts.final_hour is None: opts.final_hour = opts.start_hour
#---------------------------------------------------------------------------------------------------
# build list of dates and times from input arguments
beg_date = datetime.datetime.strptime(f'{opts.start_date} {opts.start_hour}', '%Y%m%d %H')
end_date = datetime.datetime.strptime(f'{opts.final_date} {opts.final_hour}', '%Y%m%d %H')
datetime_list = pd.date_range(beg_date, end_date, freq=opts.data_freq)
#---------------------------------------------------------------------------------------------------
# Request all available pressure levels
lev = [  '1',  '2',  '3',  '5',  '7', '10', '20', '30', '50', '70','100','125',
       '150','175','200','225','250','300','350','400','450','500','550','600',
       '650','700','750','775','800','825','850','875','900','925','950','975','1000']
#---------------------------------------------------------------------------------------------------
# loop over dates

get_atm = True
get_sfc = True
get_lnd = False

for t in datetime_list:
  #-----------------------------------------------------------------------------
  # parse date/time information
  yr = t.strftime("%Y")
  mn = t.strftime("%m")
  dy = t.strftime("%d")
  hr = t.strftime("%H")
  hr_min = f'{hr}:00'
  #-----------------------------------------------------------------------------
  # Specify output file names
  output_file_plv = f'{opts.output_root}/ERA5.atm.{yr}-{mn}-{dy}.{hr}.nc'
  output_file_mlv = output_file_plv.replace('.atm.','.mlev.')
  output_file_sfc = output_file_plv.replace('.atm.','.sfc.')
  output_file_lnd = output_file_plv.replace('.atm.','.lnd.')
  #-----------------------------------------------------------------------------
  print()
  if get_atm: print(f'output_file_plv: {output_file_plv}')
  # if get_atm: print(f'output_file_mlv: {output_file_mlv}')
  if get_sfc: print(f'output_file_sfc: {output_file_sfc}')
  if get_lnd: print(f'output_file_lnd: {output_file_lnd}')
  print()
  #-----------------------------------------------------------------------------
  # atmossphere pressure level data
  if get_atm:
      server.retrieve('reanalysis-era5-pressure-levels',{
        'product_type'  : 'reanalysis',
        'format'        : 'netcdf',
        'pressure_level': lev,
        'time'          : hr_min,
        'day'           : dy,
        'month'         : mn,
        'year'          : yr,
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
  #-----------------------------------------------------------------------------
  # surface data
  if get_sfc:
      server.retrieve('reanalysis-era5-single-levels',{
        'product_type'  : 'reanalysis',
        'format'        : 'netcdf',
        'time'          : hr_min,
        'day'           : dy,
        'month'         : mn,
        'year'          : yr,
        'variable'      : ['surface_pressure'
                          ,'skin_temperature'
                          ,'sea_surface_temperature'
                          ,'2m_temperature'
                          ,'soil_temperature_level_1'
                          ,'soil_temperature_level_2'
                          ,'soil_temperature_level_3'
                          ,'soil_temperature_level_4'
                          ,'snow_depth'
                          ,'temperature_of_snow_layer'
                          ,'geopotential'
                          ,'sea_ice_cover'
                          ],
      }, output_file_sfc)
  #-----------------------------------------------------------------------------
  # land model data
  if get_lnd:
      server.retrieve('reanalysis-era5-land',{
        'product_type'  : 'reanalysis',
        'format'        : 'netcdf',
        'time'          : hr_min,
        'day'           : dy,
        'month'         : mn,
        'year'          : yr,
        'variable'      : ['leaf_area_index_high_vegetation'
                          ,'leaf_area_index_low_vegetation'
                          ,'skin_reservoir_content'
                          ,'snow_albedo'
                          ,'snow_cover'
                          ,'snow_density'
                          ,'snow_depth'
                          ,'snow_depth_water_equivalent'
                          ,'soil_temperature_level_1'
                          ,'soil_temperature_level_2'
                          ,'soil_temperature_level_3'
                          ,'soil_temperature_level_4'
                          ,'temperature_of_snow_layer'
                          ,'volumetric_soil_water_layer_1'
                          ,'volumetric_soil_water_layer_2'
                          ,'volumetric_soil_water_layer_3'
                          ,'volumetric_soil_water_layer_4'
                          ],
      },output_file_lnd)
  #-----------------------------------------------------------------------------
