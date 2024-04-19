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

Example of requesting multiple data files spanning a range of dates/times every 3-hours:
  {clr.GREEN}python get_validation_data.ERA5.py --start-date=<yyyymmdd> --final-date=<yyyymmdd> --output-root=<path>{clr.END}
'''
from optparse import OptionParser
parser = OptionParser(usage=usage)
parser.add_option('--start-date',  dest='start_date',  default=None,  help='date of first file [yyyymmdd]')
parser.add_option('--final-date',  dest='final_date',  default=None,  help='date of last file [yyyymmdd] (optional)')
# parser.add_option('--start-hour',  dest='start_hour',  default='00',  help='UTC hour of first file (default=00Z)')
# parser.add_option('--final-hour',  dest='final_hour',  default='00',  help='UTC hour of last file (default=00Z)')
# parser.add_option('--data-freq',   dest='data_freq',   default='3h', help='frequency of data files (default=24h)')
parser.add_option('--output-root', dest='output_root', default='./',  help='Output path for data files (default is PWD)')
(opts, args) = parser.parse_args()
#---------------------------------------------------------------------------------------------------
# check that input arguments are valid
if opts.start_date is None: raise ValueError(f'{clr.RED}start date was not specified{clr.END}')
if opts.final_date is None: raise ValueError(f'{clr.RED}final date was not specified{clr.END}')
#---------------------------------------------------------------------------------------------------
# build list of dates and times from input arguments
beg_date = datetime.datetime.strptime(f'{opts.start_date}', '%Y%m%d')
end_date = datetime.datetime.strptime(f'{opts.final_date}', '%Y%m%d')
datetime_list = pd.date_range(beg_date, end_date, freq='24h')
# --------------------------------------------------------------------------------------------------
prs_short_list,prs_era5_list,prs_lev_list = [],[],[]
sfc_short_list,sfc_era5_list = [],[]
def add_var(var_type,var_name_short,name_name_era5,lev=None):
  if var_type=='sfc':
    sfc_short_list.append(var_name_short)
    sfc_era5_list.append(name_name_era5)
  if var_type=='prs':
    if lev is None: raise ValueError(f'{clr.RED}lev needs to be provided for pressure level variables{clr.END}  var: {var_name_short}')
    tmp_lev_list = lev if type(lev) is list else [lev]
    prs_short_list.append(var_name_short)
    prs_era5_list.append(name_name_era5)
    prs_lev_list.append(tmp_lev_list)
# --------------------------------------------------------------------------------------------------

get_atm = True
get_sfc = False

file_prefix = 'ERA5_validation'

hr_mn_list = ['00:00','03:00','06:00','09:00','12:00','15:00','18:00','21:00']

# --------------------------------------------------------------------------------------------------
# buidl list of variables

add_var('prs','Z','geopotential',        lev=['100','500','700'])
add_var('prs','T','temperature',         lev=['500','850'])
add_var('prs','Q','specific_humidity',   lev=['850'])
add_var('prs','U','u_component_of_wind', lev=['200','850'])
add_var('prs','V','v_component_of_wind', lev=['200','850'])

add_var('sfc','TS','skin_temperature')
add_var('sfc','PS','surface_pressure')
# add_var('sfc','SST','sea_surface_temperature')
# add_var('sfc','SIC','sea_ice_cover')
# add_var('sfc','SNOD','snow_depth')

# --------------------------------------------------------------------------------------------------
for t in datetime_list:
  # parse date/time information
  date = t.strftime('%Y-%m-%d')
  yr = t.strftime("%Y")
  mn = t.strftime("%m")
  dy = t.strftime("%d")
  #-----------------------------------------------------------------------------
  # atmosphere pressure level variables
  for v,var in enumerate(prs_short_list):
    output_file = f'{opts.output_root}/{file_prefix}.{var}.{date}.nc'
    server.retrieve('reanalysis-era5-pressure-levels',{
      'product_type'  : 'reanalysis',
      'format'        : 'netcdf',
      'time'          : hr_mn_list,
      'day'           : dy,
      'month'         : mn,
      'year'          : yr,
      'pressure_level': prs_lev_list[v],
      'variable'      : [prs_era5_list[v]],
    }, output_file)

  # single-level variables
  for v,var in enumerate(sfc_short_list):
    output_file = f'{opts.output_root}/{file_prefix}.{var}.{date}.nc'
    server.retrieve('reanalysis-era5-single-levels',{
      'product_type'  : 'reanalysis',
      'format'        : 'netcdf',
      'time'          : hr_mn_list,
      'day'           : dy,
      'month'         : mn,
      'year'          : yr,
      'variable'      : [sfc_era5_list[v]],
    }, output_file)
# --------------------------------------------------------------------------------------------------
