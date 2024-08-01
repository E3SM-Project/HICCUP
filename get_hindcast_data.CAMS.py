#!/usr/bin/env python
#---------------------------------------------------------------------------------------------------
# See the CAMS data documentation for more information:
#   https://confluence.ecmwf.int/display/CKB/CAMS%3A+Reanalysis+data+documentation
#---------------------------------------------------------------------------------------------------
import os, yaml, cdsapi, datetime, pandas as pd
#---------------------------------------------------------------------------------------------------
class clr:END,RED,GREEN,MAGENTA,CYAN = '\033[0m','\033[31m','\033[32m','\033[35m','\033[36m'
#---------------------------------------------------------------------------------------------------
usage = f'''
This script is intended to streamline the aquisition of CAMS atmospheric composition 
data for generating E3SM atmospheric initial conditions with HICCUP

Minimum example of requesting single data file for single initialization date/time:
  {clr.GREEN}python get_hindcast_data.CAMS.py --start-date=<yyyymmdd>{clr.END}

Example of requesting multiple data files spanning a range of dates/times every 3-hours:
  {clr.GREEN}python get_hindcast_data.CAMS.py --start-date=<yyyymmdd> --final-date=<yyyymmdd> --start-hour=<hh> --final-hour=<hh> --data-freq=3h --output-root=<path>{clr.END}
'''
from optparse import OptionParser
parser = OptionParser(usage=usage)
parser.add_option('--start-date',  dest='start_date',  default=None,  help='date of first file [yyyymmdd]')
parser.add_option('--final-date',  dest='final_date',  default=None,  help='date of last file [yyyymmdd] (optional)')
parser.add_option('--start-hour',  dest='start_hour',  default='00',  help='UTC hour of first file (default=00Z)')
parser.add_option('--final-hour',  dest='final_hour',  default='00',  help='UTC hour of last file (default=00Z)')
parser.add_option('--data-freq',   dest='data_freq',   default='24h', help='frequency of data files (default=24h)')
parser.add_option('--output-root', dest='output_root', default='./',  help='Output path for data files (default is PWD)')
(opts, args) = parser.parse_args()
#---------------------------------------------------------------------------------------------------
# check that input arguments are valid
if opts.start_date is None: raise ValueError(f'{clr.RED}initialization date was not specified{clr.END}')
if opts.final_date is None: opts.final_date = opts.start_date
#---------------------------------------------------------------------------------------------------
# build list of dates and times from input arguments
beg_date = datetime.datetime.strptime(f'{opts.start_date} {opts.start_hour}', '%Y%m%d %H')
end_date = datetime.datetime.strptime(f'{opts.final_date} {opts.final_hour}', '%Y%m%d %H')
datetime_list = pd.date_range(beg_date, end_date, freq=opts.data_freq)
#---------------------------------------------------------------------------------------------------
# Request all available pressure levels
lev = [  '1',  '2',  '3',  '5',  '7', '10', '20', '30', '50', '70','100',
       '200','300','400','500','600','700','800','850','900','925','950','1000']
#---------------------------------------------------------------------------------------------------
# loop over dates

get_atm = True
get_sfc = True

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
  output_file_atm = f'{opts.output_root}/CAMS.atm.{yr}-{mn}-{dy}.{hr}.nc'
  output_file_sfc = f'{opts.output_root}/CAMS.sfc.{yr}-{mn}-{dy}.{hr}.nc'
  #-----------------------------------------------------------------------------
  # atmossphere pressure level data
  if get_atm:
      with open(os.getenv('HOME')+'/.adsapirc', 'r') as f: credentials = yaml.safe_load(f)
      ads_server = cdsapi.Client(url=credentials['url'], key=credentials['key'])
      ads_server.retrieve('cams-global-reanalysis-eac4',{
        'type'          : 'an',
        'stream'        : 'oper',
        # 'levtype'       : 'pl',
        'format'        : 'netcdf',
        'pressure_level': lev,
        'time'          : hr_min,
        'date': f'{yr}-{mn}-{dy}/{yr}-{mn}-{dy}',
        'variable'      : ['dust_aerosol_0.03-0.55um_mixing_ratio'
                          ,'dust_aerosol_0.55-0.9um_mixing_ratio'
                          ,'dust_aerosol_0.9-20um_mixing_ratio'
                          ],
      }, output_file_atm)

  #-----------------------------------------------------------------------------
# surface data
if get_sfc:
    cds_server = cdsapi.Client()
    cds_server.retrieve('reanalysis-era5-single-levels',{
        'product_type'  : 'reanalysis',
        'grid'          : [0.75, 0.75],
        'time'          : time,
        'day'           : dy,
        'month'         : mn,
        'year'          : yr,
        'format'        : 'netcdf',
        'variable'      : ['surface_pressure'],
    }, output_file_sfc)
    del cds_server
  #-----------------------------------------------------------------------------
