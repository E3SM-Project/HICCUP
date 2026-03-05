#!/usr/bin/env python
#---------------------------------------------------------------------------------------------------
# links to CDS web interface:
#   https://cds.climate.copernicus.eu/datasets/reanalysis-era5-pressure-levels
#   https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels
# A list of available variables can also be found in the ERA5 documentation:
#   https://confluence.ecmwf.int/display/CKB/ERA5%3A+data+documentation
# NOTE: requires cdsapi >= 0.7.0 and updated ~/.cdsapirc with personal access token
#---------------------------------------------------------------------------------------------------
import os, cdsapi, datetime, pandas as pd
server = cdsapi.Client()
#---------------------------------------------------------------------------------------------------
class clr:END,RED,GREEN,MAGENTA,CYAN = '\033[0m','\033[31m','\033[32m','\033[35m','\033[36m'
#---------------------------------------------------------------------------------------------------
usage = f'''
This script is intended to streamline the aquisition of ERA5 pressure level 
and surface data for processing data for E3SM atmospheric nudging with HICCUP

Minimum example of requesting single data file for single initialization date/time:
  {clr.GREEN}python get_nudging_data.ERA5.py --start-date=<yyyymmdd>{clr.END}

Example of requesting multiple data files spanning a range of dates/times every 3-hours:
  {clr.GREEN}python get_nudging_data.ERA5.py --start-date=<yyyymmdd> --final-date=<yyyymmdd> --start-hour=<hh> --final-hour=<hh> --data-freq=3h --output-root=<path>{clr.END}
'''
from optparse import OptionParser
parser = OptionParser(usage=usage)
parser.add_option('--start-date',  dest='start_date',  default=None,   help='date of first file [yyyymmdd]')
parser.add_option('--final-date',  dest='final_date',  default=None,   help='date of last file [yyyymmdd] (optional)')
parser.add_option('--start-hour',  dest='start_hour',  default='00',   help='UTC hour of first file (default=00Z)')
parser.add_option('--final-hour',  dest='final_hour',  default='00',   help='UTC hour of last file (default=00Z)')
parser.add_option('--data-freq',   dest='data_freq',   default='3h',   help='frequency of data files (default=3h)')
parser.add_option('--batch-mode',  dest='batch_mode',  default='month',help='batch requests by: month (default) or day')
parser.add_option('--output-root', dest='output_root', default='./',   help='Output path for data files (default is PWD)')
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
lev = [  '1',  '2',  '3',  '5',  '7', '10', '20', '30', '50', '70','100','125',
       '150','175','200','225','250','300','350','400','450','500','550','600',
       '650','700','750','775','800','825','850','875','900','925','950','975','1000']
#---------------------------------------------------------------------------------------------------
# Group timestamps into batches to minimize API calls
def get_batch_key(t):
    if opts.batch_mode == 'day':   return (t.strftime("%Y"), t.strftime("%m"), t.strftime("%d"))
    if opts.batch_mode == 'month': return (t.strftime("%Y"), t.strftime("%m"))
    raise ValueError(f'{clr.RED}invalid batch_mode: {opts.batch_mode}{clr.END}')

batches = {}
for t in datetime_list:
    batches.setdefault(get_batch_key(t), []).append(t)
#---------------------------------------------------------------------------------------------------
get_atm = True
get_sfc = True

# loop over batches - set by --batch-mode argument
for key, times in batches.items():
    yr = times[0].strftime("%Y")
    mn = times[0].strftime("%m")

    days  = sorted(set(t.strftime("%d")    for t in times))
    hours = sorted(set(t.strftime("%H:%M") for t in times))

    tag = f'{yr}-{mn}' if opts.batch_mode == 'month' else f'{yr}-{mn}-{times[0].strftime("%d")}'
    output_file_plv = f'{opts.output_root}/ERA5.atm.{tag}.nc'
    output_file_sfc = output_file_plv.replace('.atm.', '.sfc.')

    print()
    if get_atm: print(f'  output_file_plv : {output_file_plv}')
    if get_sfc: print(f'  output_file_sfc : {output_file_sfc}')
    print(f'  days  : {days}')
    print(f'  hours : {hours}')
    print()

    if get_atm and not os.path.exists(output_file_plv):
        server.retrieve('reanalysis-era5-pressure-levels', {
            'product_type'  : ['reanalysis'],
            'data_format'   : 'netcdf',
            'pressure_level': lev,
            'year'          : [yr],
            'month'         : [mn],
            'day'           : days,
            'time'          : hours,
            'variable'      : ['u_component_of_wind', 'v_component_of_wind'],
        }, output_file_plv)

    if get_sfc and not os.path.exists(output_file_sfc):
        server.retrieve('reanalysis-era5-single-levels', {
            'product_type'  : ['reanalysis'],
            'data_format'   : 'netcdf',
            'year'          : [yr],
            'month'         : [mn],
            'day'           : days,
            'time'          : hours,
            'variable'      : ['surface_pressure'],
        }, output_file_sfc)
