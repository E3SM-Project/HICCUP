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
parser.add_option('--start-hour',  dest='start_hour',  default='00',  help='UTC hour of first file (default=00Z)')
parser.add_option('--final-hour',  dest='final_hour',  default=None,  help='UTC hour of last file (default=00Z)')
parser.add_option('--data-freq',   dest='data_freq',   default='3h',  help='frequency of data files (default=3h)')
parser.add_option('--output-root', dest='output_root', default='./',  help='Output path for data files (default is PWD)')
(opts, args) = parser.parse_args()
#---------------------------------------------------------------------------------------------------
# check that input arguments are valid
if opts.start_date is None: raise ValueError(f'{clr.RED}start date was not specified{clr.END}')
# if opts.final_date is None: raise ValueError(f'{clr.RED}final date was not specified{clr.END}')
if opts.final_date is None: opts.final_date = opts.start_date
if opts.final_date is None: opts.final_date = opts.start_date
if opts.final_hour is None: opts.final_hour = opts.start_hour
#---------------------------------------------------------------------------------------------------
# build list of dates and times from input arguments
# beg_date = datetime.datetime.strptime(f'{opts.start_date}', '%Y%m%d')
# end_date = datetime.datetime.strptime(f'{opts.final_date}', '%Y%m%d')
beg_date = datetime.datetime.strptime(f'{opts.start_date} {opts.start_hour}', '%Y%m%d %H')
end_date = datetime.datetime.strptime(f'{opts.final_date} {opts.final_hour}', '%Y%m%d %H')
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

# --------------------------------------------------------------------------------------------------
# def generate_hour_list(hr_freq):
#     return pd.date_range('00:00', '23:59', freq=opts.data_freq).strftime('%H:%M').tolist()
# --------------------------------------------------------------------------------------------------
# hr_mn_list = ['00:00','03:00','06:00','09:00','12:00','15:00','18:00','21:00']
hr_mn_list = pd.date_range('00:00', '23:59', freq=opts.data_freq).strftime('%H:%M').tolist()
#---------------------------------------------------------------------------------------------------
lev_all = [   '1',  '2',  '3',  '5',  '7', '10', '20', '30', '50', '70','100','125',
            '150','175','200','225','250','300','350','400','450','500','550','600',
            '650','700','750','775','800','825','850','875','900','925','950','975','1000']
# --------------------------------------------------------------------------------------------------
# buidl list of variables

if True:
    add_var('prs','Z','geopotential',        lev=lev_all)
    add_var('prs','T','temperature',         lev=lev_all)
    add_var('prs','Q','specific_humidity',   lev=lev_all)
    add_var('prs','U','u_component_of_wind', lev=lev_all)
    add_var('prs','V','v_component_of_wind', lev=lev_all)
else:
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
    #---------------------------------------------------------------------------
    # parse date/time information
    yr = t.strftime("%Y")
    mn = t.strftime("%m")
    dy = t.strftime("%d")
    hr = t.strftime("%H")
    # hr_min = f'{hr}:00'
    #---------------------------------------------------------------------------
    # Specify output file names
    output_file_plv = f'{opts.output_root}/{file_prefix}.atm.{yr}-{mn}-{dy}.nc'
    # output_file_plv = f'{opts.output_root}/{file_prefix}.atm.{yr}-{mn}-{dy}.{hr}.nc'
    output_file_sfc = output_file_plv.replace('.atm.','.sfc.')
    #---------------------------------------------------------------------------
    print()
    print(f'  output_file_plv: {output_file_plv}')
    print(f'  output_file_sfc: {output_file_sfc}')
    print()
    # exit()
    #---------------------------------------------------------------------------
    # atmosphere pressure level variables
    for v,var in enumerate(prs_short_list):
        output_file_plv_tmp = output_file_plv.replace('.nc',f'.{var}.nc')
        server.retrieve('reanalysis-era5-pressure-levels',{
            'product_type'  : 'reanalysis',
            'format'        : 'netcdf',
            'time'          : hr_mn_list,
            # 'time'          : hr_min,
            'day'           : dy,
            'month'         : mn,
            'year'          : yr,
            'pressure_level': prs_lev_list[v],
            'variable'      : [prs_era5_list[v]],
        }, output_file_plv_tmp)

    # single-level variables
    for v,var in enumerate(sfc_short_list):
        output_file_sfc_tmp = output_file_sfc.replace('.nc',f'.{var}.nc')
        server.retrieve('reanalysis-era5-single-levels',{
            'product_type'  : 'reanalysis',
            'format'        : 'netcdf',
            'time'          : hr_mn_list,
            # 'time'          : hr_min,
            'day'           : dy,
            'month'         : mn,
            'year'          : yr,
            'variable'      : [sfc_era5_list[v]],
        }, output_file_sfc_tmp)
# --------------------------------------------------------------------------------------------------
