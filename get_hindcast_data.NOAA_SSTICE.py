#!/usr/bin/env python
#---------------------------------------------------------------------------------------------------
# This script is intended to streamline the aquisition of NOAA SST 
# and sea ice data for hindcasts using HICCUP. The NOAA data is organized
# into yearly files, that are downloaded separately.
# for more information about the dataset:
# https://climatedataguide.ucar.edu/climate-data/sst-data-noaa-high-resolution-025x025-blended-analysis-daily-sst-and-ice-oisstv2
#---------------------------------------------------------------------------------------------------
import ftplib, os, datetime, pandas as pd
#---------------------------------------------------------------------------------------------------
class clr:END,RED,GREEN,MAGENTA,CYAN = '\033[0m','\033[31m','\033[32m','\033[35m','\033[36m'
#---------------------------------------------------------------------------------------------------
usage = f'''
This script is intended to streamline the aquisition of ERA5 pressure level 
and surface data for generating E3SM atmospheric initial conditions with HICCUP

Example of requesting single data file for single initialization date/time:
  {clr.GREEN}python get_hindcast_data.NOAA_SSTICE.py --start-year=<yyyy> --output-root=<path>{clr.END}
'''
from optparse import OptionParser
parser = OptionParser(usage=usage)
parser.add_option('--start-year',  dest='start_year',  default=None,  help='year of first file [yyyymmdd]')
parser.add_option('--final-year',  dest='final_year',  default=None,  help='year of last file [yyyymmdd] (optional)')
parser.add_option('--output-root', dest='output_root', default='./',  help='Output path for data files (default is PWD)')
(opts, args) = parser.parse_args()
#---------------------------------------------------------------------------------------------------
# check that input arguments are valid
if opts.start_year is None: raise ValueError(f'{clr.RED}initialization date was not specified{clr.END}')
if opts.final_year is None: opts.final_year = opts.start_year
#---------------------------------------------------------------------------------------------------
# build list of years from input argument
beg_date = datetime.datetime.strptime(f'{opts.start_year}', '%Y')
end_date = datetime.datetime.strptime(f'{opts.final_year}', '%Y')
yr_list = pd.date_range(beg_date, end_date, freq='YS')
datetime_list = pd.date_range(beg_date, end_date, freq='YS')
#---------------------------------------------------------------------------------------------------
host = 'ftp.cdc.noaa.gov'
path = 'Datasets/noaa.oisst.v2.highres/'
os.makedirs(opts.output_root, exist_ok=True)  # create output_path if it doesn't exist
#---------------------------------------------------------------------------------------------------
# loop over years
for t in datetime_list:
  yr = t.strftime("%Y")

  sst_file_name = f'sst.day.mean.{yr}.nc'
  ice_file_name = f'icec.day.mean.{yr}.nc'

  ftp = ftplib.FTP(host)
  ftp.login()
  ftp.cwd(path)

  for file_name in [sst_file_name, ice_file_name]:
    with open(f'{opts.output_root}/{file_name}', 'wb') as file_pointer:
      print(f'Retrieving file: {file_name}')
      ftp.retrbinary(f'RETR {file_name}', file_pointer.write)
  ftp.quit()
#---------------------------------------------------------------------------------------------------
print('\ndone.')
