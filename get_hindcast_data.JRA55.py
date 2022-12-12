#!/usr/bin/env python
# Script for downloading ERA5 pressure level and surface data
# links to CDS web interface:
#   https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-pressure-levels
#   https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels
# A list of available variables can also be found in the ERA5 documentation:
#   https://confluence.ecmwf.int/display/CKB/ERA5%3A+data+documentation

import os
import urllib

# single day of data
yr_list,mn_list,dy_list = ['1970'],['01'],['01']

# # list of one date acros many years
# yr_list = [str(y+2000) for y in range(8,19)]
# mn_list = ['10']*len(yr_list)
# dy_list = ['01']*len(yr_list)

# output_path = os.getenv('PWD')+'/data'
output_path = os.getenv('PWD')+'/data_scratch'

for i in range(len(yr_list)):
  yr = yr_list[i]
  mn = mn_list[i]
  dy = dy_list[i]
  time = time_list[i]

  #-----------------------------------------------------------------------------
  dst_file = f'{output_path}/JRA55.atm.{yr}-{mn}-{dy}.nc'
  
  urllib.request.urlretrieve(url, filename=dst_file)

  #-----------------------------------------------------------------------------