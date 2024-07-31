#!/usr/bin/env python
# Script for downloading NOAA SST and sea ice data for hindcasts
# Only one year is downloaded at a time
# for more information about the dataset:
# https://climatedataguide.ucar.edu/climate-data/sst-data-noaa-high-resolution-025x025-blended-analysis-daily-sst-and-ice-oisstv2

import ftplib
import os
import shutil

yr_list = [2008]

host = 'ftp.cdc.noaa.gov'
path = 'Datasets/noaa.oisst.v2.highres/'

output_path = os.getenv('HOME')+'/HICCUP/test_data'
os.makedirs(output_path, exist_ok=True)  # create output_path if it doesn't exist

for year in yr_list:

  sst_file_name = f'sst.day.mean.{year}.nc'
  ice_file_name = f'icec.day.mean.{year}.nc'

  ftp = ftplib.FTP(host)
  ftp.login()
  ftp.cwd(path)

  for file_name in [sst_file_name, ice_file_name]:
    with open(f'{output_path}/{file_name}', 'wb') as file_pointer:
      print(f'Retrieving file: {file_name}')
      ftp.retrbinary(f'RETR {file_name}', file_pointer.write)

  ftp.quit()

# This renaming will happen in test_scripts/remap.test_data.NOAA.py 
# shutil.move(f'{output_path}/sst.day.mean.{year}.nc', f'{output_path}/HICCUP_TEST.NOAA.sst.nc')
# shutil.move(f'{output_path}/icec.day.mean.{year}.nc', f'{output_path}/HICCUP_TEST.NOAA.ice.nc')

print('\ndone.')
