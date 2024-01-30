import os, cdsapi, datetime, pandas as pd
from datetime import date, timedelta
server = cdsapi.Client()

get_atm = True
get_sfc = True

bdate_str,edate_str = '2020-01-20 00','2020-01-26 21' # 2024 SCREAM autocalibration (DYAMOND2)

# Create list of date and time 
beg_date = datetime.datetime.strptime(bdate_str, '%Y-%m-%d %H')
end_date = datetime.datetime.strptime(edate_str, '%Y-%m-%d %H')
datetime_list = pd.date_range(beg_date, end_date, freq='3H')

# All available pressure levels
lev = [  '1',  '2',  '3',  '5',  '7', '10', '20', '30', '50', '70','100','125','150','175','200'
      ,'225','250','300','350','400','450','500','550','600','650','700','750','775','800','825'
      ,'850','875','900','925','950','975','1000']

output_path = os.getenv('HOME')+'/HICCUP/data_scratch/nudge_data'

os.makedirs(output_path, exist_ok=True)

# ------------------------------------------------------------------------------
for t in datetime_list:
   yr = t.strftime("%Y")
   mn = t.strftime("%m")
   dy = t.strftime("%d")
   hr = t.strftime("%H")
   # ---------------------------------------------------------------------------
   output_file_plv = f'{output_path}/ERA5.atm.{yr}-{mn}-{dy}-{hr}:00.nc'
   output_file_sfc = f'{output_path}/ERA5.sfc.{yr}-{mn}-{dy}-{hr}:00.nc'
   # ---------------------------------------------------------------------------
   # atmosphere pressure level data
   if get_atm:
      server.retrieve('reanalysis-era5-pressure-levels',{
         'product_type'  : 'reanalysis',
         'pressure_level': lev,
         'time'          : f'{hr}:00',
         'day'           : dy,
         'month'         : mn,
         'year'          : yr,
         'format'        : 'netcdf',
         'variable'      : ['u_component_of_wind', 'v_component_of_wind']
      }, output_file_plv)
   # ---------------------------------------------------------------------------
   # surface data
   if get_sfc:
      server.retrieve('reanalysis-era5-single-levels',{
         'product_type'  : 'reanalysis',
         'time'          : f'{hr}:00',
         'day'           : dy,
         'month'         : mn,
         'year'          : yr,
         'format'        : 'netcdf',
         'variable'      : ['surface_pressure']
      }, output_file_sfc)
   # ----------------------------------------------------------------------------
