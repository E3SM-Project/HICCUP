#!/usr/bin/env python
import os, glob, subprocess as sp, datetime, pandas as pd
#---------------------------------------------------------------------------------------------------
class clr:END,RED,GREEN,MAGENTA,CYAN = '\033[0m','\033[31m','\033[32m','\033[35m','\033[36m'
def run_cmd(cmd): print('\n'+clr.GREEN+cmd+clr.END) ; os.system(cmd); return
#---------------------------------------------------------------------------------------------------
usage = f'''
This script is intended to facilitate regridding ERA5 data for hindcast
validation to the native pg2 output grid of E3SM atmosphere data

Example of multiple data files spanning a range of dates/times every 3-hours:
  {clr.GREEN}python remap.validation_data_native.ERA5.py --start-date=<yyyymmdd> --final-date=<yyyymmdd> --start-hour=<hh> --final-hour=<hh> --data-freq=3h --scratch-root=<path> --input-root=<path> --output-root=<path>{clr.END}
'''
from optparse import OptionParser
parser = OptionParser(usage=usage)
parser.add_option('--output-grid-ne',dest='output_grid_ne',default=30,    help='number of elements (ne) for output grid')
parser.add_option('--start-date',    dest='start_date',    default=None,  help='date of first file [yyyymmdd]')
parser.add_option('--final-date',    dest='final_date',    default=None,  help='date of last file [yyyymmdd] (optional)')
parser.add_option('--start-hour',    dest='start_hour',    default='00',  help='UTC hour of first file (default=00Z)')
parser.add_option('--final-hour',    dest='final_hour',    default='00',  help='UTC hour of last file (default=00Z)')
parser.add_option('--data-freq',     dest='data_freq',     default='24h', help='frequency of data files (default=24h)')
parser.add_option('--input-root',    dest='input_root',    default=None,  help='Input path for data files')
parser.add_option('--output-root',   dest='output_root',   default=None,  help='Output path for data files')
parser.add_option('--scratch-root',  dest='scratch_root',  default=None,  help='Root path for grid and map files')
parser.add_option('--skip-map', dest='skip_map', action='store_true', default=False, help='Switch to disable grid and map creation')
(opts, args) = parser.parse_args()
#---------------------------------------------------------------------------------------------------
# check that input arguments are valid
if opts.start_date   is None: raise ValueError(f'{clr.RED}initialization date was not specified{clr.END}')
if opts.final_date   is None: opts.final_date = opts.start_date
if opts.input_root   is None: raise ValueError(f'{clr.RED}input root was not specified{clr.END}')
if opts.output_root  is None: raise ValueError(f'{clr.RED}output root was not specified{clr.END}')
if opts.scratch_root is None: raise ValueError(f'{clr.RED}scratch root was not specified{clr.END}')
if not os.path.exists(opts.input_root) : raise ValueError(f'{clr.RED}input root does not exist{clr.END}: {opts.input_root}') 
if not os.path.exists(opts.output_root): raise ValueError(f'{clr.RED}output root does not exist{clr.END}: {opts.output_root}') 
#---------------------------------------------------------------------------------------------------
# build list of dates and times from input arguments
beg_date = datetime.datetime.strptime(f'{opts.start_date} {opts.start_hour}', '%Y%m%d %H')
end_date = datetime.datetime.strptime(f'{opts.final_date} {opts.final_hour}', '%Y%m%d %H')
datetime_list = pd.date_range(beg_date, end_date, freq=opts.data_freq)
# --------------------------------------------------------------------------------------------------

var_list = ['Z','T','Q','U','V','TS','PS']

nlat_src,nlon_src = 721,1440

src_grid_name = f'{nlat_src}x{nlon_src}_n2s'
dst_grid_name = f'ne{opts.output_grid_ne}pg2'

grid_file_root = f'{opts.scratch_root}/files_grid'
map_file_root  = f'{opts.scratch_root}/files_map'

src_grid_file = f'{grid_file_root}/scrip_{src_grid_name}.nc'
dst_grid_file = f'{grid_file_root}/scrip_{dst_grid_name}.nc'

map_file = f'{map_file_root}/map_{src_grid_name}_to_{dst_grid_name}.nc'

regrid_data = True
create_map  = True
unpack      = True
clean       = True

if opts.skip_map: create_map  = False

# ------------------------------------------------------------------------------
# print some informative stuff
print(f'''
  start_date   : {opts.start_date}
  final_date   : {opts.final_date}
  start_hour   : {opts.start_hour}
  final_hour   : {opts.final_hour}
  data_freq    : {opts.data_freq}

  src_grid_file: {src_grid_file}
  dst_grid_file: {dst_grid_file}
  map_file     : {map_file}

  regrid_data  : {regrid_data}
  create_map   : {create_map}
  unpack       : {unpack}
  clean        : {clean}
''')
# ------------------------------------------------------------------------------
if create_map:

   if not os.path.exists(grid_file_root): os.mkdir(grid_file_root)
   if not os.path.exists(map_file_root): os.mkdir(map_file_root)

   # remove old grid and mapping files
   if clean :
      if os.path.exists(src_grid_file) : run_cmd(f'rm {src_grid_file}')
      if os.path.exists(dst_grid_file) : run_cmd(f'rm {dst_grid_file}')
      if os.path.exists(map_file)      : run_cmd(f'rm {map_file}')
   
   # Generate source grid file:
   cmd  = f'ncremap -G '
   cmd += f'ttl=\'Equi-Angular grid {nlat_src}x{nlon_src}\''
   cmd += f'#latlon={nlat_src},{nlon_src}'
   cmd += f'#lat_typ=uni'
   cmd += f'#lat_drc=n2s'
   cmd += f'#lon_typ=grn_ctr '
   cmd += f'-g {src_grid_file} '
   run_cmd(cmd)

   # Generate target grid file:
   ne = opts.output_grid_ne
   run_cmd(f'GenerateCSMesh --alt --res {ne} --file {grid_file_root}/exodus_ne{ne}.g')
   run_cmd(f'GenerateVolumetricMesh --in {grid_file_root}/exodus_ne{ne}.g --out {grid_file_root}/exodus_ne{ne}pg2.g --np 2 --uniform')
   run_cmd(f'ConvertMeshToSCRIP --in {grid_file_root}/exodus_ne{ne}pg2.g --out {grid_file_root}/scrip_ne{ne}pg2.nc')

   # Generate mapping file:
   cmd  = f'ncremap -6 --alg_typ=traave '
   cmd += f' --src_grd={src_grid_file}'
   cmd += f' --dst_grd={dst_grid_file}'
   cmd += f' -m {map_file} '
   run_cmd(cmd)

# --------------------------------------------------------------------------------------------------
# remap the data
for var in var_list:
   for t in datetime_list:
      # parse date/time information
      yr = t.strftime("%Y")
      mn = t.strftime("%m")
      dy = t.strftime("%d")
      hr = t.strftime("%H")
      hr_min = f'{hr}:00'
      
   src_file_name = f'{opts.input_root}/ERA5_validation.{var}.{date}.nc'
   dst_file_name = src_file_name.replace('.nc',f'.remap_{dst_grid_name}.nc')
   
   # --------------------------------------------------------------------------
   if unpack: run_cmd(f'ncpdq --ovr -U {src_file_name} {src_file_name}')
   # --------------------------------------------------------------------------
   # remap the data
   if regrid_data:
      run_cmd(f'ncremap {alg} -m {map_file} -i {src_file_name} -o {dst_file_name}  ')

# --------------------------------------------------------------------------------------------------
