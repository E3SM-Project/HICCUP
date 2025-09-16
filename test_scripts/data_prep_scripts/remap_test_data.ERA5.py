#!/usr/bin/env python3
import os, glob, subprocess as sp
pwd = os.getenv('PWD')
# ------------------------------------------------------------------------------
# The user needs to "unpack" the data if downloaded directly from CDS
'''
ncpdq -U test_data/HICCUP_TEST.ERA5.atm.nc test_data/HICCUP_TEST.ERA5.atm.upack.nc
ncpdq -U test_data/HICCUP_TEST.ERA5.sfc.nc test_data/HICCUP_TEST.ERA5.sfc.upack.nc
'''

sfc_src_file_name = 'test_data/HICCUP_TEST.ERA5.sfc.upack.nc'
sfc_dst_file_name = 'test_data/HICCUP_TEST.ERA5.sfc.low-res.nc'

atm_src_file_name = 'test_data/HICCUP_TEST.ERA5.atm.upack.nc'
atm_dst_file_name = 'test_data/HICCUP_TEST.ERA5.atm.low-res.nc'

clean       = True
create_map  = True
regrid_plv  = True
regrid_sfc  = True

# ------------------------------------------------------------------------------
nlat_src,nlon_src = 721,1440
nlat_dst,nlon_dst = 180,360

src_grid_file = f'{pwd}/files_grid/scrip_{nlat_src}x{nlon_src}.nc'
dst_grid_file = f'{pwd}/files_grid/scrip_{nlat_dst}x{nlon_dst}.nc'

all_grid_opts = 'lat_typ=uni#lat_drc=n2s#lon_typ=grn_ctr'
src_grid_opts = f'-G ttl=\'Equi-Angular grid {nlat_src}x{nlon_src}\'#latlon={nlat_src},{nlon_src}#{all_grid_opts}'
dst_grid_opts = f'-G ttl=\'Equi-Angular grid {nlat_dst}x{nlon_dst}\'#latlon={nlat_dst},{nlon_dst}#{all_grid_opts}'

map_file = f'{pwd}/files_map/map_{nlat_src}x{nlon_src}_to_{nlat_dst}x{nlon_dst}.nc'

os.makedirs(f'{pwd}/files_grid', exist_ok=True)  # create files_grid path if it doesn't exist
os.makedirs(f'{pwd}/files_map', exist_ok=True)  # create files_map path if it doesn't exist

# ------------------------------------------------------------------------------
# Set up simple class for coloring terminal text
class tcolor:
    ENDC, BLACK, RED     = '\033[0m','\033[30m','\033[31m'
    GREEN, YELLOW, BLUE  = '\033[32m','\033[33m','\033[34m'
    MAGENTA, CYAN, WHITE = '\033[35m','\033[36m','\033[37m'
# ------------------------------------------------------------------------------
# routine for running shell commands
def run_cmd(cmd,prepend_line=True,use_color=True,shell=True,execute=True):
    """
    Common method for printing and running commands
    """
    prefix='  '
    suffix=''
    if prepend_line : prefix = '\n'+prefix
    msg = f'{prefix}{cmd}{suffix}'
    if use_color : msg = tcolor.GREEN + msg + tcolor.ENDC
    print(msg)
    if not execute : return
    if shell:
        sp.check_call(cmd,shell=True)
    else:
        sp.check_call(cmd.split())
    return
# ------------------------------------------------------------------------------
if create_map:

    # remove old grid and mapping files
    if clean :
        if src_grid_file in glob.glob(src_grid_file) : run_cmd(f'rm {src_grid_file}')
        if dst_grid_file in glob.glob(dst_grid_file) : run_cmd(f'rm {dst_grid_file}')
        if map_file in glob.glob(map_file) : run_cmd(f'rm {map_file}')
    
    # Generate source and target grid files:
    run_cmd(f'ncremap {src_grid_opts} -g {src_grid_file} ')
    run_cmd(f'ncremap {dst_grid_opts} -g {dst_grid_file} ')
    
    # Need to make sure the 'grid_imask' variable is an integer for TempestRemap
    run_cmd(f'ncap2 -s \'grid_imask=int(grid_imask)\' {src_grid_file} {src_grid_file} --ovr')
    run_cmd(f'ncap2 -s \'grid_imask=int(grid_imask)\' {dst_grid_file} {dst_grid_file} --ovr')
    
    # Generate mapping file:
    run_cmd(f'ncremap -a traave --src_grd={src_grid_file} --dst_grd={dst_grid_file} -m {map_file} ')
# ------------------------------------------------------------------------------
if regrid_sfc:
    # remap the 2D "surface" data
    run_cmd(f'ncremap -m {map_file} -i {sfc_src_file_name} -o {sfc_dst_file_name}  ')
    run_cmd(f'ncrename -v lat,latitude -v lon,longitude {sfc_dst_file_name} ')

if regrid_plv:
    # remap the 3D atmosphere data
    run_cmd(f'ncremap -m {map_file} -i {atm_src_file_name} -o {atm_dst_file_name}  ')
    run_cmd(f'ncrename -v lat,latitude -v lon,longitude {atm_dst_file_name} ')
# ------------------------------------------------------------------------------
print()
print(f'sfc data:  {sfc_src_file_name}  >  {sfc_dst_file_name}')
print(f'atm data:  {atm_src_file_name}  >  {atm_dst_file_name}')
print()
# ------------------------------------------------------------------------------
