#!/usr/bin/env python
import os
import glob
import subprocess as sp
pwd = os.getenv('PWD')

class tcolor:
    ENDC     = '\033[0m'
    BLACK    = '\033[30m'
    RED      = '\033[31m'
    GREEN    = '\033[32m'
    YELLOW   = '\033[33m'
    BLUE     = '\033[34m'
    MAGENTA  = '\033[35m'
    CYAN     = '\033[36m'
    WHITE    = '\033[37m'

# The user likely needs to "unpack" the data if downloaded directly from CDS
# examples: 
# ncpdq -U HICCUP_TEST.ERA5.atm.nc HICCUP_TEST.ERA5.atm.upack.nc
# ncpdq -U HICCUP_TEST.ERA5.sfc.nc HICCUP_TEST.ERA5.sfc.upack.nc

src_file_name = 'data/HICCUP_TEST.ERA5.sfc.upack.nc'
dst_file_name = src_file_name.replace('.upack.nc','.low-res.nc')

nlat_src,nlon_src = 721,1440
nlat_dst,nlon_dst = 180,360

src_grid_file = f'{pwd}/files_grid/scrip_{nlat_src}x{nlon_src}.nc'
dst_grid_file = f'{pwd}/files_grid/scrip_{nlat_dst}x{nlon_dst}.nc'

map_file = f'{pwd}/map_files/map_{nlat_src}x{nlon_src}_to_{nlat_dst}x{nlon_dst}.nc'

alg = '-a tempest'

clean       = True
create_map  = True
regrid_data = True

# ------------------------------------------------------------------------------
def run_cmd(cmd,prepend_line=True,use_color=True,shell=True,execute=False):
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
    
    # Generate source grid file:
    cmd  = f'ncremap {alg} -G ttl=\'Equi-Angular grid {nlat_src}x{nlon_src}\'#latlon={nlat_src},{nlon_src}'+\
           f'#lat_typ=uni#lat_drc=n2s#lon_typ=grn_ctr -g {src_grid_file} '
    run_cmd(cmd)
    
    # Generate target grid file:
    cmd  = f'ncremap {alg} -G ttl=\'Equi-Angular grid {nlat_dst}x{nlon_dst}\'#latlon={nlat_dst},{nlon_dst}'+\
           f'#lat_typ=uni#lat_drc=n2s#lon_typ=grn_ctr -g {dst_grid_file} '
    run_cmd(cmd)
    
    # Need to make sure the 'grid_imask' variable is an integer for TempestRemap
    cmd = f'ncap2 -s \'grid_imask=int(grid_imask)\' {src_grid_file} {src_grid_file} --ovr'
    run_cmd(cmd)

    cmd = f'ncap2 -s \'grid_imask=int(grid_imask)\' {dst_grid_file} {dst_grid_file} --ovr'
    run_cmd(cmd)
    
    # Generate mapping file:
    cmd  = f'ncremap {alg} -a fv2fv --src_grd={src_grid_file} --dst_grd={dst_grid_file} -m {map_file} '
    run_cmd(cmd)
# ------------------------------------------------------------------------------
if regrid_data:
    # remap the data
    cmd = f'ncremap {alg} -m {map_file} -i {src_file_name} -o {dst_file_name}  '
    run_cmd(cmd)

    cmd = f'ncrename -v lat,latitude -v lon,longitude {dst_file_name} '
    run_cmd(cmd)
# ------------------------------------------------------------------------------
print(f'\n\nsrc file: {src_file_name}\ndst file: {dst_file_name}\n')
# ------------------------------------------------------------------------------