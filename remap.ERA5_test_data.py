#!/usr/bin/env python
import os
import glob
import subprocess
pwd = os.getenv('PWD')

# The user might need to "unpack" the data if downloaded directly from CDS
# example: ncpdq -U ERA5.HICCUP_TEST.sfc.nc ERA5.HICCUP_TEST.sfc.upack.nc

src_file_name = 'ERA5.HICCUP_TEST.atm.upack.nc'
dst_file_name = src_file_name.replace('.upack.nc','.remap.nc')

nlat_src,nlon_src = 721,1440
nlat_dst,nlon_dst = 180,360

src_grid_file = f'{pwd}/scrip_{nlat_src}x{nlon_src}.nc'
dst_grid_file = f'{pwd}/scrip_{nlat_dst}x{nlon_dst}.nc'

map_file = f'{pwd}//map_{nlat_src}x{nlon_src}_to_{nlat_dst}x{nlon_dst}.nc'

alg = '-a tempest'

clean       = True
execute     = True

create_map  = False
regrid_data = True

#-------------------------------------------------------------------------------
if create_map:

    # remove old grid and mapping files
    if clean :
        cmd = f'rm {src_grid_file} {dst_grid_file} {map_file}'
        print('\n'+cmd+'\n')
        if execute:  subprocess.call(cmd, shell=True)
    
    # Generate source grid file:
    cmd  = f'ncremap {alg} -G ttl=\'Equi-Angular grid {nlat_src}x{nlon_src}\'#latlon={nlat_src},{nlon_src}'+\
           f'#lat_typ=uni#lat_drc=s2n#lon_typ=grn_ctr -g {src_grid_file} '
    print('\n'+cmd)
    if execute:  subprocess.call(cmd, shell=True)
    
    # Generate target grid file:
    cmd  = f'ncremap {alg} -G ttl=\'Equi-Angular grid {nlat_dst}x{nlon_dst}\'#latlon={nlat_dst},{nlon_dst}'+\
           f'#lat_typ=uni#lat_drc=s2n#lon_typ=grn_ctr -g {dst_grid_file} '
    print('\n'+cmd)
    if execute:  subprocess.call(cmd, shell=True)
    
    # Need to make sure the 'grid_imask' variable is an integer for TempestRemap
    cmd = f'ncap2 -s \'grid_imask=int(grid_imask)\' {src_grid_file} {src_grid_file} --ovr'
    print('\n'+cmd+'\n')
    if execute:  subprocess.call(cmd, shell=True)
    cmd = f'ncap2 -s \'grid_imask=int(grid_imask)\' {dst_grid_file} {dst_grid_file} --ovr'
    print('\n'+cmd+'\n')
    if execute:  subprocess.call(cmd, shell=True)
    
    # Generate mapping file:
    cmd  = f'ncremap {alg} -a fv2fv --src_grd={src_grid_file} --dst_grd={dst_grid_file} -m {map_file} '
    print('\n'+cmd+'\n')
    if execute:  subprocess.call(cmd, shell=True)
#-------------------------------------------------------------------------------
if regrid_data:
    # remap the data
    cmd = f'ncremap {alg} -m {map_file} -i {src_file_name} -o {dst_file_name}  '
    print('\n'+cmd+'\n')
    if execute: subprocess.call(cmd, shell=True)

    cmd = f'ncrename -v lat,latitude -v lon,longitude {dst_file_name} '
    print('\n'+cmd+'\n')
    if execute: subprocess.call(cmd, shell=True)
#-------------------------------------------------------------------------------
print(f'\n\nsrc file: {src_file_name}\ndst file: {dst_file_name}\n')
#-------------------------------------------------------------------------------