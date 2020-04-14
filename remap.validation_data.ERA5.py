#!/usr/bin/env python
import os
import glob
import subprocess as sp

# Generate map file for E3SM ne30np4: 
# ncremap --alg_typ=tempest --src_grd=./files_grid/exodus_ne30.g --dst_grd=./files_grid/scrip_90x180_s2n.nc --map_file=/global/homes/w/whannah/maps/map_ne30np4_to_90x180.nc --wgt_opt='--in_type cgll --in_np 4 --out_type fv --out_np 2 --out_double'

# Generate map file for E3SM ne30pg2: 
# ncremap --alg_typ=tempest --src_grd=./files_grid/exodus_ne30pg2.nc --dst_grd=./files_grid/scrip_90x180_s2n.nc --map_file=/global/homes/w/whannah/maps/map_ne30pg2_to_90x180.nc --wgt_opt='--in_type fv --in_np 2 --out_type fv --out_np 2 --out_double'

# Remap the data: ncremap -m <map file> -i <input file> -o <output file>
# obs example: FILE=data_scratch/ERA5_validation.Z.2016-08-01 ; ncremap -m ./map_files/map_721x1440_n2s_to_90x180_s2n.nc -i $FILE.nc -o $FILE.remap_90x180.nc
# hindcast example:
# export MSCRATCH=/global/cscratch1/sd/whannah/e3sm_scratch/cori-knl/
# CASE=E3SM_HINDCAST-TEST_2016-08-01_ne30_FC5AV1C-L_00 ; FILE=$MSCRATCH/$CASE/run/$CASE.cam.h1.2016-08-01-00000 ; ncremap -m $HOME/maps/map_ne30np4_to_90x180.nc -i $FILE.nc -o $FILE.remap_90x180.nc
# CASE=E3SM_HINDCAST-TEST_2016-08-01_ne30pg2_FC5AV1C-L_00 ; FILE=$MSCRATCH/$CASE/run/$CASE.cam.h1.2016-08-01-00000 ; ncremap -m $HOME/maps/map_ne30pg2_to_90x180.nc -i $FILE.nc -o $FILE.remap_90x180.nc

hiccup_root = os.getenv('HOME')+'/HICCUP/'

nlat_src,nlon_src = 721,1440
nlat_dst,nlon_dst = 90,180

date = '2016-08-01'

var_list = ['Z','T','Q','U']

alg = '-a tempest'

def main():

    unpack      = True
    clean       = True
    create_map  = True
    regrid_data = True

    src_grid_name = f'{nlat_src}x{nlon_src}_n2s'
    dst_grid_name = f'{nlat_dst}x{nlon_dst}_s2n'

    src_grid_file = f'{hiccup_root}/files_grid/scrip_{src_grid_name}.nc'
    dst_grid_file = f'{hiccup_root}/files_grid/scrip_{dst_grid_name}.nc'

    map_file = f'{hiccup_root}/map_files/map_{nlat_src}x{nlon_src}_to_{nlat_dst}x{nlon_dst}.nc'

    # --------------------------------------------------------------------------
    for var in var_list:

        src_file_name = f'{hiccup_root}/data_scratch/ERA5_validation.{var}.{date}.nc'
        dst_file_name = src_file_name.replace('.nc',f'.remap_{nlat_dst}x{nlon_dst}.nc')
        
        if unpack:
            cmd = f'ncpdq --ovr -U {src_file_name} {src_file_name}'
            run_cmd(cmd)
        # ----------------------------------------------------------------------
        if create_map:

            # remove old grid and mapping files
            if clean :
                if os.path.exists(src_grid_file) : run_cmd(f'rm {src_grid_file}')
                if os.path.exists(dst_grid_file) : run_cmd(f'rm {dst_grid_file}')
                if os.path.exists(map_file)      : run_cmd(f'rm {map_file}')
            
            # Generate source grid file:
            cmd  = f'ncremap {alg} -G '
            cmd += f'ttl=\'Equi-Angular grid {nlat_src}x{nlon_src}\''
            cmd += f'#latlon={nlat_src},{nlon_src}'
            cmd += f'#lat_typ=uni'
            cmd += f'#lat_drc=n2s'
            cmd += f'#lon_typ=grn_ctr '
            cmd += f'-g {src_grid_file} '
            run_cmd(cmd)
            
            # Generate target grid file:
            cmd  = f'ncremap {alg} -G '
            cmd += f'ttl=\'Equi-Angular grid {nlat_dst}x{nlon_dst}\''
            cmd += f'#latlon={nlat_dst},{nlon_dst}'
            cmd += f'#lat_typ=uni'
            cmd += f'#lat_drc=s2n'
            cmd += f'#lon_typ=grn_ctr '
            cmd += f'-g {dst_grid_file} '
            run_cmd(cmd)
            
            # Need to make sure the 'grid_imask' variable is an integer for TempestRemap
            # run_cmd(f'ncap2 -s \'grid_imask=int(grid_imask)\' {src_grid_file} {src_grid_file} --ovr')
            # run_cmd(f'ncap2 -s \'grid_imask=int(grid_imask)\' {dst_grid_file} {dst_grid_file} --ovr')
            
            # Generate mapping file:
            cmd  = f'ncremap {alg} -a fv2fv'
            cmd += f' --src_grd={src_grid_file}'
            cmd += f' --dst_grd={dst_grid_file}'
            cmd += f' -m {map_file} '
            run_cmd(cmd)
        # ----------------------------------------------------------------------
        # remap the data
        run_cmd(f'ncremap {alg} -m {map_file} -i {src_file_name} -o {dst_file_name}  ')

        # run_cmd(f'ncrename -v lat,latitude -v lon,longitude {dst_file_name} ')
        # --------------------------------------------------------------------------
        print(f'\n\nsrc file: {src_file_name}\ndst file: {dst_file_name}\n')

# --------------------------------------------------------------------------------------------------
class tcolor:
    """ simple class for coloring terminal text """
    ENDC, BLACK, RED     = '\033[0m','\033[30m','\033[31m'
    GREEN, YELLOW, BLUE  = '\033[32m','\033[33m','\033[34m'
    MAGENTA, CYAN, WHITE = '\033[35m','\033[36m','\033[37m'
# --------------------------------------------------------------------------------------------------
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
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__': 
   # # Parse the command line options
   # parser = OptionParser()
   # parser.add_option('-i',dest='ifile',default=None,help='input file name')
   # (opts, args) = parser.parse_args()

   main()