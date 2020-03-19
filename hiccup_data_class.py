# Class structure for defining parameters specific to the 
# input data used for creating initial condition files
# 
# NOTE: Variable name dictionaries are defined with the key as the model's 
# variable name and the value as the reanalysis data variable name

import numpy as np
import xarray as xr
import pandas as pd
import subprocess as sp
import datetime
# import cftime
import shutil 
import re
import glob
import os

# default output paths
default_output_dir  = './data/'
default_grid_dir    = './grid_files/'
default_map_dir     = './map_files/'
default_tmp_dir     = './tmp/'

# algorithm flag for ncremap
ncremap_alg         = ' --alg_typ=tempest '    

# log file for Tempest output
tempest_log_file    = 'TempestRemap.log'

# override the xarray default netcdf format of 
# NETCDF4 to avoid file permission issue
hiccup_sst_nc_format = 'NETCDF3_64BIT'
ncremap_fl_fmt = '64bit_data' # netcdf4_classic / 64bit_data / 64bit_offset

# Global verbosity default
hiccup_verbose = False

# Set numpy to ignore overflow errors
np.seterr(over='ignore')

# Numeric parameters
tk_zero = 273.15 # value for converting between celsius and Kelvin

# Set up terminal colors
class tcolor:
    ENDC, BLACK, RED     = '\033[0m','\033[30m','\033[31m'
    GREEN, YELLOW, BLUE  = '\033[32m','\033[33m','\033[34m'
    MAGENTA, CYAN, WHITE = '\033[35m','\033[36m','\033[37m'



# ------------------------------------------------------------------------------
# Common method for printing and running commands
# ------------------------------------------------------------------------------
def run_cmd(cmd,verbose=None,prepend_line=True,use_color=True,shell=False):
    """
    Method to encapsulate running system commands and checking for failures
    """
    prefix='  '
    suffix=''
    if prepend_line : prefix = '\n'+prefix
    if verbose is None : verbose = hiccup_verbose
    msg = f'{prefix}{cmd}{suffix}'
    if use_color : msg = tcolor.GREEN + msg + tcolor.ENDC
    if verbose : print(msg)
    if shell:
        sp.check_call(cmd,shell=True)
    else:
        sp.check_call(cmd.split())
    return
# ------------------------------------------------------------------------------
# Method for checking if required software is installed
# ------------------------------------------------------------------------------
def check_dependency(cmd):
    """ 
    Check for required system commands 
    """
    if shutil.which(cmd) is None : raise OSError(f'{cmd} is not in system path')
    return
# ------------------------------------------------------------------------------
# Method for returning class object
# ------------------------------------------------------------------------------
def create_hiccup_data(name,atm_file,sfc_file,dst_horz_grid,dst_vert_grid,
                       output_dir=default_output_dir,grid_dir=default_grid_dir,
                       map_dir=default_map_dir,tmp_dir=default_tmp_dir,
                       sstice_combined_file=None,sstice_name=None,
                       sst_file=None,ice_file=None,
                       lev_type='',verbose=False):
    """ 
    Create HICCUP data class object, check for required input variables and 
    create specified output directories if they do not exist
    """
    global hiccup_verbose
    hiccup_verbose = verbose
    for subclass in hiccup_data.__subclasses__():
        if subclass.is_name_for(name):
            # Create the object
            obj = subclass(name
                      ,atm_file=atm_file
                      ,sfc_file=sfc_file
                      ,sstice_name=sstice_name
                      ,sst_file=sst_file
                      ,ice_file=ice_file
                      ,sstice_combined_file=sstice_combined_file
                      ,dst_horz_grid=dst_horz_grid
                      ,dst_vert_grid=dst_vert_grid
                      ,output_dir=output_dir
                      ,grid_dir=grid_dir
                      ,map_dir=map_dir
                      ,tmp_dir=tmp_dir
                      ,lev_type=lev_type)
            
            # Check input files for for required variables
            obj.check_file_vars()

            # Create the output, grid, and map folders if they do not exist
            if not os.path.exists(output_dir) : os.makedirs(output_dir)
            if not os.path.exists(grid_dir)   : os.makedirs(grid_dir)
            if not os.path.exists(map_dir)    : os.makedirs(map_dir)

            # Return the object if everything checks out
            return obj
    raise ValueError(f'{name} is not a valid HICCUP dataset name')
# ------------------------------------------------------------------------------
# Base Class
# ------------------------------------------------------------------------------
class hiccup_data(object):
    """ 
    Base class for HICCUP data object 
    """
    def __init__(self,atm_file,sfc_file,dst_horz_grid,dst_vert_grid,
                 output_dir=default_output_dir,grid_dir=default_grid_dir,
                 map_dir=default_map_dir,tmp_dir=default_tmp_dir,
                 sstice_combined_file=None,sstice_name=None,
                 sst_file=None,ice_file=None,lev_type=''):
        self.lev_type = lev_type
        self.atm_file = atm_file
        self.sfc_file = sfc_file
        self.atm_var_name_dict = {}
        self.sfc_var_name_dict = {}
        self.lnd_var_name_dict = {}
        self.nlat = -1
        self.nlon = -1
        self.dst_horz_grid = dst_horz_grid
        self.dst_vert_grid = dst_vert_grid
        self.src_grid_name = ''
        self.src_grid_file = None
        self.dst_grid_file = None
        self.map_file = None

        self.sstice_name = sstice_name
        self.sst_file = sst_file
        self.ice_file = ice_file
        self.sstice_combined_file = sstice_combined_file
        self.sstice_nlat_src = None
        self.sstice_nlon_src = None
        self.sstice_nlat_dst = None
        self.sstice_nlon_dst = None

        # Set output paths for data, grid, and map files
        if output_dir=='' or output_dir==None : output_dir = './'
        if grid_dir=='' or grid_dir==None : grid_dir = default_grid_dir
        if map_dir=='' or map_dir==None : map_dir = default_map_dir
        if tmp_dir=='' or tmp_dir==None : tmp_dir = default_tmp_dir

        # Make sure directory strings are formatted with trailing slash
        if not output_dir.endswith('/'): output_dir += '/'
        if not grid_dir.endswith('/'): grid_dir += '/'
        if not map_dir.endswith('/'): map_dir += '/'
        if not tmp_dir.endswith('/'): tmp_dir += '/'

        self.output_dir = output_dir
        self.grid_dir = grid_dir
        self.map_dir = map_dir
        self.tmp_dir = tmp_dir

        # Check if sst/ice dataset is supported
        # if self.sstice_name not in [None,'NOAA']: 
        #     err_msg = f'sstice_name={self.sstice_name} is not currently supported'
        #     if self.sstice_name=='ERA5': 
        #         err_msg += ' due to an issue handling missing values during remapping'
        #     raise ValueError(err_msg)
            
        # set input variable names for SST and sea ice 
        if self.sstice_name=='NOAA': self.sst_name,self.ice_name = 'sst','icec'
        if self.sstice_name=='ERA5': self.sst_name,self.ice_name = 'sst','siconc'

        # Load input files into xarray datasets
        self.ds_atm = xr.open_dataset(self.atm_file)
        self.ds_sfc = xr.open_dataset(self.sfc_file)
    # --------------------------------------------------------------------------
    def __str__(self):
        str_out = ''
        fmt_key_len = 18
        for key in self.__dict__.keys(): 
            attribute = getattr(self,key)
            if isinstance(attribute,dict):
                str_out += f'  {key:{fmt_key_len}}:\n'
                for k in attribute.keys(): 
                    str_out += ' '*(fmt_key_len+4)+f'{k:8}  {attribute[k]}\n'
            elif isinstance(attribute, xr.Dataset) : 
                str_out += f'  {key:{fmt_key_len}}:'
                ds_str = attribute.__str__().replace('\n','\n'+' '*(fmt_key_len+4))
                str_out += f'  {ds_str}\n'
            else:
                if attribute!='' :
                    str_out += f'  {key:{fmt_key_len}}:  {attribute}\n'

        return str_out
    # --------------------------------------------------------------------------
    def check_file_vars(self):
        """ 
        Check that required variables are in the input files 
        """

        # Create list of variables in the files
        atm_file_vars = []
        sfc_file_vars = []
        for key in self.ds_atm.variables.keys(): atm_file_vars.append(key)
        for key in self.ds_sfc.variables.keys(): sfc_file_vars.append(key)

        # Check that all required data exists in the atm file
        for key in self.atm_var_name_dict : 
            if self.atm_var_name_dict[key] not in atm_file_vars: 
                raise ValueError(f'{self.atm_var_name_dict[key]} is not in ATM dataset: ({self.atm_file})')

        # Check that all required data exists in the sfc file
        for key in self.sfc_var_name_dict : 
            if self.sfc_var_name_dict[key] not in sfc_file_vars: 
                raise ValueError(f'{self.sfc_var_name_dict[key]} is not in SFC dataset: ({self.sfc_file})')

        return
    # --------------------------------------------------------------------------
    def unpack_data_files(self,verbose=None):
        """
        Make sure data files are unpacked
        """
        if verbose is None : verbose = hiccup_verbose
        if verbose : print('\nUnpacking data files...')

        check_dependency('ncpdq')

        for f in [ self.atm_file, self.sfc_file, 
                   self.sst_file, self.ice_file,
                   self.sstice_combined_file ]:
            if f is not None :
                run_cmd(f'ncpdq -U --ovr {f} {f}',verbose,prepend_line=False)
        return
    # --------------------------------------------------------------------------
    def create_dst_grid_file(self,verbose=None):
        """ 
        Generate destination model grid file 
        """
        if verbose is None : verbose = hiccup_verbose
        if verbose : print('\nGenerating dst grid file...')

        if 'ne' in self.dst_horz_grid and 'np' in self.dst_horz_grid : 
            
            # Spectral element grid with physics on GLL nodes
            ne = re.search('ne(.*)np', self.dst_horz_grid).group(1)
            self.dst_grid_file = self.grid_dir+f'exodus_ne{ne}.g'
            
            check_dependency('GenerateCSMesh')
            cmd = f'GenerateCSMesh --res {ne} --file {self.dst_grid_file}'
            cmd += f' >> {tempest_log_file}'
            run_cmd(cmd,verbose,shell=True)

        elif 'ne' in self.dst_horz_grid and 'pg' in self.dst_horz_grid : 
            
            # Spectral element grid with FV physics grid (ex. ne30pg2)
            ne  = re.search('ne(.*)pg', self.dst_horz_grid).group(1)
            npg = re.search('pg(.*)', self.dst_horz_grid).group(1)
            exodus_file = self.grid_dir+f'exodus_ne{ne}.g'

            # First create exodus file
            check_dependency('GenerateCSMesh')
            cmd = f'GenerateCSMesh --res {ne} --file {exodus_file}'
            cmd += f' >> {tempest_log_file}'
            run_cmd(cmd,verbose,shell=True)
            
            # Next switch to volumetric mesh that matches the physgrid
            self.dst_grid_file = self.grid_dir+f'exodus_{self.dst_horz_grid}.nc'
            check_dependency('GenerateVolumetricMesh')
            cmd = 'GenerateVolumetricMesh'
            cmd += f' --in {exodus_file} '
            cmd += f' --out {self.dst_grid_file} '
            cmd += f' --np {npg} --uniform'
            cmd += f' >> {tempest_log_file}'
            run_cmd(cmd,verbose,shell=True)

            # # Create scrip file while we're at it (can be slow)
            check_dependency('ConvertExodusToSCRIP')
            scrip_file = self.grid_dir+f'scrip_{self.dst_horz_grid}.nc'
            cmd = 'ConvertExodusToSCRIP'
            cmd += f' --in {self.dst_grid_file} '
            cmd += f' --out {scrip_file} '
            cmd += f' >> {tempest_log_file}'
            run_cmd(cmd,verbose,shell=True)

        else:
            raise ValueError(f'grid_name={self.dst_horz_grid} is not currently supported')

        return 
    # --------------------------------------------------------------------------
    def create_map_file(self,verbose=None):
        """ 
        Generate mapping file after grid files have been created 
        """
        if verbose is None : verbose = hiccup_verbose
        if verbose : print('\nGenerating mapping file...')

        check_dependency('ncremap')

        # Check that grid file fields are not empty
        if self.src_grid_file == None : 
            raise ValueError('src_grid_file is not defined for hiccup_data object')
        if self.dst_grid_file == None : 
            raise ValueError('dst_grid_file is not defined for hiccup_data object')

        # Set the map options
        self.map_opts = '--in_type fv --in_np 1 --mono --out_double '

        # speciic special options depending on target atmos grid
        if 'ne' in self.dst_horz_grid and 'np' in self.dst_horz_grid : 
            self.map_opts = self.map_opts+' --out_type cgll --out_np 4 ' # options for SE grid
            ne = re.search('ne(.*)np', self.dst_horz_grid).group(1)
        elif 'ne' in self.dst_horz_grid and 'pg' in self.dst_horz_grid :
            self.map_opts = self.map_opts+' --out_type fv --out_np 1 --volumetric '
            ne = re.search('ne(.*)pg', self.dst_horz_grid).group(1)
        else:
            raise ValueError(f'dst_horz_grid={self.dst_horz_grid} does not seem to be valid')
        

        cmd = f'ncremap {ncremap_alg} '
        cmd += f' --src_grd={self.src_grid_file}'
        cmd += f' --dst_grd={self.dst_grid_file}'
        cmd += f' --map_file={self.map_file}'
        cmd += f' --wgt_opt=\'{self.map_opts}\' '
        # Add special flag for "very fine" grids
        if int(ne) > 100 : cmd += ' --lrg2sml '
        run_cmd(cmd,verbose,shell=True)

        return
    # --------------------------------------------------------------------------
    def rename_vars(self,file_name,verbose=None):
        """ 
        Rename variables in file according to variable name dictionaries 
        """
        if verbose is None : verbose = hiccup_verbose
        if verbose : print('\nRenaming variables to match model variable names...')

        check_dependency('ncrename')

        def rename_proc(key,var_name_dict):
            # Skip over entries that are blank
            if key == '' : return
            if var_name_dict[key] == '' : return
            # set up the ncrename command
            cmd = f'ncrename --hst --variable {var_name_dict[key]},{key} {file_name}'
            # coords are often already renamed by the remapping step, 
            # so make them optional by adding a preceeding dot
            if key in ['lat','lon']: cmd = cmd.replace(f'{var_name_dict[key]}',f'.{var_name_dict[key]}')
            # print the command and execute
            run_cmd(cmd,verbose,prepend_line=False,shell=True)

        if verbose : print('\n Renaming ATM vars... \n')
        for key in self.atm_var_name_dict : rename_proc(key,self.atm_var_name_dict)

        if verbose : print('\n Renaming SFC vars... \n')
        for key in self.sfc_var_name_dict : rename_proc(key,self.sfc_var_name_dict)

        self.rename_vars_special(file_name,verbose)

        return
    # --------------------------------------------------------------------------
    def add_reference_pressure(self,file_name,verbose=None):
        """ 
        Add P0 variable 
        """
        if verbose is None : verbose = hiccup_verbose
        if verbose : print('\nAdding reference pressure (P0)...')

        check_dependency('ncap2')
        check_dependency('ncatted')

        # Add the variable
        run_cmd(f"ncap2 --hst -A -s 'P0=100000.' {file_name} {file_name}",
                verbose,prepend_line=False,shell=True)

        # add long_name attribute
        run_cmd(f"ncatted --hst -A -a long_name,P0,a,c,'reference pressure' {file_name}",
                verbose,prepend_line=False,shell=True)
        
        # add units attribute
        run_cmd(f"ncatted --hst -A -a units,P0,a,c,'Pa' {file_name}",
                verbose,prepend_line=False,shell=True)

        return
    # --------------------------------------------------------------------------
    def remap_horizontal(self,output_file_name,verbose=None):
        """  
        Horizontally remap data and combine into single file 
        """
        if verbose is None : verbose = hiccup_verbose
        if verbose : print('\nHorizontally remapping the data to temporary files...')

        if self.map_file is None: raise ValueError('map_file cannot be None!')
        if self.atm_file is None: raise ValueError('atm_file cannot be None!')
        if self.sfc_file is None: raise ValueError('sfc_file cannot be None!')

        # Define temporary files that will be deleted at the end
        atm_tmp_file_name = f'{self.tmp_dir}tmp_atm_data.nc'
        sfc_tmp_file_name = f'{self.tmp_dir}tmp_sfc_data.nc'

        # Remove temporary files if they exist
        if atm_tmp_file_name in glob.glob(atm_tmp_file_name):
            run_cmd(f'rm {atm_tmp_file_name} ',verbose)
        if sfc_tmp_file_name in glob.glob(sfc_tmp_file_name):
            run_cmd(f'rm {sfc_tmp_file_name} ',verbose)

        check_dependency('ncremap')
        check_dependency('ncks')

        # Horzontally remap atmosphere data
        var_list = ','.join(self.atm_var_name_dict.values())
        cmd =  f'ncremap {ncremap_alg} '
        cmd += f' --map_file={self.map_file} '
        cmd += f' --in_file={self.atm_file} '
        cmd += f' --out_file={atm_tmp_file_name} '
        cmd += f' --var_lst={var_list} '
        cmd += f' --fl_fmt={ncremap_fl_fmt} '
        run_cmd(cmd,verbose)

        # Horzontally remap surface data
        var_list = ','.join(self.sfc_var_name_dict.values())
        cmd =  f'ncremap {ncremap_alg} '
        cmd += f' --map_file={self.map_file} '
        cmd += f' --in_file={self.sfc_file} '
        cmd += f' --out_file={sfc_tmp_file_name} '
        cmd += f' --var_lst={var_list} '
        cmd += f' --fl_fmt={ncremap_fl_fmt} '
        run_cmd(cmd,verbose)

        # Remove output file if it already exists
        if output_file_name in glob.glob(output_file_name) : run_cmd(f'rm {output_file_name} ',verbose)

        if verbose : print('\nCombining temporary remapped files...')

        # Add atmosphere temporary file data into the final output file
        run_cmd(f'ncks -A {atm_tmp_file_name} {output_file_name} ',
                verbose,prepend_line=False)

        # Add surface temporary file data into the final output file
        run_cmd(f'ncks -A {sfc_tmp_file_name} {output_file_name} ',
                verbose,prepend_line=False)

        # delete the temporary files
        run_cmd(f'rm {sfc_tmp_file_name} {atm_tmp_file_name} ',
                verbose,prepend_line=False)

        return
    # --------------------------------------------------------------------------
    def remap_vertical(self,input_file_name,output_file_name,
                       vert_file_name,vert_remap_var_list=None,
                       verbose=None):
        """  
        Vertically remap data and combine into single file 
        """
        if verbose is None : verbose = hiccup_verbose
        if verbose : print('\nVertically remapping the data...')

        check_dependency('ncremap')

        # Specify temporary file for vertically interpolated output
        # This allows for input and output files to be the same
        vert_tmp_file_name = output_file_name.replace('.nc',f'.vert_remap_tmp.nc')

        # Build variable list from the input file if not supplied
        if vert_remap_var_list is None :
            vert_remap_var_list = []
            ds = xr.open_dataset(input_file_name)
            for key in ds.variables.keys(): 
                vert_remap_var_list.append(key)
                # only remap variables with lev coord in order to
                # ignore other variables (i.e. TS, PS) - not necessary?
                # if self.lev_name in ds[key].dims and key!=self.lev_name :
                #     vert_remap_var_list.append(key)
        vert_remap_var_list = ','.join(vert_remap_var_list)

        # Perform the vertical remapping
        cmd  = 'ncremap'
        cmd += f' --vrt_fl={vert_file_name}'
        cmd += f' --var_lst={vert_remap_var_list}'
        cmd += f' --in_fl={input_file_name}'
        cmd += f' --out_fl={vert_tmp_file_name}'
        cmd += f' --fl_fmt={ncremap_fl_fmt} '
        run_cmd(cmd,verbose,shell=True)

        # Overwrite the output file with the vertically interpolated data
        run_cmd(f'mv {vert_tmp_file_name} {output_file_name} ')

        # The command above was causing a strange issue where the file was
        # not completely written when the next step happened. Adding the lines
        # below to fixes the problem, but I do not understand why...
        ds = xr.open_dataset(output_file_name)
        ds.close()

        return

    # --------------------------------------------------------------------------
    def add_time_date_variables(self,ds,verbose=None):
        """
        Check final output file and add necessary time and date information
        """

        if verbose is None : verbose = hiccup_verbose
        if verbose : print('\nEditing time and date variables...')

        # time coordinate information
        time_shape = ( len(ds['time']) )
        time_coord = {'time':ds['time'].values}
        time_dim   = ['time']

        # Use the pandas library to parse the time/date information
        datetimeindex = pd.DatetimeIndex( ds['time'].values )
        year  = datetimeindex.year.values
        month = datetimeindex.month.values
        day   = datetimeindex.day.values
        hour  = datetimeindex.hour.values
        sec   = datetimeindex.second.values

        # Change time attributes
        ds['time'].attrs['calendar'] = 'noleap'       # this might cause xarray to throw an error
        # ds['time'].attrs['calendar'] = 'gregorian'
        ds['time'].attrs['bounds'] = 'time_bnds'

        # add time_bnds variable
        time_bnds = np.full( (len(ds['time']),2), 0. )
        for t in range(len(ds['time'])) : time_bnds[t,:] = ds['time'].values
        ds['time_bnds'] = xr.DataArray( time_bnds, coords=time_coord, dims=['time','nbnd'] )
        ds['time_bnds'].attrs['long_name'] = 'time interval endpoints'

        # string representation of date
        # date_str = str(year).zfill(4)+str(month).zfill(2)+str(day).zfill(2)

        # Integer representation of date
        date_int = year*1e4 + month*1e2 + day

        # miscellaneous time/date variables
        if 'date' not in ds.variables :
            ds['date'] = xr.DataArray( np.array(date_int,dtype=np.int), 
                                        coords=time_coord, dims=time_dim )
            ds['date'].attrs['long_name'] = 'current date (YYYYMMDD)'
        if 'datesec' not in ds.variables :
            ds['datesec'] = xr.DataArray( np.array(sec,dtype=np.int), 
                                        coords=time_coord, dims=time_dim )
            ds['datesec'].attrs['long_name'] = 'current seconds of current date'
        if 'ndcur' not in ds.variables :
            ds['ndcur'] = xr.DataArray( np.full(time_shape,0.,dtype=np.int), 
                                        coords=time_coord, dims=time_dim )
            ds['ndcur'].attrs['long_name'] = 'current day (from base day)'
        if 'nscur' not in ds.variables :
            ds['nscur'] = xr.DataArray( np.array(sec,dtype=np.int), 
                                        coords=time_coord, dims=time_dim )
            ds['nscur'].attrs['long_name'] = 'current seconds of current day'
        if 'nsteph' not in ds.variables :
            ds['nsteph'] = xr.DataArray( np.full(time_shape,0.,dtype=np.int), 
                                        coords=time_coord, dims=time_dim )
            ds['nsteph'].attrs['long_name'] = 'current timestep'

        # Base day time - not sure what this is for
        if 'nbdate' not in ds.variables :
            ds['nbdate'] = xr.DataArray( np.array(date_int,dtype=np.int), 
                                        coords=time_coord, dims=time_dim )
            ds['nbdate'].attrs['long_name'] = 'base date (YYYYMMDD)'
        if 'ndbase' not in ds.variables :
            ds['ndbase'] = xr.DataArray( np.array(day,dtype=np.int), 
                                        coords=time_coord, dims=time_dim )
            ds['ndbase'].attrs['long_name'] = 'base day'
        if 'nsbase' not in ds.variables :
            ds['nsbase'] = xr.DataArray( np.array(sec,dtype=np.int), 
                                        coords=time_coord, dims=time_dim )
            ds['nsbase'].attrs['long_name'] = 'seconds of base day'
        if 'nbsec' not in ds.variables :
            ds['nbsec'] = xr.DataArray( np.array(sec,dtype=np.int), 
                                        coords=time_coord, dims=time_dim )
            ds['nbsec'].attrs['long_name'] = 'seconds of base date'

        # (Not sure how to handle string/char variables, but are they needed?)

        # if 'date_written' not in ds.variables :
        #     ds['date_written'] = xr.DataArray( np.full(time_shape,0.), coords=time_coord, dims=time_dim )
        #     ds['date_written'].attrs['long_name'] = ''

        # if 'time_written' not in ds.variables :
        #     ds['time_written'] = xr.DataArray( np.full(time_shape,0.), coords=time_coord, dims=time_dim )
        #     ds['time_written'].attrs['long_name'] = ''

        return

    # --------------------------------------------------------------------------
    def add_extra_data_variables(self,ds) :
        """
        Check final output file and add other data variables needed by the model
        """
        shape_2D = ( len(ds['time']), len(ds['ncol']) )
        shape_3D = ( len(ds['time']), len(ds['lev']), len(ds['ncol']) )
        coord_2D = {'time':ds['time'], 'ncol':ds['ncol'] }
        coord_3D = {'time':ds['time'], 'lev':ds['lev'], 'ncol':ds['ncol'] }
        dims_2D  = ['time','ncol']
        dims_3D  = ['time','lev','ncol']

        if 'NUMICE' not in ds.variables :
            ds['NUMICE'] = xr.DataArray( np.full(shape_3D,0.), coords=coord_3D, dims=dims_3D )
            ds['NUMICE'].attrs['units'] = 'kg/kg'
            ds['NUMICE'].attrs['long_name'] = 'Grid box averaged cloud ice number'

        if 'NUMLIQ' not in ds.variables :
            ds['NUMLIQ'] = xr.DataArray( np.full(shape_3D,0.), coords=coord_3D, dims=dims_3D )
            ds['NUMLIQ'].attrs['units'] = 'kg/kg'
            ds['NUMLIQ'].attrs['long_name'] = 'Grid box averaged cloud liquid number'
        
        # Other variables that might be needed:
        # H2SO4
        # H2O2
        # DMS
        # SOAG
        # SO2

        return
    # --------------------------------------------------------------------------
    def clean_global_attributes(self,file_name,verbose=None):
        """ 
        Remove messy global attributes of the file 
        """
        if verbose is None : verbose = hiccup_verbose
        if verbose: print('\nCleaning up excessive global attributes...')
        
        global_att_list = ['history_of_appended_files', 'nco_openmp_thread_number', 
                           'input_file', 'map_file', 'remap_version', 'remap_hostname', 
                           'remap_command', 'remap_script', 'NCO' ]
        
        check_dependency('ncatted')

        # Remove the attributes listed in global_att_list
        for att in global_att_list:
            run_cmd(f'ncatted -O -a {att},global,d,, {file_name} {file_name}',
                    verbose,prepend_line=False)

        # Also reset the history attribute
        run_cmd(f'ncatted -h -O -a history,global,o,c, {file_name} {file_name}',
                verbose,prepend_line=False)
    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------
    # SST and sea ice methods
    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------
    def sstice_create_src_grid_file(self,diagnose_grid=True,nlat=None,nlon=None,
                                    force_overwrite=False,verbose=None):
        """
        Create a source grid file to use for remapping the SST and sea ice data.
        The SST and ice data are assumed to exist on the same grid.
        """
        if verbose is None : verbose = hiccup_verbose

        if diagnose_grid:
            # Load the SST file as xarray datasets to read the grid dimensions
            sst_file = self.sstice_combined_file
            if sst_file is None : sst_file = self.sst_file
            if sst_file is None : raise ValueError('No valid sst data file found!')

            ds_sst = xr.open_dataset(sst_file)

            # Determine input grid 
            if self.sstice_name=='NOAA': lat_name, lon_name = 'lat', 'lon'
            if self.sstice_name=='ERA5': lat_name, lon_name = 'latitude', 'longitude'
            self.sstice_nlat_src = len( ds_sst[lat_name].values )
            self.sstice_nlon_src = len( ds_sst[lon_name].values )

            # Close the dataset
            ds_sst.close()
        else:
            # Use supplied grid dimensions
            if nlat is None : raise ValueError('nlat can not be None if diagnose_grid is False!')
            if nlon is None : raise ValueError('nlon can not be None if diagnose_grid is False!')
            self.sstice_nlat_src = nlat
            self.sstice_nlon_src = nlon

        src_grid = f'{self.sstice_nlat_src}x{self.sstice_nlon_src}'

        # Define the source grid file to be created
        self.sstice_src_grid_file = f'{default_grid_dir}scrip_{src_grid}.nc'

        # Create the source grid file
        file_does_not_exist = self.sstice_src_grid_file not in glob.glob(self.sstice_src_grid_file) 
        if file_does_not_exist or force_overwrite :
            if verbose : print(f'\nCreating source grid file for SST and sea ice data...')
            cmd  = f'ncremap {ncremap_alg} --tmp_dir={self.tmp_dir}'
            cmd += f' -G ttl=\'Equi-Angular grid {src_grid}\'' 
            cmd += f'#latlon={self.sstice_nlat_src},{self.sstice_nlon_src}'
            cmd +=  '#lat_typ=uni'
            cmd +=  '#lon_typ=grn_ctr'
            if self.sstice_name=='NOAA': cmd += '#lat_drc=s2n'
            if self.sstice_name=='ERA5': cmd += '#lat_drc=n2s'
            cmd += f' -g {self.sstice_src_grid_file} '
            run_cmd(cmd,verbose,shell=True,prepend_line=False)

        return
    # --------------------------------------------------------------------------
    def get_sst_file(self):
        """
        Return file name that contains input SST data
        """
        sst_file = self.sstice_combined_file
        if sst_file is None : sst_file = self.sst_file
        if sst_file is None : raise ValueError('No valid sst data file found!')
        return sst_file
    # --------------------------------------------------------------------------
    def open_combined_sstice_dataset(self):
        """
        Return xarray dataset of combined SST and sea ice dataset
        """
        if self.sstice_combined_file is None :
            # if sstice_combined_file is not defined then we need to 
            # combine the individual SST and sea ice data files
            if self.sst_file is None and self.ice_file is None: 
                raise ValueError('sst_file and ice_file must be set if sstice_combined_file is not set!')
            # Combine individual SST and sea ice data files
            ds_sst = xr.open_dataset(self.sst_file).load()
            ds_ice = xr.open_dataset(self.ice_file).load()
            ds_out = xr.merge([ds_sst,ds_ice])
            ds_sst.close(); ds_ice.close()
        else:
            # Open combined sst/ice data file
            ds_out = xr.open_dataset(self.sstice_combined_file)

        return ds_out
    # --------------------------------------------------------------------------
    def sstice_create_dst_grid_file(self,output_grid_spacing=1,force_overwrite=False,
                                    verbose=None):
        """
        Create a target grid file to use for remapping the SST and sea ice data.
        The SST and ice data are assumed to exist on the same grid.
        """
        if verbose is None : verbose = hiccup_verbose

        # Define output grid dimensions
        self.sstice_nlat_dst = int( 180/output_grid_spacing )
        self.sstice_nlon_dst = int( 360/output_grid_spacing )

        dst_grid = f'{self.sstice_nlat_dst}x{self.sstice_nlon_dst}'

        # Define the destination grid file to be created
        # (add 's2n' in case input data is same grid with opposite orientation)
        self.sstice_dst_grid_file = f'{default_grid_dir}scrip_{dst_grid}_s2n.nc'

        # Create the destination grid file
        file_does_not_exist = self.sstice_dst_grid_file not in glob.glob(self.sstice_dst_grid_file)
        if file_does_not_exist or force_overwrite :
            if verbose : print(f'\nCreating target grid file for SST and sea ice data...')
            cmd  = f'ncremap {ncremap_alg} --tmp_dir={self.tmp_dir}'
            cmd += f' -G ttl=\'Equi-Angular grid {dst_grid}\'' 
            cmd += f'#latlon={self.sstice_nlat_dst},{self.sstice_nlon_dst}'
            cmd +=  '#lat_typ=uni'
            cmd +=  '#lon_typ=grn_ctr'
            cmd +=  '#lat_drc=s2n'
            cmd += f' -g {self.sstice_dst_grid_file} '
            run_cmd(cmd,verbose,shell=True)

        return
    # --------------------------------------------------------------------------
    def sstice_create_map_file(self,force_overwrite=False,verbose=None):
        """
        Create a mapping file to be used for SST and sea ice data
        """
        if verbose is None : verbose = hiccup_verbose

        src_grid = f'{self.sstice_nlat_src}x{self.sstice_nlon_src}'
        dst_grid = f'{self.sstice_nlat_dst}x{self.sstice_nlon_dst}'

        self.sstice_map_file = f'{default_map_dir}map_{src_grid}_to_{dst_grid}_s2n.nc'

        # Generate mapping file
        file_does_not_exist = self.sstice_map_file not in glob.glob(self.sstice_map_file)
        if file_does_not_exist or force_overwrite :
            if verbose : print(f'\nCreating mapping file for SST and sea ice data...')
            cmd  = f'ncremap {ncremap_alg} '
            cmd +=  ' -a fv2fv '
            cmd += f' --src_grd={self.sstice_src_grid_file}'
            cmd += f' --dst_grd={self.sstice_dst_grid_file}'
            cmd += f' --map_file={self.sstice_map_file}'
            run_cmd(cmd,verbose,shell=True,prepend_line=False)

        return
    # --------------------------------------------------------------------------
    def sstice_slice_and_remap(self,output_file_name,time_slice_method='initial',
                               atm_file=None,verbose=None):
        """
        Horizontally remap the SST and sea ice data after time slicing and 
        combining SST and sea ice files (if they are separate)
        Extract a temporal subset of the SST and sea ice data
        and combine into new temporary file for later regridding
        method options (only initial is currently supported):
        - initial           just use first time
        - match_atmos       find time that corresponds to time in atmos data
        - date_range_avg    produce an average over a specified date range
        - use_all           one output file per time slice
        """
        if verbose is None : verbose = hiccup_verbose
        if verbose : print(f'\nTime slicing {self.sstice_name} SST and sea ice data...')

        check_dependency('ncatted')
        check_dependency('ncremap')

        # Define temporary file to hold the time sliced data for regridding
        sstice_tmp_file_name = f'{self.tmp_dir}tmp_sstice_timeslice_data.nc'

        # Check that the time_slice_method is supported
        if time_slice_method not in ['initial','match_atmos'] :
            raise ValueError(f'time_slice_method: {time_slice_method} is not a supported method')

        if time_slice_method == 'initial' :
            time_slice = 0
        if time_slice_method == 'match_atmos' :
            if atm_file is None : 
                raise ValueError('atm_file can not be None for time_slice_method= \'match_atmos\'')
            sst_file = self.get_sst_file()
            ds_sst = xr.open_dataset(sst_file)
            ds_atm = xr.open_dataset(atm_file)
            time_chk = ds_atm['time'].values[0] == ds_sst['time'].values
            ds_sst.close(); ds_atm.close()
            time_slice = np.flatnonzero(time_chk)
            if time_slice.size == 0 :
                raise ValueError('No matching time slice found between SST and atmos data!')
            # If there are multiple matching times, just use the first instance
            if time_slice.size > 1 : time_slice = time_slice[0]

        # Load time sliced data
        ds_out = self.open_combined_sstice_dataset().isel(time=time_slice).load()
        # Drop any extra variables
        for var in ds_out.variables :
            if var not in [self.sst_name,self.ice_name] and var not in ds_out.coords : 
                ds_out = ds_out.drop_vars(var)
        # write out to temporary file
        ds_out.to_netcdf(sstice_tmp_file_name,format=hiccup_sst_nc_format)
        ds_out.close()

        # Replace nan values with missing_value to avoid remapping issues
        if self.sstice_name=='NOAA':
            missing_value = 1e36
            # ds = xr.open_dataset(sstice_tmp_file_name).load()
            # ds[self.sst_name] = ds[self.sst_name].where( ds[self.sst_name].values != np.nan, missing_value )
            # ds.to_netcdf(sstice_tmp_file_name,format=hiccup_sst_nc_format)
            # ds.close()
            cmd = f'ncatted -h -O -a   xxxx,o,f,{missing_value} {sstice_tmp_file_name}'
            run_cmd(cmd.replace('xxxx',   f'_FillValue,{self.sst_name}'),verbose,shell=True,prepend_line=False)
            run_cmd(cmd.replace('xxxx',   f'_FillValue,{self.ice_name}'),verbose,shell=True,prepend_line=False)
            run_cmd(cmd.replace('xxxx',f'missing_value,{self.sst_name}'),verbose,shell=True,prepend_line=False)
            run_cmd(cmd.replace('xxxx',f'missing_value,{self.ice_name}'),verbose,shell=True,prepend_line=False)
        # if self.sstice_name=='ERA5':
        #     missing_value = 1e36
        #     # ds = xr.open_dataset(sstice_tmp_file_name).load()
        #     # ds.fillna(value=missing_value)
        #     # ds.to_netcdf(sstice_tmp_file_name,format=hiccup_sst_nc_format)
        #     # ds.close()
        #     cmd = f'ncatted -h -O -a   xxxx,o,f,{missing_value} {sstice_tmp_file_name}'
        #     run_cmd(cmd.replace('xxxx',   f'_FillValue,{self.sst_name}'),verbose,shell=True,prepend_line=False)
        #     run_cmd(cmd.replace('xxxx',   f'_FillValue,{self.ice_name}'),verbose,shell=True,prepend_line=False)
        #     run_cmd(cmd.replace('xxxx',f'missing_value,{self.sst_name}'),verbose,shell=True,prepend_line=False)
        #     run_cmd(cmd.replace('xxxx',f'missing_value,{self.ice_name}'),verbose,shell=True,prepend_line=False)

        # run_cmd(f'ncdump -h {sstice_tmp_file_name}')
        # exit()

        if verbose : print(f'\nRemapping {self.sstice_name} SST and sea ice data...')

        # remap the SST data onto the target grid for the model
        cmd =  f'ncremap {ncremap_alg} '
        cmd += f' --vars={self.sst_name},{self.ice_name} '
        cmd += f' --map_file={self.sstice_map_file} '
        cmd += f' --in_file={sstice_tmp_file_name} '
        cmd += f' --out_file={output_file_name} '
        run_cmd(cmd,verbose,shell=True)

        # delete the temporary file
        run_cmd(f'rm {sstice_tmp_file_name}',verbose)

        return
    # --------------------------------------------------------------------------
    def sstice_rename_vars(self, output_file_name, new_sst_name='SST_cpl',
                           new_ice_name='ice_cov', verbose=None):
        """
        Rename sst and sea icea variables and remove unnecessary variables
        """
        if verbose is None : verbose = hiccup_verbose
        if verbose : print('\nRenaming SST and sea ice variables...')

        check_dependency('ncrename')
        check_dependency('ncks')
        check_dependency('ncatted')

        # rename variables
        cmd = f'ncrename --hst --variable {self.sst_name},{new_sst_name} {output_file_name}'
        run_cmd(cmd,verbose,shell=True)

        cmd = f'ncrename --hst --variable {self.ice_name},{new_ice_name} {output_file_name}'
        run_cmd(cmd,verbose,shell=True)

        # Make sure dimensions names are correct
        run_cmd(f'ncrename --hst --dimension .longitude,lon {output_file_name}',verbose,shell=True)
        run_cmd(f'ncrename --hst --dimension .latitude,lat  {output_file_name}',verbose,shell=True)

        # Drop unnecessary variables
        run_cmd(f'ncks -C -O -x -v area,gw,lat_bnds,lon_bnds {output_file_name} {output_file_name}',
                verbose,prepend_line=False,shell=True)

        # Drop unnecessary attributes
        global_att_list = ['history_of_appended_files']
        for att in global_att_list:
            run_cmd(f'ncatted -O -a {att},global,d,, {output_file_name} {output_file_name}',
                    verbose,prepend_line=False)

        # Reset the history attribute
        run_cmd(f'ncatted -h -O -a history,global,o,c, {output_file_name} {output_file_name}',
                verbose,prepend_line=False)

        return
    # --------------------------------------------------------------------------
    def sstice_adjustments(self, output_file_name, verbose=None):
        """
        Perform miscellaneous adjustments to the final SST/ice file
        - make sure SST units are Celsius
        - limit sea ice fraction  
        - interpolate to fill in missing SST data, extrapolate where necessary (i.e. Antarctica)
        - make sure time is a coordinate
        - add date and datesec variables
        """
        if verbose is None : verbose = hiccup_verbose
        if verbose : print('\nAdjusting SST and sea ice data values...')

        # Open remapped and combined data file
        ds = xr.open_dataset(output_file_name).load()

        # Convert units to Celsius
        if 'units' not in ds['SST_cpl'].attrs:
            ds['SST_cpl'].attrs['units'] = 'degrees_C'
        if ds['SST_cpl'].attrs['units'] in ['K','degrees_K','Kelvin'] :
            ds['SST_cpl'] = ds['SST_cpl']-tk_zero
            ds['SST_cpl'].attrs['units'] = 'degrees_C'

        # Set invalid values to np.nan before using interpolate_na()
        ds['SST_cpl'] = xr.where( np.fabs(ds['SST_cpl'].values) < 999, ds['SST_cpl'], np.nan)

        # fill in missing SST values over continents by linearly interpolating
        ds['SST_cpl'] = ds['SST_cpl'].interpolate_na(dim='lon',period=360)

        # extrapolate in latitude direction to deal with poles (requires scipy)
        ds['SST_cpl'] = ds['SST_cpl'].interpolate_na(dim='lat',method='nearest'
                                                    ,fill_value='extrapolate')

        # Fill missing ice values with zero
        ds['ice_cov'] = ds['ice_cov'].where( ds['ice_cov']>=0, 0)
        ds['ice_cov'] = ds['ice_cov'].where( ds['ice_cov']<=2, 0)

        # Also limit ice values slightly above 1 that might result from remap
        ds['ice_cov'] = ds['ice_cov'].where( ds['ice_cov']<=1, 1)

        # Make sure time is a coordinate and dimension of the dataset
        if 'time' not in ds.dims : ds = ds.expand_dims('time',axis=0)
        if 'time' not in ds.coords : ds = ds.assign_coords(coords={'time':ds['time']})

        # # For single time value, append extra dummy time for temporal interpolation
        # # STILL NOT SURE IF WE NEED THIS OR NOT
        # if ds['time'].size == 1:
        #     dt = np.timedelta64(10,'D')
        #     ds_dummy = ds.copy(deep=True)
        #     ds_dummy['time'] = ds_dummy['time']-dt
        #     ds = xr.concat([ds_dummy,ds],dim='time')
        #     ds_dummy['time'] = ds_dummy['time']+dt*2
        #     ds = xr.concat([ds,ds_dummy],dim='time')
        #     ds_dummy.close()

        # Create date and datesec variables
        time_index = pd.DatetimeIndex( ds['time'].values )
        date = np.array( time_index.year*1e4+time_index.month*1e2+time_index.day, dtype=np.int )
        datesec = np.array( time_index.second, dtype=np.int )

        # # Create alternate time index - doesn't work
        # cftime_index = pd.to_datetime(time_index).astype('cftime.DatetimeNoLeap')
        # cftime_index = pd.to_datetime(time_index)

        # # Create cftime index - doesn't work
        # cftime_index = cftime.datetime(year=time_index.year[0]
        #                               ,month=time_index.month[0]
        #                               ,day=time_index.day[0]
        #                               ,hour=time_index.hour[0]
        #                               ,minute=time_index.minute[0]
        #                               ,second=time_index.second [0] )


        # Add date and datesec variables to dataset
        ds['date'] = xr.DataArray( date, coords={'time':ds['time'].values}, dims=['time'] )
        ds['date'].attrs['long_name'] = 'current date (YYYYMMDD)'
        ds['datesec'] = xr.DataArray( datesec, coords={'time':ds['time'].values}, dims=['time'] )
        ds['datesec'].attrs['long_name'] = 'current seconds of current date'

        # ----------------------------------------------------------------------
        # Experimental - Open climatological data and use it to replace time coordinate
        # ----------------------------------------------------------------------
        # ds2 = xr.open_dataset('/global/cfs/cdirs/acme/inputdata/atm/cam/sst/sst_HadOIBl_bc_1x1_clim_c101029.nc')
        # # ds.drop('time')
        # ds = ds.assign_coords(coords={'time':ds2['time'][0:2]})
        # # ds['time'] = ds2['time'][0:1]
        # if 'time' not in ds.dims : ds = ds.expand_dims('time',axis=0)
        # if 'time' not in ds.coords : ds = ds.assign_coords(coords={'time':ds['time']})
        # ds['date'] = ds2['date'][0:2]
        # ds['datesec'] = ds2['datesec'][0:2]
        # ds2.close()
        # ----------------------------------------------------------------------
        # ----------------------------------------------------------------------

        # Write back to final file
        ds.to_netcdf(output_file_name
                    ,unlimited_dims=['time'] 
                    ,encoding={'time':{'dtype':'float64'}}
                    ,format='NETCDF4_CLASSIC' )
        ds.close()

        run_cmd(f'ncatted --hst -A -a calendar,time,m,c,\'365_day\' {output_file_name}',
                verbose,prepend_line=False,shell=True)

        run_cmd(f'ncatted -O -a _FillValue,time,d,, {output_file_name}',
                verbose,prepend_line=False)

        return
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Subclasses
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
class ERA5(hiccup_data):
    @classmethod
    def is_name_for(cls,name) : return name == 'ERA5'
    def __init__(self,name,atm_file,sfc_file,dst_horz_grid,dst_vert_grid,
                 output_dir=default_output_dir,grid_dir=default_grid_dir,
                 map_dir=default_map_dir,tmp_dir=default_tmp_dir,
                 sstice_name=None,sst_file=None,ice_file=None,
                 sstice_combined_file=None,lev_type=''):
        super().__init__(atm_file=atm_file
                        ,sfc_file=sfc_file
                        ,dst_horz_grid=dst_horz_grid
                        ,dst_vert_grid=dst_vert_grid
                        ,sstice_name=sstice_name
                        ,sst_file=sst_file
                        ,ice_file=ice_file
                        ,sstice_combined_file=sstice_combined_file
                        ,output_dir=output_dir
                        ,grid_dir=grid_dir
                        ,map_dir=map_dir
                        ,tmp_dir=tmp_dir
                        ,lev_type=lev_type)
        
        self.name = 'ERA5'
        self.lev_name = 'level'

        # Atmospheric variables
        self.atm_var_name_dict.update({'lat':'latitude'})
        self.atm_var_name_dict.update({'lon':'longitude'})
        self.atm_var_name_dict.update({'T':'t'})            # temperature
        self.atm_var_name_dict.update({'Q':'q'})            # specific humidity
        # self.atm_var_name_dict.update({'Z3':'z'})           # geopotential (not sure we need this)
        self.atm_var_name_dict.update({'U':'u'})            # zonal wind
        self.atm_var_name_dict.update({'V':'v'})            # meridional wind 
        self.atm_var_name_dict.update({'CLDLIQ':'clwc'})    # specific cloud liq water 
        self.atm_var_name_dict.update({'CLDICE':'ciwc'})    # specific cloud ice water 
        self.atm_var_name_dict.update({'O3':'o3'})          # ozone mass mixing ratio 

        # Surface variables
        self.sfc_var_name_dict.update({'PS':'sp'})         # sfc pressure 
        self.sfc_var_name_dict.update({'TS':'skt'})        # skin temperature 
        self.sfc_var_name_dict.update({'PHIS':'z'})        # surface geopotential
        # self.sfc_var_name_dict.update({'SST':'sst'})       # sea sfc temperature 
        self.sfc_var_name_dict.update({'TS1':'stl1'})      # Soil temperature level 1 
        self.sfc_var_name_dict.update({'TS2':'stl2'})      # Soil temperature level 2 
        self.sfc_var_name_dict.update({'TS3':'stl3'})      # Soil temperature level 3 
        self.sfc_var_name_dict.update({'TS4':'stl4'})      # Soil temperature level 4 
        self.sfc_var_name_dict.update({'ICEFRAC':'siconc'})    # Sea ice area fraction
        self.sfc_var_name_dict.update({'SNOWHICE':'sd'})     # Snow depth 
        # self.sfc_var_name_dict.update({'':'asn'})          # Snow albedo 
        # self.sfc_var_name_dict.update({'':'rsn'})          # Snow density 
        # self.sfc_var_name_dict.update({'':'tsn'})          # Temperature of snow layer 
        # self.sfc_var_name_dict.update({'':'lai_hv'})       # Leaf area index, high vegetation 
        # self.sfc_var_name_dict.update({'':'lai_lv'})       # Leaf area index, low vegetation 
        # self.sfc_var_name_dict.update({'':'src'})          # Skin reservoir content 
        # self.sfc_var_name_dict.update({'':'swvl1'})        # Volumetric soil water level 1 
        # self.sfc_var_name_dict.update({'':'swvl2'})        # Volumetric soil water level 2 
        # self.sfc_var_name_dict.update({'':'swvl3'})        # Volumetric soil water level 3 
        # self.sfc_var_name_dict.update({'':'swvl4'})        # Volumetric soil water level 4 

        self.nlat = len( self.ds_atm[ self.atm_var_name_dict['lat'] ].values )
        self.nlon = len( self.ds_atm[ self.atm_var_name_dict['lon'] ].values )

        self.src_grid_name = f'{self.nlat}x{self.nlon}'
        self.src_grid_file = self.grid_dir+f'scrip_{self.name}_{self.src_grid_name}.nc'

        self.map_file = self.map_dir+f'map_{self.src_grid_name}_to_{self.dst_horz_grid}.nc'

    # --------------------------------------------------------------------------
    def create_src_grid_file(self,verbose=None):
        """ 
        Generate source grid file 
        """
        if verbose is None : verbose = hiccup_verbose
        if verbose : print('\nGenerating src grid file...')

        # Remove the file here to prevent the warning message when ncremap overwrites it
        if self.src_grid_file in glob.glob(self.src_grid_file) : 
            run_cmd(f'rm {self.src_grid_file} ',verbose)

        check_dependency('ncremap')

        cmd  = f'ncremap {ncremap_alg} ' 
        cmd += f' --tmp_dir={self.tmp_dir}'
        cmd += f' -G ttl=\'Equi-Angular grid {self.src_grid_name}\'' 
        cmd += f'#latlon={self.nlat},{self.nlon}'                    
        cmd +=  '#lat_typ=uni'
        cmd +=  '#lat_drc=n2s'
        cmd +=  '#lon_typ=grn_ctr '
        cmd += f' -g {self.src_grid_file} '
        run_cmd(cmd,verbose,shell=True)

        return 
    # --------------------------------------------------------------------------
    def rename_vars_special(self,file_name,verbose=None):
        """ 
        Rename file vars specific to this subclass 
        """
        if verbose is None : verbose = hiccup_verbose
        
        new_lev_name = 'plev'

        check_dependency('ncrename')
        check_dependency('ncap2')
        check_dependency('ncatted')
        check_dependency('ncks')

        # Rename pressure variable (needed for vertical remap)
        run_cmd(f'ncrename -d {self.lev_name},{new_lev_name} -v level,{new_lev_name} {file_name}',
            verbose,shell=True)

        # Reset the level variable name 
        self.lev_name = new_lev_name

        # change pressure variable type to double and units to Pascals (needed for vertical remap)
        run_cmd(f"ncap2 -O -s '{new_lev_name}={new_lev_name}.convert(NC_DOUBLE)*100' {file_name} {file_name}",
                verbose,prepend_line=False,shell=True)

        # change units attribute
        run_cmd(f"ncatted --hst -A -a units,{new_lev_name},a,c,'Pa' {file_name}",
                verbose,prepend_line=False,shell=True)

        # Remove lat/lon vertices variables since they are not needed
        run_cmd(f'ncks -C -O  -x -v lat_vertices,lon_vertices {file_name} {file_name}',
                verbose,prepend_line=False,shell=True)

        # also remove "bounds" attribute
        run_cmd(f'ncatted -O -a bounds,lat,d,, {file_name} {file_name}',
                verbose,prepend_line=False,shell=True)
        run_cmd(f'ncatted -O -a bounds,lon,d,, {file_name} {file_name}',
                verbose,shell=True)

        return
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
