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
import os, sys, re, shutil
from time import perf_counter
from hiccup import hiccup_state_adjustment as hsa
from hiccup import hiccup_utilities as hu

# default output paths
default_output_dir  = './data/'
default_grid_dir    = './files_grid/'
default_map_dir     = './files_mapping/'
default_tmp_dir     = './files_tmp'

# algorithm flag for ncremap
ncremap_alg         = ' --alg_typ=tempest '    

# log file for Tempest output
tempest_log_file    = 'TempestRemap.log'

# override the xarray default netcdf format of 
# NETCDF4 to avoid file permission issue
# NETCDF4 / NETCDF4_CLASSIC / NETCDF3_64BIT
hiccup_sst_nc_format = 'NETCDF3_64BIT'
hiccup_atm_nc_format = 'NETCDF4'

# set the ncremap file type 
# netcdf4 / netcdf4_classic / 64bit_data / 64bit_offset
ncremap_file_fmt = '64bit_data' 

# use 100k header padding for improved performance when editing metadata
# see http://nco.sourceforge.net/nco.html#hdr_pad
hdr_pad = 100000

# Global verbosity default
hiccup_verbose = False

# Set numpy to ignore overflow errors
np.seterr(over='ignore')

# Disable HDF file locking to prevent permission 
# errors when writing data to files
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

# Enable timers by default
do_timers = True
timer_msg_all = []
timer_start_total = None

# Numeric parameters
tk_zero = 273.15 # value for converting between celsius and Kelvin

# Set up terminal colors
class tcolor:
    """ simple class for coloring terminal text """
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
# Print individual timer information
# ------------------------------------------------------------------------------
def print_timer(timer_start,use_color=True,prefix='\n',caller=None,print_msg=True):
    """
    Print the final timer result based on input start time
    Also update timer_msg_all for use in print_timer_summary
    """
    # if caller is not provider get name of parent routine
    if caller is None: caller = sys._getframe(1).f_code.co_name
    # calculate elapsed time
    etime = perf_counter()-timer_start
    time_str = f'{etime:10.1f} sec'
    # add minutes if longer than 60 sec or 2 hours
    if etime>60       : time_str += f' ({(etime/60):4.1f} min)'
    # if etime>(2*3600) : time_str += f' ({(etime/3600):.1f} hr)'
    # create the timer result message
    msg = f'{caller:35} elapsed time: {time_str}'
    # add message to list of messages for print_timer_summary
    timer_msg_all.append(msg)
    # Apply color
    if use_color : msg = tcolor.YELLOW + msg + tcolor.ENDC
    # print the message
    print(prefix+msg)
    return
# ------------------------------------------------------------------------------
# Print a summary of timer information
# ------------------------------------------------------------------------------
def print_timer_summary():
    """
    Print timer summary based on information compiled by print_timer()
    """
    # Add timer info for all if timer_start_total was set
    if timer_start_total is not None: 
        print_timer(timer_start_total,caller=f'Total',print_msg=False)
    if do_timers:
        print('\nHICCUP Timer results:')
        for msg in timer_msg_all:
            print(f'  {msg}')
    return
# ------------------------------------------------------------------------------
# Get default topography file name
# ------------------------------------------------------------------------------
def get_default_topo_file_name(grid,topo_file_root=None):
    """
    return default topo file associated with input grid name
    """
    root_err_msg = 'No default for topo_file_root! Topo file root must be manually specified.'
    file_err_msg = 'No default for topo_file_name! Topo file path must be manually specified.'

    # Set root directory path if not provided
    if topo_file_root is None:
        nersc_inputdata_path = '/global/cfs/projectdirs/e3sm/inputdata'
        olcf_inputdata_path  = '/gpfs/alpine/cli115/world-shared/e3sm/inputdata'
        if os.path.exists(nersc_inputdata_path): topo_file_root = f'{nersc_inputdata_path}/atm/cam/topo'
        if os.path.exists(olcf_inputdata_path) : topo_file_root = f'{ olcf_inputdata_path}/atm/cam/topo'
        if topo_file_root is None: raise ValueError(root_err_msg)
    
    # Set default topo file
    topo_file_name = None
    if grid=='ne1024np4': topo_file_name = f'{topo_file_root}/USGS-gtopo30_ne1024np4_16xconsistentSGH_20190528.nc'
    if grid=='ne256np4' : topo_file_name = f'{topo_file_root}/USGS-gtopo30_ne256np4pg2_16xdel2_20200213.nc'
    if grid=='ne120np4' : topo_file_name = f'{topo_file_root}/USGS-gtopo30_ne120np4_16xdel2-PFC-consistentSGH.nc'
    if grid=='ne45np4'  : topo_file_name = f'{topo_file_root}/USGS-gtopo30_ne45np4pg2_16xdel2.c20200615.nc'
    if grid=='ne30np4'  : topo_file_name = f'{topo_file_root}/USGS-gtopo30_ne30np4pg2_16xdel2.c20200108.nc'
    if grid=='ne16np4'  : topo_file_name = f'{topo_file_root}/USGS-gtopo30_ne16np4pg2_16xdel2_20200527.nc'
    if grid=='ne11np4'  : topo_file_name = f'{topo_file_root}/USGS-gtopo30_ne11np4_16xconsistentSGH.c20160612.nc'
    if grid=='ne4np4'   : topo_file_name = f'{topo_file_root}/USGS-gtopo30_ne4pg2_16xdel2-PFC-consistentSGH.c20190618.nc'
    
    if topo_file_name is None: raise ValueError(file_err_msg)
    return topo_file_name
# ------------------------------------------------------------------------------
# Method for returning class object
# ------------------------------------------------------------------------------
def create_hiccup_data(name,atm_file,sfc_file,dst_horz_grid,dst_vert_grid,
                       output_dir=default_output_dir,grid_dir=default_grid_dir,
                       map_dir=default_map_dir,tmp_dir=default_tmp_dir,
                       sstice_combined_file=None,sstice_name=None,
                       sst_file=None,ice_file=None,topo_file=None,
                       lev_type='',verbose=False,
                       check_input_files=True):
    """ 
    Create HICCUP data class object, check for required input variables and 
    create specified output directories if they do not exist
    """
    global hiccup_verbose
    hiccup_verbose = verbose
    hu.check_nco_version()
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
                      ,topo_file=topo_file
                      ,dst_horz_grid=dst_horz_grid
                      ,dst_vert_grid=dst_vert_grid
                      ,output_dir=output_dir
                      ,grid_dir=grid_dir
                      ,map_dir=map_dir
                      ,tmp_dir=tmp_dir
                      ,lev_type=lev_type)

            # Check input files for for required variables
            if check_input_files: obj.check_file_vars()

            # Create the output, grid, and map folders if they do not exist
            if not os.path.exists(output_dir) : os.makedirs(output_dir)
            if not os.path.exists(grid_dir)   : os.makedirs(grid_dir)
            if not os.path.exists(map_dir)    : os.makedirs(map_dir)

            global timer_start_total
            timer_start_total = perf_counter()

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
                 sstice_combined_file=None,sstice_name=None,topo_file=None,
                 sst_file=None,ice_file=None,lev_type=''):
        self.lev_type = lev_type
        self.atm_file = atm_file
        self.sfc_file = sfc_file
        self.topo_file = topo_file
        self.atm_var_name_dict = {}
        self.sfc_var_name_dict = {}
        self.lnd_var_name_dict = {}
        self.src_nlat = -1
        self.src_nlon = -1
        self.dst_horz_grid = dst_horz_grid
        self.dst_vert_grid = dst_vert_grid
        self.src_horz_grid = ''
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

        # Check that input files exist
        for file_name in [self.atm_file,self.sfc_file,self.sst_file
                         ,self.ice_file,self.topo_file]:
            if file_name is not None:
                if not os.path.exists(file_name):
                    raise ValueError(f'input file does not exist: {file_name}')

        # Load input files into xarray datasets
        self.ds_atm = xr.open_dataset(self.atm_file)
        self.ds_sfc = xr.open_dataset(self.sfc_file)
    # --------------------------------------------------------------------------
    def __str__(self):
        indent = '    '
        str_out = '\nHICCUP data object:\n'
        fmt_key_len = 18
        for key in self.__dict__.keys(): 
            attribute = getattr(self,key)
            if isinstance(attribute,dict):
                str_out += f'{indent}{key:{fmt_key_len}}:\n'
                for k in attribute.keys(): 
                    str_out += f'{indent}'*(fmt_key_len+4)+f'{k:8}  {attribute[k]}\n'
            # elif isinstance(attribute, xr.Dataset) : 
            #     # Print details of xarray dataset
            #     str_out += f'  {key:{fmt_key_len}}:'
            #     ds_str = attribute.__str__().replace('\n','\n'+' '*(fmt_key_len+4))
            #     str_out += f'  {ds_str}\n'
            else:
                if attribute!='' and not isinstance(attribute, xr.Dataset) :
                    str_out += f'{indent}{key:{fmt_key_len}}:  {attribute}\n'

        return str_out
    # --------------------------------------------------------------------------
    def get_src_grid_ne(self):
        """
        Return number of elements of source grid (if starting from model data)
        """
        if hasattr(self, 'src_horz_grid_np'):
            result = re.search('ne(.*)np', self.src_horz_grid_np)
            return result.group(1) if result else 0
        else:
            raise AttributeError('src_horz_grid_np not found in HICCUP object')
    # --------------------------------------------------------------------------
    def get_src_grid_npg(self):
        """
        Return number of FV physgrid cells (npg) of source grid (if starting from model data)
        """
        if hasattr(self, 'src_horz_grid_pg'):
            result = re.search('pg(.*)', self.src_horz_grid_pg)
            return result.group(1) if result else 0
        else:
            raise AttributeError('src_horz_grid_pg not found!')
    # --------------------------------------------------------------------------
    def get_dst_grid_ne(self):
        """
        Return number of elements of target model grid
        """
        result = re.search('ne(.*)np', self.dst_horz_grid)
        return result.group(1) if result else 0
    # --------------------------------------------------------------------------
    def get_dst_grid_npg(self):
        """
        Return number of FV physgrid cells (npg) of target model grid
        """
        if hasattr(self, 'dst_horz_grid_pg'):
            result = re.search('pg(.*)', self.dst_horz_grid_pg)
            return result.group(1) if result else 0
        else:
            result = re.search('pg(.*)', self.dst_horz_grid)
            return result.group(1) if result else 0
    # --------------------------------------------------------------------------
    def get_chunks(self):
        """
        Return chunk number to use dask for certain special cases
        """
        # By default we do not want to use chunk, but for very fine grids
        # it is useful to load into a dask array by setting the chunk size
        chunks = None
        ne,npg = int(self.get_dst_grid_ne()),int(self.get_dst_grid_npg())
        # divide total physics column count by 2 for large grids
        if ne>120 and npg==0: chunks = {'ncol':int((ne*ne*54+2)/2)} 
        if ne>120 and npg>0 : chunks = {'ncol':int((ne*ne*6*npg)/2)} 
        return chunks
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
        if do_timers: timer_start = perf_counter()
        if verbose is None : verbose = hiccup_verbose
        if verbose : print('\nUnpacking data files...')

        hu.check_dependency('ncpdq')

        for f in [ self.atm_file, self.sfc_file, 
                   self.sst_file, self.ice_file,
                   self.sstice_combined_file ]:
            if f is not None :
                run_cmd(f'ncpdq -U --ovr {f} {f}',verbose,prepend_line=False)
        if do_timers: print_timer(timer_start)
        return
    # --------------------------------------------------------------------------
    def create_dst_grid_file(self,verbose=None):
        """ 
        Generate destination model grid file. Normally, we only care about 
        mapping to the GLL/np4 grid, unless the source data is an EAM file
        """
        if do_timers: timer_start = perf_counter()
        if verbose is None : verbose = hiccup_verbose
        if verbose : print('\nGenerating dst grid file...')

        if 'ne' in self.dst_horz_grid and 'np' in self.dst_horz_grid : 
            
            # Spectral element grid with physics on GLL nodes
            ne = self.get_dst_grid_ne()
            self.dst_grid_file = f'{self.grid_dir}/exodus_ne{ne}.g'
            
            hu.check_dependency('GenerateCSMesh')
            cmd = f'GenerateCSMesh --alt --res {ne} --file {self.dst_grid_file}'
            cmd += f' >> {tempest_log_file}'
            run_cmd(cmd,verbose,shell=True)

        else:
            raise ValueError(f'dst_horz_grid={self.dst_horz_grid} is not currently supported')

        if do_timers: print_timer(timer_start)
        return 
    # --------------------------------------------------------------------------
    def create_map_file(self,verbose=None,src_type=None,dst_type=None):
        """ 
        Generate mapping file after grid files have been created.
        This routine assumes that the destination is always GLL/np4.
        For mapping EAM to EAM data this method is overloaded below. 
        """
        if do_timers: timer_start = perf_counter()
        if verbose is None : verbose = hiccup_verbose
        if verbose : print('\nGenerating mapping file...')

        hu.check_dependency('ncremap')

        # Check that grid file fields are not empty
        if self.src_grid_file is None : raise ValueError('src_grid_file is not defined!')
        if self.dst_grid_file is None : raise ValueError('dst_grid_file is not defined!')

        if src_type is None: src_type = 'FV' # assume input is FV
        if dst_type is None: dst_type = 'GLL' # assume dst grid is GLL/np4
        
        ne = self.get_dst_grid_ne()

        # Set the map options (do we need the --mono flag?)
        self.map_opts = ''
        self.map_opts = self.map_opts+' --in_type cgll --in_np 4 ' 
        self.map_opts = self.map_opts+' --out_type cgll --out_np 4 '
        self.map_opts = self.map_opts+' --out_double '
        
        # Create the map file
        cmd = f'ncremap {ncremap_alg} '
        cmd += f' --src_grd={self.src_grid_file}'
        cmd += f' --dst_grd={self.dst_grid_file}'
        cmd += f' --map_file={self.map_file}'
        cmd += f' --wgt_opt=\'{self.map_opts}\' '
        # Add special flag for "very fine" grids
        # (this assumes that source data has ~0.5 degree spacing like ERA5)
        if int(ne)>100 : cmd += ' --lrg2sml ' 
        run_cmd(cmd,verbose,shell=True)

        if do_timers: print_timer(timer_start)
        return
    # --------------------------------------------------------------------------
    def get_multifile_dict(self,verbose=None,timestamp=None):
        """
        Create dict of temporary file names associated with each variable
        """
        if do_timers: timer_start = perf_counter()
        if verbose is None : verbose = hiccup_verbose
        if verbose : print('\nCreating list of temporary files...')

        # define file list to be returned
        tmp_file_dict = {}

        if 'lat' in self.atm_var_name_dict: lat_var = self.atm_var_name_dict['lat']
        if 'lon' in self.atm_var_name_dict: lon_var = self.atm_var_name_dict['lon']

        var_dict_all = self.atm_var_name_dict.copy()
        var_dict_all.update(self.sfc_var_name_dict)

        # use timestamp to ensure these files are distinct from 
        # other instances of HICCUP that might be running concurrently
        if timestamp is None: timestamp = datetime.datetime.utcnow().strftime('%Y%m%d.%H%M%S')

        # Horzontally remap atmospher and surface data to individual files
        for key,var in var_dict_all.items() :
            if var not in [lat_var,lon_var]:
                tmp_file_name = None
                if key in self.sfc_var_name_dict.keys(): tmp_file_name = f'{self.tmp_dir}tmp_sfc_data'
                if key in self.atm_var_name_dict.keys(): tmp_file_name = f'{self.tmp_dir}tmp_atm_data'
                if tmp_file_name is not None: 
                    tmp_file_name += f'.{self.dst_horz_grid}'
                    tmp_file_name += f'.{self.dst_vert_grid}'
                    tmp_file_name += f'.{key}'
                    tmp_file_name += f'.{timestamp}'
                    tmp_file_name += f'.nc'
                    tmp_file_dict.update({key:tmp_file_name})
                    if verbose: print(f'  {key:10}   {tmp_file_name}')

        return tmp_file_dict
    # --------------------------------------------------------------------------
    def get_multifile_dict_eam(self,var_name_dict,verbose=None,timestamp=None,dst_horz_grid=None):
        """
        Create dict of temporary file names associated with each variable
        """
        if do_timers: timer_start = perf_counter()
        if verbose is None : verbose = hiccup_verbose
        if verbose : print('\nCreating list of temporary files...')

        # define file list to be returned
        tmp_file_dict = {}

        if dst_horz_grid is None: dst_horz_grid  = self.dst_horz_grid

        # if 'lat_d' in self.atm_var_name_dict: lat_var = self.atm_var_name_dict['lat_d']
        # if 'lon_d' in self.atm_var_name_dict: lon_var = self.atm_var_name_dict['lon_d']
        # if 'lat'   in self.atm_var_name_dict: lat_var = self.atm_var_name_dict['lat']
        # if 'lon'   in self.atm_var_name_dict: lon_var = self.atm_var_name_dict['lon']

        # use timestamp to ensure these files are distinct from 
        # other instances of HICCUP that might be running concurrently
        if timestamp is None: timestamp = datetime.datetime.utcnow().strftime('%Y%m%d.%H%M%S')

        # Horzontally remap atmospher and surface data to individual files
        for key,var in var_name_dict.items() :
            # if var not in [lat_var,lon_var]:
            tmp_file_name = f'{self.tmp_dir}/tmp_data'
            tmp_file_name += f'.{dst_horz_grid}'
            tmp_file_name += f'.{self.dst_vert_grid}'
            tmp_file_name += f'.{key}'
            tmp_file_name += f'.{timestamp}'
            tmp_file_name += f'.nc'
            tmp_file_dict.update({key:tmp_file_name})
            if verbose: print(f'  {key:10}   {tmp_file_name}')

        return tmp_file_dict
    # --------------------------------------------------------------------------
    def rename_vars(self,file_name,verbose=None):
        """ 
        Rename variables in file according to variable name dictionaries 
        """
        if do_timers: timer_start = perf_counter()
        if verbose is None : verbose = hiccup_verbose
        if verbose : print('\nRenaming variables to match model variable names...')

        hu.check_dependency('ncrename')

        # Alternate approach - build a single large command to rename all at once
        var_dict_all = self.atm_var_name_dict.copy()
        var_dict_all.update(self.sfc_var_name_dict)
        cmd = f'ncrename --hst'
        for key in var_dict_all : 
            if key != var_dict_all[key]:
                tmp_cmd = f' -v {var_dict_all[key]},{key} '
                if tmp_cmd not in cmd :
                    # coords are often already renamed by the remapping step, 
                    # so make them optional by adding a preceeding dot
                    if key in ['lat','lon']: 
                        tmp_cmd = tmp_cmd.replace(f'{var_dict_all[key]}',f'.{var_dict_all[key]}')
                    cmd += tmp_cmd
        run_cmd(f'{cmd} {file_name}',verbose,prepend_line=False,shell=True)

        # Stop timer here before calling rename_vars_special
        if do_timers: print_timer(timer_start)

        # Do additional variable/attribute renaming specific to the input data
        self.rename_vars_special(file_name,verbose)

        return
    # --------------------------------------------------------------------------
    def rename_vars_multifile(self,file_dict,verbose=None):
        """ 
        Rename variables in file list according to variable name dictionaries 
        This approach was developed specifically for very fine grids like ne1024
        """
        if do_timers: timer_start = perf_counter()
        if verbose is None : verbose = hiccup_verbose
        if verbose : print('\nRenaming variables to match model variable names...')

        hu.check_dependency('ncrename')

        if 'lat' in self.atm_var_name_dict: lat_var = self.atm_var_name_dict['lat']
        if 'lon' in self.atm_var_name_dict: lon_var = self.atm_var_name_dict['lon']

        var_dict_all = self.atm_var_name_dict.copy()
        var_dict_all.update(self.sfc_var_name_dict)

        new_lev_name = None
        if self.lev_name=='level': new_lev_name = 'plev'

        for var, file_name in file_dict.items():

            with xr.open_dataset(file_name) as ds_data:
                ds_data.load()
                # Rename the variable
                ds_data = ds_data.rename({var_dict_all[var]:var})
                # Rename the vertical coordinate
                if new_lev_name is not None and '_sfc_' not in file_name:
                    ds_data = ds_data.rename({self.lev_name:new_lev_name})
                # Do additional variable/attribute renaming specific to the input data
                adjust_pressure_units = False
                if new_lev_name is not None and '_sfc_' not in file_name: adjust_pressure_units = True
                self.rename_vars_special(ds_data,verbose,do_timers=False
                                        ,adjust_pressure_units=adjust_pressure_units
                                        ,change_pressure_name=False
                                        ,new_lev_name=new_lev_name)
            ds_data.to_netcdf(file_name,format=hiccup_atm_nc_format,mode='w')
            ds_data.close()

        # Reset the level variable name 
        if new_lev_name is not None: self.lev_name = new_lev_name

        if do_timers: print_timer(timer_start)

        return
    # --------------------------------------------------------------------------
    def add_reference_pressure(self,file_name,verbose=None):
        """ 
        Add P0 variable 
        """
        if do_timers: timer_start = perf_counter()
        if verbose is None : verbose = hiccup_verbose
        if verbose : print('\nAdding reference pressure (P0)...')

        hu.check_dependency('ncap2')
        hu.check_dependency('ncatted')

        # Add the variable
        run_cmd(f"ncap2 --hst -A -s 'P0=100000.' {file_name} {file_name}",
                verbose,prepend_line=False,shell=True)

        # add long_name and units attributes
        run_cmd(f"ncatted --hst -A -a long_name,P0,a,c,'reference pressure' -a units,P0,a,c,'Pa' {file_name}",
                verbose,prepend_line=False,shell=True)
        
        if do_timers: print_timer(timer_start)
        return
    # --------------------------------------------------------------------------
    def remap_horizontal(self,output_file_name,verbose=None):
        """  
        Horizontally remap data and combine into single file 
        """
        if do_timers: timer_start = perf_counter()
        if verbose is None : verbose = hiccup_verbose
        if verbose : print('\nHorizontally remapping the data to temporary files...')

        if self.map_file is None: raise ValueError('map_file cannot be None!')
        if self.atm_file is None: raise ValueError('atm_file cannot be None!')
        if self.sfc_file is None: raise ValueError('sfc_file cannot be None!')

        # Define temporary files that will be deleted at the end
        atm_tmp_file_name = f'{self.tmp_dir}tmp_atm_data.nc'
        sfc_tmp_file_name = f'{self.tmp_dir}tmp_sfc_data.nc'

        # Remove temporary files if they exist
        if os.path.isfile(atm_tmp_file_name): run_cmd(f'rm {atm_tmp_file_name} ',verbose)
        if os.path.isfile(sfc_tmp_file_name): run_cmd(f'rm {sfc_tmp_file_name} ',verbose)

        hu.check_dependency('ncremap')
        hu.check_dependency('ncks')

        # Horzontally remap atmosphere data
        var_list = ','.join(self.atm_var_name_dict.values())
        cmd =  f'ncremap {ncremap_alg} '
        cmd += f' --map_file={self.map_file} '
        cmd += f' --in_file={self.atm_file} '
        cmd += f' --out_file={atm_tmp_file_name} '
        cmd += f' --var_lst={var_list} '
        cmd += f' --fl_fmt={ncremap_file_fmt} '
        run_cmd(cmd,verbose)

        # Horzontally remap surface data
        var_list = ','.join(self.sfc_var_name_dict.values())
        cmd =  f'ncremap {ncremap_alg} '
        cmd += f' --map_file={self.map_file} '
        cmd += f' --in_file={self.sfc_file} '
        cmd += f' --out_file={sfc_tmp_file_name} '
        cmd += f' --var_lst={var_list} '
        cmd += f' --fl_fmt={ncremap_file_fmt} '
        run_cmd(cmd,verbose)

        # Remove output file if it already exists
        if os.path.isfile(output_file_name): run_cmd(f'rm {output_file_name} ',verbose)

        if verbose : print('\nCombining temporary remapped files...')

        # Add atmosphere temporary file data into the final output file
        run_cmd(f'ncks -A --hdr_pad={hdr_pad} {atm_tmp_file_name} {output_file_name} ',
                verbose,prepend_line=False)

        # Add surface temporary file data into the final output file
        run_cmd(f'ncks -A --hdr_pad={hdr_pad} {sfc_tmp_file_name} {output_file_name} ',
                verbose,prepend_line=False)

        # delete the temporary files
        run_cmd(f'rm {sfc_tmp_file_name} {atm_tmp_file_name} ',verbose,prepend_line=False)

        if do_timers: print_timer(timer_start)
        return
    # --------------------------------------------------------------------------
    def remap_horizontal_multifile(self,file_dict,verbose=None):
        """  
        Horizontally remap data into seperate files for each variable
        This approach was developed specifically for very fine grids like ne1024
        """
        if do_timers: timer_start = perf_counter()
        if verbose is None : verbose = hiccup_verbose
        if verbose : print('\nHorizontally remapping the multi-file data to temporary files...')

        if self.map_file is None: raise ValueError('map_file cannot be None!')
        if self.atm_file is None: raise ValueError('atm_file cannot be None!')
        if self.sfc_file is None: raise ValueError('sfc_file cannot be None!')

        hu.check_dependency('ncremap')

        if 'lat' in self.atm_var_name_dict: lat_var = self.atm_var_name_dict['lat']
        if 'lon' in self.atm_var_name_dict: lon_var = self.atm_var_name_dict['lon']

        # Horzontally remap atmosphere and surface data to individual files
        for var,tmp_file_name in file_dict.items():
            if var in self.sfc_var_name_dict.keys(): 
                in_var = self.sfc_var_name_dict[var]
                in_file = self.sfc_file
            if var in self.atm_var_name_dict.keys(): 
                in_var = self.atm_var_name_dict[var]
                in_file = self.atm_file
            # Remove temporary files if they exist
            if os.path.isfile(tmp_file_name): run_cmd(f'rm {tmp_file_name}',verbose)
            # Remap the data
            in_var_list = f'{in_var}'
            if 'lat' in self.atm_var_name_dict: in_var_list += f',{lat_var}'
            if 'lon' in self.atm_var_name_dict: in_var_list += f',{lon_var}'
            cmd  = f'ncremap {ncremap_alg} '
            cmd += f" --nco_opt='-O --no_tmp_fl --hdr_pad={hdr_pad}' "
            cmd += f' --map_file={self.map_file}'
            cmd += f' --in_file={in_file}'
            cmd += f' --out_file={tmp_file_name}'
            # cmd += f' --var_lst={in_var},{lat_var},{lon_var}'
            cmd += f' --var_lst={in_var_list}'
            cmd += f' --fl_fmt={ncremap_file_fmt}'
            run_cmd(cmd,verbose,shell=True)

        if do_timers: print_timer(timer_start)

        return 
    # --------------------------------------------------------------------------
    def remap_horizontal_multifile_eam(self,file_dict,verbose=None):
        """  
        Horizontally remap data into seperate files for each variable
        This approach was developed specifically for very fine grids like ne1024
        """
        if do_timers: timer_start = perf_counter()
        if verbose is None : verbose = hiccup_verbose
        if verbose : print('\nHorizontally remapping the multi-file data to temporary files...')

        if self.map_file is None: raise ValueError('map_file cannot be None!')
        if self.atm_file is None: raise ValueError('atm_file cannot be None!')
        if self.sfc_file is None: raise ValueError('sfc_file cannot be None!')

        hu.check_dependency('ncremap')

        # Horzontally remap atmosphere and surface data to individual files
        for var,tmp_file_name in file_dict.items():
            in_var = var
            in_file = self.atm_file
            # Remove temporary files if they exist
            if os.path.isfile(tmp_file_name): run_cmd(f'rm {tmp_file_name}',verbose)
            # Remap the data
            in_var_list = f'{in_var}'
            cmd  = f'ncremap {ncremap_alg} '
            cmd += f" --nco_opt='-O --no_tmp_fl --hdr_pad={hdr_pad}' "
            cmd += f' --map_file={self.map_file}'
            cmd += f' --in_file={in_file}'
            cmd += f' --out_file={tmp_file_name}'
            cmd += f' --var_lst={in_var_list}'
            cmd += f' --fl_fmt={ncremap_file_fmt}'
            run_cmd(cmd,verbose,shell=True)

        # get rid of bounds and vertices variables
        # (normally done in rename_vars_special for non-EAM cases)
        for var,tmp_file_name in file_dict.items():
            with xr.open_dataset(tmp_file_name) as ds:
                ds.load()
                if 'bounds' in ds['lat'].attrs : del ds['lat'].attrs['bounds']
                if 'bounds' in ds['lon'].attrs : del ds['lon'].attrs['bounds']
                if 'lat_vertices' in ds.variables: ds = ds.drop('lat_vertices')
                if 'lon_vertices' in ds.variables: ds = ds.drop('lon_vertices')
            ds.to_netcdf(tmp_file_name,format=hiccup_atm_nc_format,mode='w')
            ds.close()

        if do_timers: print_timer(timer_start)

        return 
    # --------------------------------------------------------------------------
    def surface_adjustment_multifile(self,file_dict,verbose=None):
        """
        Perform surface temperature and pressure adjustments 
        using a multifile xarray dataset
        """
        if do_timers: timer_start = perf_counter()
        if verbose is None : verbose = hiccup_verbose
        if verbose: print('\nPerforming surface adjustments...')

        # update lev name in case it has not been updated previously
        self.lev_name = self.new_lev_name

        # build list of file names for variables needed for adjustment
        file_list = []
        var_list = ['TS','PS','PHIS','T']
        for var,file_name in file_dict.items():
            if var in var_list: file_list.append(file_name)

        # Load topo data for surface adjustment - use same chunking
        ds_topo = xr.open_dataset(self.topo_file,chunks=self.get_chunks())

        # Adjust surface temperature to match new surface height
        with xr.open_mfdataset(file_list,combine='by_coords',chunks=self.get_chunks()) as ds_data:
            hsa.adjust_surface_temperature( ds_data, ds_topo, verbose=verbose )
            ds_data = ds_data.compute()
        ds_data['TS'].to_netcdf(file_dict['TS'],format=hiccup_atm_nc_format,mode='w')
        ds_data.close()

        # Adjust surface pressure to match new surface height
        with xr.open_mfdataset(file_list,combine='by_coords',chunks=self.get_chunks()) as ds_data:
            hsa.adjust_surface_pressure( ds_data, ds_topo, pressure_var_name=self.lev_name
                                        ,lev_coord_name=self.lev_name, verbose=verbose )
            ds_data = ds_data.compute()
        ds_data['PS'].to_netcdf(file_dict['PS'],format=hiccup_atm_nc_format,mode='w')
        ds_data.close()

        if do_timers: print_timer(timer_start)

        return

    # --------------------------------------------------------------------------
    def remap_vertical(self,input_file_name,output_file_name,
                       vert_file_name,vert_remap_var_list=None,
                       verbose=None):
        """  
        Vertically remap data and combine into single file 
        """
        if do_timers: timer_start = perf_counter()
        if verbose is None : verbose = hiccup_verbose
        if verbose : print('\nVertically remapping the data...')

        hu.check_dependency('ncremap')

        # Specify temporary file for vertically interpolated output
        # This allows for input and output files to be the same
        if input_file_name==output_file_name:
            vert_tmp_file_name = output_file_name.replace('.nc',f'.vert_remap_tmp.nc')
        else:
            vert_tmp_file_name = output_file_name

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

        # Delete tmp file if it exists
        if os.path.isfile(vert_tmp_file_name): run_cmd(f'rm {vert_tmp_file_name} ',verbose)

        # Perform the vertical remapping
        cmd  = 'ncremap '
        cmd += f" --nco_opt='-O --no_tmp_fl --hdr_pad={hdr_pad}' " # doesn't work with vertical regridding?
        cmd += f' --vrt_fl={vert_file_name}'
        cmd += f' --var_lst={vert_remap_var_list}'
        cmd += f' --in_fl={input_file_name}'
        cmd += f' --out_fl={vert_tmp_file_name}'
        cmd += f' --fl_fmt={ncremap_file_fmt} '
        run_cmd(cmd,verbose,shell=True)

        # Overwrite the output file with the vertically interpolated data
        if input_file_name==output_file_name:
            run_cmd(f'mv {vert_tmp_file_name} {output_file_name} ')

        # The command above was causing a strange issue where the file was
        # not completely written when the next step happened. Adding the lines
        # below to fixes the problem, but I do not understand why...
        ds = xr.open_dataset(output_file_name)
        ds.close()

        if do_timers: print_timer(timer_start)
        return
    # --------------------------------------------------------------------------
    def remap_vertical_multifile(self,file_dict,vert_file_name,verbose=None):
        """
        wrapper around remap_vertical to support multi-file workflow
        specifically needed for very fine grids like ne1024
        """
        if verbose is None : verbose = hiccup_verbose
        if verbose : print('\nVertically remapping the multi-file data...')

        # temporarily disable timers and put a timer around the vertical remap loop
        global do_timers
        if do_timers: timer_start = perf_counter()
        prev_do_timers = do_timers
        do_timers = False

        ps_file_name = file_dict['PS']

        for var,file_name in file_dict.items() :
            if '_sfc_' not in file_name :
                # Append surface pressure for vertical interpolation
                run_cmd(f'ncks -A --hdr_pad={hdr_pad} {ps_file_name} {file_name}'
                        ,verbose,prepend_line=False)
                # Do the vertical interpolation for this file
                self.remap_vertical(input_file_name=file_name
                                   ,output_file_name=file_name
                                   ,vert_file_name=vert_file_name
                                   ,vert_remap_var_list=[var])
        
        # Re-set do_timers to previous value
        do_timers = prev_do_timers

        if do_timers: print_timer(timer_start)
        return
    # --------------------------------------------------------------------------
    def atmos_state_adjustment_multifile(self,file_dict,verbose=None):
        """
        Perform post-remapping atmospheric state adjustments 
        for the multifile workflow
        """
        if do_timers: timer_start = perf_counter()
        if verbose is None : verbose = hiccup_verbose
        if verbose: print('\nPerforming state adjustments...')

        # build list of file names for variables needed for adjustment
        file_list = []
        var_list = ['Q','T','PS','CLDLIQ','CLDICE']
        for var,file_name in file_dict.items():
            if var in var_list: file_list.append(file_name)
        
        with xr.open_mfdataset(file_list,combine='by_coords',chunks=self.get_chunks()) as ds_data:

            # adjust water vapor to eliminate supersaturation
            timer_start_adj = perf_counter()
            hsa.remove_supersaturation( ds_data, hybrid_lev=True, verbose=verbose )
            ds_data.compute()
            print_timer(timer_start_adj,caller='remove_supersaturation')

        ds_data['Q'].to_netcdf(file_dict['Q'],format=hiccup_atm_nc_format,mode='w')

        with xr.open_mfdataset(file_list,combine='by_coords',chunks=self.get_chunks()) as ds_data:

            # adjust cloud water to remove negative values
            timer_start_adj = perf_counter()
            hsa.adjust_cld_wtr( ds_data, verbose=verbose )
            ds_data.compute()
            print_timer(timer_start_adj,caller='adjust_cld_wtr')

        # Write adjusted data back to the individual data files
        for var in ['CLDLIQ','CLDICE']:
            if var in self.atm_var_name_dict.keys():
                ds_data[var].to_netcdf(file_dict[var],format=hiccup_atm_nc_format,mode='w')
        
        ds_data.close()

        if do_timers: print_timer(timer_start)
        return
    # --------------------------------------------------------------------------
    def add_time_date_variables(self,ds,verbose=None,do_timers=do_timers):
        """
        Check final output file and add necessary time and date information
        """
        if do_timers: timer_start = perf_counter()
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

        if do_timers: print_timer(timer_start)
        return
    # --------------------------------------------------------------------------
    def add_time_date_variables_multifile(self,file_dict,verbose=None):
        """
        """
        if do_timers: timer_start = perf_counter()
        if verbose is None : verbose = hiccup_verbose
        if verbose : print('\nEditing time and date variables...')

        for file_name in file_dict.values() :
            with xr.open_dataset(file_name) as ds_data:
                ds_data.load()
                self.add_time_date_variables(ds_data,verbose=False,do_timers=False)
            ds_data.to_netcdf(file_name,format=hiccup_atm_nc_format,mode='w')
            ds_data.close()
        if do_timers: print_timer(timer_start)
        return
    # --------------------------------------------------------------------------
    def add_time_date_variables_multifile_eam(self,file_dict,verbose=None):
        """
        copy necessary time and date information
        """
        if do_timers: timer_start = perf_counter()
        if verbose is None : verbose = hiccup_verbose
        if verbose : print('\nEditing time and date variables...')

        ds_in = xr.open_dataset(self.atm_file)

        time_vars = []
        time_vars.append('date')
        time_vars.append('datesec')
        time_vars.append('ndcur')
        time_vars.append('nscur')
        time_vars.append('nsteph')
        time_vars.append('nbdate')
        time_vars.append('ndbase')
        time_vars.append('nsbase')
        time_vars.append('nbsec')

        for file_name in file_dict.values() :
            with xr.open_dataset(file_name) as ds_data:
                ds_data.load()
                # copy time variables
                ds_data['time'] = ds_in['time']
                ds_data['time_bnds'] = ds_in['time_bnds']
                for time_var in time_vars:
                    ds_data[time_var] = ds_in[time_var]
                # # update attributes
                # ds_data['time'].attrs['calendar'] = 'noleap'
                # ds_data['time'].attrs['bounds'] = 'time_bnds'
                # ds_data['time_bnds'].attrs['long_name'] = 'time interval endpoints'

            ds_data.to_netcdf(file_name,format=hiccup_atm_nc_format,mode='w')
            ds_data.close()
        if do_timers: print_timer(timer_start)
        return
    # --------------------------------------------------------------------------
    def clean_global_attributes(self,file_name,verbose=None):
        """ 
        Remove messy global attributes of the file 
        """
        if do_timers: timer_start = perf_counter()
        if verbose is None : verbose = hiccup_verbose
        if verbose: print('\nCleaning up excessive global attributes...')
        
        global_att_list = ['history_of_appended_files', 'nco_openmp_thread_number', 
                           'input_file', 'map_file', 'remap_version', 'remap_hostname', 
                           'remap_command', 'remap_script', 'NCO' ]
        
        hu.check_dependency('ncatted')

        # Remove the attributes listed in global_att_list
        cmd = 'ncatted -O '
        for att in global_att_list: cmd += f' -a {att},global,d,, '
        cmd += f' {file_name} {file_name} '
        run_cmd(cmd,verbose)

        # Also reset the history attribute
        run_cmd(f'ncatted -h -O -a history,global,o,c, {file_name} {file_name}',
                verbose,prepend_line=False)

        if do_timers: print_timer(timer_start)
        return
    # --------------------------------------------------------------------------
    def combine_files(self,file_dict,output_file_name,delete_files=False,verbose=None):
        """
        Combine files in file_dict into single output file
        """
        if do_timers: timer_start = perf_counter()
        if verbose is None : verbose = hiccup_verbose
        if verbose: print('\nCombining temporary files into new file...')

        hu.check_dependency('ncks')

        if os.path.isfile(output_file_name): 
            run_cmd(f'rm {output_file_name} ',verbose)

        # Append each file to the output file - use netcdf4 for performance, then convert afterwards
        for var,file_name in file_dict.items() :
            # if do_timers: timer_start_combine = perf_counter()
            cmd = f'ncks -A --hdr_pad={hdr_pad} --no_tmp_fl --fl_fmt=netcdf4 {file_name} {output_file_name} '
            run_cmd(cmd,verbose,prepend_line=False)
            # if do_timers: print_timer(timer_start_combine,caller=f'  append file - {var}')

        # convert file to desired format
        if ncremap_file_fmt != 'netcdf4':
            # if do_timers: timer_start_convert = perf_counter()
            cmd = f'ncks --hdr_pad={hdr_pad} --fl_fmt={ncremap_file_fmt} {output_file_name} {output_file_name} '
            # if do_timers: print_timer(timer_start_combine,caller=f'  convert file format from netcdf4 to {ncremap_file_fmt}')

        # Delete temp files
        if delete_files:
            if verbose: print('\nDeleting temporary files...')
            for file_name in file_dict.values() :
                run_cmd(f'rm {file_name}',verbose,prepend_line=False)

        if do_timers: print_timer(timer_start)
        return
    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------
    # SST and sea ice methods
    # --------------------------------------------------------------------------
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
    def sstice_create_src_grid_file(self,diagnose_grid=True,nlat=None,nlon=None,
                                    force_overwrite=False,verbose=None):
        """
        Create a source grid file to use for remapping the SST and sea ice data.
        The SST and ice data are assumed to exist on the same grid.
        """
        if do_timers: timer_start = perf_counter()
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
        self.sstice_src_grid_file = f'{self.grid_dir}scrip_{src_grid}.nc'

        # Create the source grid file
        if force_overwrite or not os.path.isfile(self.sstice_src_grid_file) :
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

        if do_timers: print_timer(timer_start)
        return
    # --------------------------------------------------------------------------
    def open_combined_sstice_dataset(self):
        """
        Return xarray dataset of combined SST and sea ice dataset
        """
        if do_timers: timer_start = perf_counter()
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

        if do_timers: print_timer(timer_start)
        return ds_out
    # --------------------------------------------------------------------------
    def sstice_create_dst_grid_file(self,output_grid_spacing=1,force_overwrite=False,
                                    verbose=None):
        """
        Create a target grid file to use for remapping the SST and sea ice data.
        The SST and ice data are assumed to exist on the same grid.
        """
        if do_timers: timer_start = perf_counter()
        if verbose is None : verbose = hiccup_verbose

        # Define output grid dimensions
        self.sstice_nlat_dst = int( 180/output_grid_spacing )
        self.sstice_nlon_dst = int( 360/output_grid_spacing )

        dst_grid = f'{self.sstice_nlat_dst}x{self.sstice_nlon_dst}'

        # Define the destination grid file to be created
        # (add 's2n' in case input data is same grid with opposite orientation)
        self.sstice_dst_grid_file = f'{self.grid_dir}scrip_{dst_grid}_s2n.nc'

        # Create the destination grid file
        if force_overwrite or not os.path.isfile(self.sstice_dst_grid_file) :
            if verbose : print(f'\nCreating target grid file for SST and sea ice data...')
            cmd  = f'ncremap {ncremap_alg} --tmp_dir={self.tmp_dir}'
            cmd += f' -G ttl=\'Equi-Angular grid {dst_grid}\'' 
            cmd += f'#latlon={self.sstice_nlat_dst},{self.sstice_nlon_dst}'
            cmd +=  '#lat_typ=uni'
            cmd +=  '#lon_typ=grn_ctr'
            cmd +=  '#lat_drc=s2n'
            cmd += f' -g {self.sstice_dst_grid_file} '
            run_cmd(cmd,verbose,shell=True)

        if do_timers: print_timer(timer_start)
        return
    # --------------------------------------------------------------------------
    def sstice_create_map_file(self,force_overwrite=False,verbose=None):
        """
        Create a mapping file to be used for SST and sea ice data
        """
        if do_timers: timer_start = perf_counter()
        if verbose is None : verbose = hiccup_verbose

        src_grid = f'{self.sstice_nlat_src}x{self.sstice_nlon_src}'
        dst_grid = f'{self.sstice_nlat_dst}x{self.sstice_nlon_dst}'

        self.sstice_map_file = f'{self.map_dir}map_{src_grid}_to_{dst_grid}_s2n.nc'

        # Generate mapping file
        if force_overwrite or not os.path.isfile(self.sstice_map_file) :
            if verbose : print(f'\nCreating mapping file for SST and sea ice data...')
            cmd  = f'ncremap {ncremap_alg} '
            cmd +=  ' -a fv2fv '
            cmd += f' --src_grd={self.sstice_src_grid_file}'
            cmd += f' --dst_grd={self.sstice_dst_grid_file}'
            cmd += f' --map_file={self.sstice_map_file}'
            run_cmd(cmd,verbose,shell=True,prepend_line=False)

        if do_timers: print_timer(timer_start)
        return
    # --------------------------------------------------------------------------
    def sstice_slice_and_remap(self,output_file_name,
                               time_slice_method='match_atmos',
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
        if do_timers: timer_start = perf_counter()
        if verbose is None : verbose = hiccup_verbose
        if verbose : print(f'\nTime slicing {self.sstice_name} SST and sea ice data...')

        hu.check_dependency('ncatted')
        hu.check_dependency('ncremap')

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

        # make sure output file is deleted (overwrite flag not working?)
        if os.path.isfile(output_file_name): 
            run_cmd(f'rm {output_file_name}',verbose)

        # remap the SST data onto the target grid for the model
        cmd =  f'ncremap {ncremap_alg} '
        cmd += f" --nco_opt='-O --no_tmp_fl --hdr_pad={hdr_pad}' "
        cmd += f' --vars={self.sst_name},{self.ice_name} '
        cmd += f' --map_file={self.sstice_map_file} '
        cmd += f' --in_file={sstice_tmp_file_name} '
        cmd += f' --out_file={output_file_name} '
        cmd += f' --fl_fmt={ncremap_file_fmt} '
        run_cmd(cmd,verbose,shell=True)

        # delete the temporary file
        run_cmd(f'rm {sstice_tmp_file_name}',verbose)

        if do_timers: print_timer(timer_start)
        return
    # --------------------------------------------------------------------------
    def sstice_rename_vars(self, output_file_name, new_sst_name='SST_cpl',
                           new_ice_name='ice_cov', verbose=None):
        """
        Rename sst and sea icea variables and remove unnecessary variables
        """
        if do_timers: timer_start = perf_counter()
        if verbose is None : verbose = hiccup_verbose
        if verbose : print('\nRenaming SST and sea ice variables...')

        hu.check_dependency('ncrename')
        hu.check_dependency('ncks')
        hu.check_dependency('ncatted')

        # rename variables
        cmd = f'ncrename --hst'
        cmd+= f' --variable {self.sst_name},{new_sst_name}'
        cmd+= f' --variable {self.ice_name},{new_ice_name}'
        cmd+= f' {output_file_name}'
        run_cmd(cmd,verbose,shell=True)

        # Make sure dimensions names are correct
        cmd = f'ncrename --hst'
        cmd+= f' --dimension .longitude,lon'
        cmd+= f' --dimension .latitude,lat'
        cmd+= f' {output_file_name}'
        run_cmd(cmd,verbose,shell=True)

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

        if do_timers: print_timer(timer_start)
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
        if do_timers: timer_start = perf_counter()
        if verbose is None : verbose = hiccup_verbose
        if verbose : print('\nAdjusting SST and sea ice data values...')

        # Open remapped and combined data file
        ds = xr.open_dataset(output_file_name).load()

        # Convert units to Celsius
        if 'units' not in ds['SST_cpl'].attrs:
            ds['SST_cpl'].attrs['units'] = 'deg_C'
        if ds['SST_cpl'].attrs['units'] in ['K','degrees_K','deg_K','Kelvin'] :
            ds['SST_cpl'] = ds['SST_cpl']-tk_zero
            ds['SST_cpl'].attrs['units'] = 'deg_C'

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

        # reset attributes to avoid error reading file during model run
        ds['SST_cpl'].attrs = {'long_name':'Sea Surface Temperature','units':'deg_C'}
        ds['ice_cov'].attrs = {'long_name':'Sea Ice Fraction','units':'fraction'}
        ds['lat'].attrs = {'long_name':'latitude','units':'degrees_north'}
        ds['lon'].attrs = {'long_name':'longitude','units':'degrees_east'}

        # Write back to final file
        ds.to_netcdf(output_file_name
                    ,unlimited_dims=['time'] 
                    ,encoding={'time':{'dtype':'float64'}}
                    ,format=hiccup_sst_nc_format
                    )
        ds.close()

        run_cmd(f'ncatted --hst -A -a calendar,time,m,c,\'365_day\' {output_file_name}',
                verbose,prepend_line=False,shell=True)

        run_cmd(f'ncatted -O -a _FillValue,time,d,, {output_file_name}',
                verbose,prepend_line=False)

        # A confusing issue occurs with the dimension ordering 
        # the commands below provide a clunky workaround
        tmp_file = output_file_name.replace('.nc','.ncpdq_tmp.nc')
        cmd1 = f'ncpdq --permute=lon,lat,time --overwrite {output_file_name} {tmp_file}'
        cmd2 = f'ncpdq --permute=time,lat,lon --overwrite {tmp_file} {output_file_name}'
        run_cmd(cmd1,verbose)
        run_cmd(cmd2,verbose,prepend_line=False)
        run_cmd(f'rm {tmp_file} ',verbose,prepend_line=False)

        # Reset the history attribute
        run_cmd(f'ncatted -h -O -a history,global,o,c, {output_file_name} {output_file_name}',
                verbose,prepend_line=False)

        if do_timers: print_timer(timer_start)
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
                 sstice_name=None,sst_file=None,ice_file=None,topo_file=None,
                 sstice_combined_file=None,lev_type=''):
        super().__init__(atm_file=atm_file
                        ,sfc_file=sfc_file
                        ,dst_horz_grid=dst_horz_grid
                        ,dst_vert_grid=dst_vert_grid
                        ,sstice_name=sstice_name
                        ,sst_file=sst_file
                        ,ice_file=ice_file
                        ,sstice_combined_file=sstice_combined_file
                        ,topo_file = topo_file
                        ,output_dir=output_dir
                        ,grid_dir=grid_dir
                        ,map_dir=map_dir
                        ,tmp_dir=tmp_dir
                        ,lev_type=lev_type)
        
        self.name = 'ERA5'
        self.lev_name = 'level'
        self.new_lev_name = 'plev'

        # Atmospheric variables
        self.atm_var_name_dict.update({'lat':'latitude'})
        self.atm_var_name_dict.update({'lon':'longitude'})
        self.atm_var_name_dict.update({'T':'t'})            # temperature
        self.atm_var_name_dict.update({'Q':'q'})            # specific humidity
        self.atm_var_name_dict.update({'U':'u'})            # zonal wind
        self.atm_var_name_dict.update({'V':'v'})            # meridional wind 
        self.atm_var_name_dict.update({'CLDLIQ':'clwc'})    # specific cloud liq water 
        self.atm_var_name_dict.update({'CLDICE':'ciwc'})    # specific cloud ice water 
        self.atm_var_name_dict.update({'O3':'o3'})          # ozone mass mixing ratio 
        # self.atm_var_name_dict.update({'Z3':'z'})           # geopotential (not sure we need this)

        # Surface variables
        self.sfc_var_name_dict.update({'PS':'sp'})         # sfc pressure 
        self.sfc_var_name_dict.update({'TS':'skt'})        # skin temperature 
        self.sfc_var_name_dict.update({'PHIS':'z'})        # surface geopotential

        self.sfc_var_name_dict.update({'TS1':'stl1'})      # Soil temperature level 1 
        self.sfc_var_name_dict.update({'TS2':'stl2'})      # Soil temperature level 2 
        self.sfc_var_name_dict.update({'TS3':'stl3'})      # Soil temperature level 3 
        self.sfc_var_name_dict.update({'TS4':'stl4'})      # Soil temperature level 4 
        self.sfc_var_name_dict.update({'ICEFRAC':'siconc'})# Sea ice area fraction
        self.sfc_var_name_dict.update({'SNOWHICE':'sd'})   # Snow depth 
        
        # self.sfc_var_name_dict.update({'':'asn'})          # Snow albedo 
        # self.sfc_var_name_dict.update({'':'rsn'})          # Snow density 
        # self.sfc_var_name_dict.update({'':'tsn'})          # Temperature of snow layer 
        # self.sfc_var_name_dict.update({'':'lai_hv'})       # Leaf area index, high vegetation 
        # self.sfc_var_name_dict.update({'':'lai_lv'})       # Leaf area index, low vegetation 
        # self.sfc_var_name_dict.update({'':'src'})          # Skin reservoir content 
        # self.sfc_var_name_dict.update({'SST':'sst'})       # sea sfc temperature 
        # self.sfc_var_name_dict.update({'':'swvl1'})        # Volumetric soil water level 1 
        # self.sfc_var_name_dict.update({'':'swvl2'})        # Volumetric soil water level 2 
        # self.sfc_var_name_dict.update({'':'swvl3'})        # Volumetric soil water level 3 
        # self.sfc_var_name_dict.update({'':'swvl4'})        # Volumetric soil water level 4 

        self.src_nlat = len( self.ds_atm[ self.atm_var_name_dict['lat'] ].values )
        self.src_nlon = len( self.ds_atm[ self.atm_var_name_dict['lon'] ].values )

        self.src_horz_grid = f'{self.src_nlat}x{self.src_nlon}'
        self.src_grid_file = f'{self.grid_dir}/scrip_{self.name}_{self.src_horz_grid}.nc'

        self.map_file = f'{self.map_dir}/map_{self.src_horz_grid}_to_{self.dst_horz_grid}.nc'

    # --------------------------------------------------------------------------
    def create_src_grid_file(self,verbose=None):
        """ 
        Generate source grid file 
        """
        if do_timers: timer_start = perf_counter()
        if verbose is None : verbose = hiccup_verbose
        if verbose : print('\nGenerating src grid file...')

        # Remove the file here to prevent the warning message when ncremap overwrites it
        if os.path.isfile(self.src_grid_file): run_cmd(f'rm {self.src_grid_file} ',verbose)

        hu.check_dependency('ncremap')

        cmd  = f'ncremap {ncremap_alg} ' 
        cmd += f' --tmp_dir={self.tmp_dir}'
        cmd += f' -G ttl=\'Equi-Angular grid {self.src_horz_grid}\'' 
        cmd += f'#latlon={self.src_nlat},{self.src_nlon}'                    
        cmd +=  '#lat_typ=uni'
        cmd +=  '#lat_drc=n2s'
        cmd +=  '#lon_typ=grn_ctr '
        cmd += f' -g {self.src_grid_file} '
        run_cmd(cmd,verbose,shell=True)

        if do_timers: print_timer(timer_start)
        return 
    # --------------------------------------------------------------------------
    def rename_vars_special(self,ds,verbose=None,do_timers=do_timers
                           ,new_lev_name=None,change_pressure_name=True
                           ,adjust_pressure_units=True):
        """ 
        Rename file vars specific to this subclass 
        """
        if do_timers: timer_start = perf_counter()
        if verbose is None : verbose = hiccup_verbose

        if new_lev_name is None: new_lev_name = self.new_lev_name

        # Rename pressure variable and reset the lev var name (needed for vertical remap)
        if change_pressure_name:
            ds = ds.rename({self.lev_name:new_lev_name})
            self.lev_name = new_lev_name

        # change pressure variable type to double and units to Pascals (needed for vertical remap)
        if adjust_pressure_units:
            ds[new_lev_name] = ds[new_lev_name]*100.
            ds[new_lev_name].attrs['units'] = 'Pa'

        if 'bounds' in ds['lat'].attrs : del ds['lat'].attrs['bounds']
        if 'bounds' in ds['lon'].attrs : del ds['lon'].attrs['bounds']
        if 'lat_vertices' in ds.variables: ds = ds.drop('lat_vertices')
        if 'lon_vertices' in ds.variables: ds = ds.drop('lon_vertices')
        
        if do_timers: print_timer(timer_start)
        return
# ------------------------------------------------------------------------------
# subclass for remapping EAM data
# ------------------------------------------------------------------------------
class EAM(hiccup_data):
    @classmethod
    def is_name_for(cls,name) : return name == 'EAM'
    def __init__(self,name,atm_file,sfc_file,dst_horz_grid,dst_vert_grid,
                 output_dir=default_output_dir,grid_dir=default_grid_dir,
                 map_dir=default_map_dir,tmp_dir=default_tmp_dir,
                 sstice_name=None,sst_file=None,ice_file=None,topo_file=None,
                 sstice_combined_file=None,lev_type=''):
        super().__init__(atm_file=atm_file
                        ,sfc_file=sfc_file
                        ,dst_horz_grid=dst_horz_grid
                        ,dst_vert_grid=dst_vert_grid
                        ,sstice_name=sstice_name
                        ,sst_file=sst_file
                        ,ice_file=ice_file
                        ,sstice_combined_file=sstice_combined_file
                        ,topo_file = topo_file
                        ,output_dir=output_dir
                        ,grid_dir=grid_dir
                        ,map_dir=map_dir
                        ,tmp_dir=tmp_dir
                        ,lev_type=lev_type)
        
        self.name = 'EAM'
        self.lev_name = 'lev'
        self.new_lev_name = 'plev'

        self.npg = 2

        self.dst_horz_grid_np = self.dst_horz_grid
        self.dst_horz_grid_pg = self.dst_horz_grid.replace('np4',f'pg{self.npg}')

        # Determine source grid from input atmosphere file
        ds = xr.open_dataset(atm_file)
        if 'ne' in ds.attrs:
            ne = ds.attrs['ne']
        elif 'ncol' in ds.sizes:
            # use ncol formula to solve for # elements (ne): ncol_dyn = ne^2*6*9+2
            ne = int( np.sqrt( (ds.sizes['ncol']-2)/(6*9) ) )
        else:
            raise KeyError(f'Cannot determine source grid from atm_file: {atm_file}')

        self.src_horz_grid     = f'ne{ne}np4'
        self.src_horz_grid_np = f'ne{ne}np4'
        self.src_horz_grid_pg = f'ne{ne}pg{self.npg}'

        ncol_size_np = np.square(ne)*6*9+2
        ncol_size_pg = np.square(ne)*6*np.square(self.npg)

        src_ne = self.get_src_grid_ne()
        dst_ne = self.get_dst_grid_ne()

        self.src_grid_file_np = f'{self.grid_dir}/exodus_ne{src_ne}.g'
        self.src_grid_file_pg = f'{self.grid_dir}/scrip_{ self.src_horz_grid_pg}.nc'

        self.dst_grid_file_np = f'{self.grid_dir}/exodus_ne{dst_ne}.g'
        self.dst_grid_file_pg = f'{self.grid_dir}/scrip_{ self.dst_horz_grid_pg}.nc'

        self.map_file_np = f'{self.map_dir}/map_{self.src_horz_grid_np}_to_{self.dst_horz_grid_np}.nc'
        self.map_file_pg = f'{self.map_dir}/map_{self.src_horz_grid_pg}_to_{self.dst_horz_grid_pg}.nc'


        # Atmospheric variables - need separate treatment for np4 and pgN data
        self.atm_var_name_dict_np = {}
        self.atm_var_name_dict_pg = {}
        ds = xr.open_dataset(self.atm_file)
        for key in ds.variables.keys(): 
            if key in ['lat','lon','lat_d','lon_d']:
                continue
            if 'ncol_d' in ds[key].dims: 
                self.atm_var_name_dict_np.update({key:key})
            if 'ncol' in ds[key].dims:
                if ds.sizes['ncol']==ncol_size_np:
                    self.atm_var_name_dict_np.update({key:key})
                if ds.sizes['ncol']==ncol_size_pg:
                    self.atm_var_name_dict_pg.update({key:key})

    # --------------------------------------------------------------------------
    def create_src_grid_file(self,verbose=None):
        """ 
        Generate source grid file 
        """
        if do_timers: timer_start = perf_counter()
        if verbose is None : verbose = hiccup_verbose
        if verbose : print('\nGenerating src grid files (np+pg)...')

        # Remove the file here to prevent the warning message when ncremap overwrites it
        if self.src_grid_file is not None:
            if os.path.isfile(self.src_grid_file): run_cmd(f'rm {self.src_grid_file} ',verbose)

        hu.check_dependency('ncremap')

        # Spectral element grid
        ne = self.get_src_grid_ne()
        npg = self.get_src_grid_npg()

        hu.check_dependency('GenerateCSMesh')
        cmd = f'GenerateCSMesh --alt --res {ne} --file {self.src_grid_file_np}'
        cmd += f' >> {tempest_log_file}'
        run_cmd(cmd,verbose,shell=True)

        # Next switch to volumetric mesh that matches the physgrid
        tmp_exodus_file = f'{self.grid_dir}/exodus_{self.src_horz_grid_pg}.g'
        hu.check_dependency('GenerateVolumetricMesh')
        cmd = 'GenerateVolumetricMesh'
        cmd += f' --in {self.src_grid_file_np} '
        cmd += f' --out {tmp_exodus_file} '
        cmd += f' --np {npg} --uniform'
        cmd += f' >> {tempest_log_file}'
        run_cmd(cmd,verbose,shell=True)

        # Create pgN scrip file
        hu.check_dependency('ConvertExodusToSCRIP')
        scrip_file = f'{self.grid_dir}/scrip_{self.dst_horz_grid_pg}.nc'
        cmd = 'ConvertExodusToSCRIP'
        cmd += f' --in {tmp_exodus_file} '
        cmd += f' --out {self.src_grid_file_pg} '
        cmd += f' >> {tempest_log_file}'
        run_cmd(cmd,verbose,shell=True)

        # fix grid_imask type
        run_cmd(f'ncap2 --overwrite -s \'grid_imask=int(grid_imask)\' '
                +f'{self.src_grid_file_pg} {self.src_grid_file_pg}',verbose,shell=True)

        # delete temporary exodus file
        run_cmd(f'rm {tmp_exodus_file} ',verbose)

        if do_timers: print_timer(timer_start)
        return 
    # --------------------------------------------------------------------------
    def create_dst_grid_file(self,verbose=None):
        """ 
        Generate destination model grid file. For the case where EAM data 
        is the source we need both np4 and pgN grids
        """
        if do_timers: timer_start = perf_counter()
        if verbose is None : verbose = hiccup_verbose
        if verbose : print('\nGenerating dst grid files (np+pg)...')
        
        # Spectral element grid with physics on GLL nodes
        ne  = self.get_dst_grid_ne()
        npg = self.get_dst_grid_npg()

        # First create exodus file
        hu.check_dependency('GenerateCSMesh')
        cmd = f'GenerateCSMesh --alt --res {ne} --file {self.dst_grid_file_np}'
        cmd += f' >> {tempest_log_file}'
        run_cmd(cmd,verbose,shell=True)
        
        # Next switch to volumetric mesh that matches the physgrid
        tmp_exodus_file = f'{self.grid_dir}/exodus_{self.dst_horz_grid_pg}.g'
        hu.check_dependency('GenerateVolumetricMesh')
        cmd = 'GenerateVolumetricMesh'
        cmd += f' --in {self.dst_grid_file_np} '
        cmd += f' --out {tmp_exodus_file} '
        cmd += f' --np {npg} --uniform'
        cmd += f' >> {tempest_log_file}'
        run_cmd(cmd,verbose,shell=True)

        # Create scrip file while we're at it (can be slow)
        hu.check_dependency('ConvertExodusToSCRIP')
        cmd = 'ConvertExodusToSCRIP'
        cmd += f' --in {tmp_exodus_file} '
        cmd += f' --out {self.dst_grid_file_pg} '
        cmd += f' >> {tempest_log_file}'
        run_cmd(cmd,verbose,shell=True)

        # fix grid_imask type
        run_cmd('ncap2 --overwrite -s \'grid_imask=int(grid_imask)\' '
                +f'{self.dst_grid_file_pg} {self.dst_grid_file_pg}',verbose,shell=True)

        # delete temporary exodus file
        run_cmd(f'rm {tmp_exodus_file} ',verbose)

        if do_timers: print_timer(timer_start)
        return 
    # --------------------------------------------------------------------------
    def create_map_file(self,verbose=None):
        """ 
        Generate mapping files for EAM after grid files have been created 
        (overloads default routine that assumes only one map file is needed)
        """
        if do_timers: timer_start = perf_counter()
        if verbose is None : verbose = hiccup_verbose
        if verbose : print('\nGenerating mapping files (np+pg)...')

        hu.check_dependency('ncremap')

        dst_ne = self.get_dst_grid_ne()
        src_ne = self.get_src_grid_ne()

        # Check that grid file fields are not empty
        if self.src_grid_file_np is None : raise ValueError('src_grid_file_np is not defined!')
        if self.src_grid_file_pg is None : raise ValueError('src_grid_file_pg is not defined!')
        if self.dst_grid_file_np is None : raise ValueError('dst_grid_file_np is not defined!')
        if self.dst_grid_file_pg is None : raise ValueError('dst_grid_file_pg is not defined!')

        # Set the map options
        self.map_opts_np = '  --out_double  --in_type cgll --in_np 4  --out_type cgll --out_np 4 '
        self.map_opts_pg = '  --out_double  --in_type fv --in_np 2  --out_type fv --out_np 2 --volumetric '

        # Create the np4 map file
        cmd = f'ncremap {ncremap_alg} '
        cmd += f' --src_grd={self.src_grid_file_np}'
        cmd += f' --dst_grd={self.dst_grid_file_np}'
        cmd += f' --map_file={self.map_file_np}'
        cmd += f' --wgt_opt=\'{self.map_opts_np}\' '
        if dst_ne>src_ne : cmd += ' --lrg2sml '
        run_cmd(cmd,verbose,shell=True)

        # Create the pgN map file
        cmd = f'ncremap {ncremap_alg} '
        cmd += f' --src_grd={self.src_grid_file_pg}'
        cmd += f' --dst_grd={self.dst_grid_file_pg}'
        cmd += f' --map_file={self.map_file_pg}'
        cmd += f' --wgt_opt=\'{self.map_opts_pg}\' '
        if dst_ne>src_ne : cmd += ' --lrg2sml '
        run_cmd(cmd,verbose,shell=True)

        if do_timers: print_timer(timer_start)
        return
    # --------------------------------------------------------------------------
    def rename_vars_special(self,ds,verbose=None,do_timers=do_timers
                           ,new_lev_name=None,change_pressure_name=True
                           ,adjust_pressure_units=True):
        """ 
        Rename file vars specific to this subclass
        """

        return
    # --------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
