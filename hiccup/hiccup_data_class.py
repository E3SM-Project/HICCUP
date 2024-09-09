# Class structure for defining parameters specific to the 
# input data used for creating initial condition files
# 
# NOTE: Variable name dictionaries are defined with the key as the model's 
# variable name and the value as the reanalysis data variable name
# ------------------------------------------------------------------------------
import numpy as np
import xarray as xr
import pandas as pd
import datetime
import os, re
from time import perf_counter
# ------------------------------------------------------------------------------
import hiccup.hiccup_state_adjustment as hsa
from hiccup.hiccup_utilities import check_dependency
from hiccup.hiccup_utilities import run_cmd
from hiccup.hiccup_utilities import tcolor
from hiccup.hiccup_utilities import print_mem_usage
# ------------------------------------------------------------------------------
from hiccup.hiccup_constants import MW_dryair
from hiccup.hiccup_constants import MW_ozone
# ------------------------------------------------------------------------------
enable_chunks = True
ncol_chunk_size = 1000

print_memory_usage = False

# use 100k header padding for improved performance when editing metadata
# see http://nco.sourceforge.net/nco.html#hdr_pad
hdr_pad = 100000

# override the xarray default netcdf format of 
# NETCDF4 to avoid file permission issue
# NETCDF4 / NETCDF4_CLASSIC / NETCDF3_64BIT
hiccup_sst_nc_format = 'NETCDF3_64BIT'
hiccup_atm_nc_format = 'NETCDF4'

# set the ncremap file type 
# netcdf4 / netcdf4_classic / 64bit_data / 64bit_offset
ncremap_file_fmt = '64bit_data'

# log file for Tempest output
tempest_log_file = 'TempestRemap.log'

# ------------------------------------------------------------------------------
# Base Class
# ------------------------------------------------------------------------------
class hiccup_data(object):
    """ 
    Base class for HICCUP data object 
    """
    # --------------------------------------------------------------------------
    # from hiccup import hiccup_sst_methods
    def __init__( self,
                  target_model=None,
                  atm_file=None,
                  sfc_file=None,
                  dst_horz_grid=None,
                  dst_vert_grid=None,
                  output_dir=None,
                  grid_dir=None,
                  map_dir=None,
                  tmp_dir=None,
                  sstice_combined_file=None,
                  sstice_name=None,
                  topo_file=None,
                  sst_file=None,
                  ice_file=None,
                  lev_type=None,
                  check_input_files=None,
                  RRM_grid=None,
                  do_timers=None,
                  verbose=False,
                  verbose_indent='',
                ):
        
        if target_model is None: raise ValueError('target_model can not be None')

        self.target_model = target_model
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
        self.lev_type = lev_type
        self.sstice_name = sstice_name
        self.sst_file = sst_file
        self.ice_file = ice_file
        self.sstice_combined_file = sstice_combined_file
        self.sstice_nlat_src = None
        self.sstice_nlon_src = None
        self.sstice_nlat_dst = None
        self.sstice_nlon_dst = None

        self.RRM_grid  = RRM_grid  if RRM_grid  is not None else False
        self.do_timers = do_timers if do_timers is not None else True

        self.verbose = verbose
        self.verbose_indent = verbose_indent

        if check_input_files is None: check_input_files = True

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
        if check_input_files:
            for file_name in [self.atm_file,self.sfc_file,self.sst_file
                             ,self.ice_file,self.topo_file]:
                if file_name is not None:
                    if not os.path.exists(file_name):
                        raise ValueError(f'input file does not exist: {file_name}')

        # Load input files into xarray datasets
        if self.atm_file is not None: self.ds_atm = xr.open_dataset(self.atm_file)
        if self.sfc_file is not None: self.ds_sfc = xr.open_dataset(self.sfc_file)
    # --------------------------------------------------------------------------
    # Import SST/sea-ice methods
    from hiccup.hiccup_data_class_sstice_methods import get_sst_file
    from hiccup.hiccup_data_class_sstice_methods import sstice_create_src_grid_file
    from hiccup.hiccup_data_class_sstice_methods import open_combined_sstice_dataset
    from hiccup.hiccup_data_class_sstice_methods import sstice_create_dst_grid_file
    from hiccup.hiccup_data_class_sstice_methods import sstice_create_map_file
    from hiccup.hiccup_data_class_sstice_methods import sstice_slice_and_remap
    from hiccup.hiccup_data_class_sstice_methods import sstice_rename_vars
    from hiccup.hiccup_data_class_sstice_methods import sstice_adjustments
    # --------------------------------------------------------------------------
    # Import timer methods
    from hiccup.hiccup_data_class_timer_methods import timer_msg_all
    from hiccup.hiccup_data_class_timer_methods import timer_start_total
    from hiccup.hiccup_data_class_timer_methods import print_timer
    from hiccup.hiccup_data_class_timer_methods import print_timer_summary
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
        if 'np4' in self.dst_horz_grid:
            result = re.search('ne(.*)np', self.dst_horz_grid)
        if 'pg' in self.dst_horz_grid:
            result = re.search('ne(.*)pg', self.dst_horz_grid)
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
    def get_dst_grid_ncol(self):
        """
        Return ncol for destination grid
        """
        if self.RRM_grid:
            # use map file to determine ncol
            if self.map_file is None : raise ValueError('get_dst_grid_ncol: ncol cannot be determined')
            ds_grid = xr.open_dataset(self.map_file)
            ncol = int(ds_grid['n_b'].values)
        else:
            ne  = int(self.get_dst_grid_ne())
            npg = int(self.get_dst_grid_npg())
            if npg==0: ncol = int(ne*ne*6*9+2)
            if npg>0 : ncol = int(ne*ne*6*npg)
        return ncol
    # --------------------------------------------------------------------------
    def get_chunks(self,ncol_only=True):
        """
        Return chunk number to use dask for certain special cases
        Note: a chunk size of 1000 helps avoid memory problems for ne1024 cases.
        A smaller chunk size will reduce the memory footprint, but will also
        notably decrease the speed of calculations
        """
        if enable_chunks:
            if ncol_only:
                return {'ncol':ncol_chunk_size}
            else:
                return {'ncol':ncol_chunk_size,'lev':10}
        else:
            chunk = None
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
        if self.do_timers: timer_start = perf_counter()
        if verbose is None: verbose = self.verbose
        if verbose: print(f'\n{self.verbose_indent}Unpacking data files...')

        check_dependency('ncpdq')

        for f in [ self.atm_file, self.sfc_file, 
                   self.sst_file, self.ice_file,
                   self.sstice_combined_file ]:
            if f is not None :
                run_cmd(f'ncpdq -U --ovr {f} {f}',verbose,prepend_line=False)
        if self.do_timers: self.print_timer(timer_start)
        return
    # --------------------------------------------------------------------------
    def create_dst_grid_file(self,verbose=None):
        """ 
        Generate destination model grid file. Normally, we only care about 
        mapping to the GLL/np4 grid, unless the source data is an EAM file
        """
        if self.do_timers: timer_start = perf_counter()
        if verbose is None: verbose = self.verbose
        if verbose: print(f'\n{self.verbose_indent}Generating dst grid file...')

        if 'ne' in self.dst_horz_grid and 'np' in self.dst_horz_grid : 
            
            # Spectral element grid with physics on GLL nodes
            ne = self.get_dst_grid_ne()
            self.dst_grid_file = f'{self.grid_dir}/exodus_ne{ne}.g'
            
            check_dependency('GenerateCSMesh')
            cmd = f'GenerateCSMesh --alt --res {ne} --file {self.dst_grid_file}'
            cmd += f' >> {tempest_log_file}'
            run_cmd(cmd,verbose,shell=True)
        
        elif 'ne' in self.dst_horz_grid and 'pg' in self.dst_horz_grid :

            # Spectral element grid with physics on GLL nodes
            ne  = self.get_dst_grid_ne()
            npg = self.get_dst_grid_npg()

            # First create exodus file
            check_dependency('GenerateCSMesh')
            cmd = f'GenerateCSMesh --alt --res {ne} --file {self.dst_grid_file_np}'
            cmd += f' >> {tempest_log_file}'
            run_cmd(cmd,verbose,shell=True)
            
            # Next switch to volumetric mesh that matches the physgrid
            tmp_exodus_file = f'{self.grid_dir}/exodus_{self.dst_horz_grid_pg}.g'
            check_dependency('GenerateVolumetricMesh')
            cmd = 'GenerateVolumetricMesh'
            cmd += f' --in {self.dst_grid_file_np} '
            cmd += f' --out {tmp_exodus_file} '
            cmd += f' --np {npg} --uniform'
            cmd += f' >> {tempest_log_file}'
            run_cmd(cmd,verbose,shell=True)

            # Create scrip file while we're at it (can be slow)
            check_dependency('ConvertMeshToSCRIP')
            cmd = 'ConvertMeshToSCRIP'
            cmd += f' --in {tmp_exodus_file} '
            cmd += f' --out {self.dst_grid_file_pg} '
            cmd += f' >> {tempest_log_file}'
            run_cmd(cmd,verbose,shell=True)

            # # fix grid_imask type
            # run_cmd('ncap2 --overwrite -s \'grid_imask=int(grid_imask)\' '
            #         +f'{self.dst_grid_file_pg} {self.dst_grid_file_pg}',verbose,shell=True)

            # delete temporary exodus file
            run_cmd(f'rm {tmp_exodus_file} ',verbose)

        else:
            raise ValueError(f'dst_horz_grid={self.dst_horz_grid} is not currently supported')

        if self.do_timers: self.print_timer(timer_start)
        return 
    # --------------------------------------------------------------------------
    def create_map_file(self,verbose=None,src_type=None,dst_type=None,lrg2sml=False):
        """ 
        Generate mapping file after grid files have been created.
        This routine assumes that the destination is always GLL/np4.
        For mapping EAM to EAM data this method is overloaded below. 
        """
        if self.do_timers: timer_start = perf_counter()
        if verbose is None: verbose = self.verbose
        if verbose: print(f'\n{self.verbose_indent}Generating mapping file...')

        check_dependency('ncremap')

        # Check that grid file fields are not empty
        if self.src_grid_file is None : raise ValueError('src_grid_file is not defined!')
        if self.dst_grid_file is None : raise ValueError('dst_grid_file is not defined!')

        if src_type is None: src_type = 'FV' # assume input is FV
        if dst_type is None: dst_type = 'GLL' # assume dst grid is GLL/np4

        if src_type is not None and src_type not in ['FV','GLL']:
            raise ValueError(f'The value of src_type={src_type} is not supported')
        if dst_type is not None and dst_type not in ['FV','GLL']:
            raise ValueError(f'The value of src_type={src_type} is not supported')

        # Set the mapping algorithm
        if src_type=='FV' and dst_type=='GLL': alg_flag = '-a fv2se_flx'
        if src_type=='GLL'and dst_type=='GLL': alg_flag = '-a se2se'
        if src_type=='FV' and dst_type=='FV' : alg_flag = '-a fv2fv_flx'
        if src_type=='GLL'and dst_type=='FV' : alg_flag = '-a se2fv_flx'
        
        # Create the map file
        cmd = f'ncremap {alg_flag} '
        cmd += f' --src_grd={self.src_grid_file}'
        cmd += f' --dst_grd={self.dst_grid_file}'
        cmd += f' --map_file={self.map_file}'
        if lrg2sml: cmd += ' --lrg2sml ' # special flag for "very fine" grids 
        run_cmd(cmd,verbose,shell=True)

        if self.do_timers: self.print_timer(timer_start)
        return
    # --------------------------------------------------------------------------
    def get_multifile_dict(self,verbose=None,timestamp=None):
        """
        Create dict of temporary file names associated with each variable
        """
        if self.do_timers: timer_start = perf_counter()
        if verbose is None: verbose = self.verbose
        if verbose: print(f'\n{self.verbose_indent}Creating list of temporary files...')

        # define file list to be returned
        tmp_file_dict = {}

        if 'lat' in self.atm_var_name_dict: lat_var = self.atm_var_name_dict['lat']
        if 'lon' in self.atm_var_name_dict: lon_var = self.atm_var_name_dict['lon']

        var_dict_all = self.atm_var_name_dict.copy()
        var_dict_all.update(self.sfc_var_name_dict)

        # use timestamp to ensure these files are distinct from 
        # other instances of HICCUP that might be running concurrently
        if timestamp is None: timestamp = datetime.datetime.utcnow().strftime('%Y%m%d.%H%M%S')

        max_key_len = 0
        for key in var_dict_all.keys(): max_key_len = max( max_key_len, len(key) )

        # Horzontally remap atmosphere and surface data to individual files
        for key,var in var_dict_all.items():
            if var not in [lat_var,lon_var]:
                if key in self.sfc_var_name_dict.keys(): file_prefix = 'tmp_sfc_data'
                if key in self.atm_var_name_dict.keys(): file_prefix = 'tmp_atm_data'
                if file_prefix is not None:
                    tmp_file_name = f'{self.tmp_dir}/{file_prefix}'
                    tmp_file_name += f'.{self.dst_horz_grid}'
                    tmp_file_name += f'.{self.dst_vert_grid}'
                    tmp_file_name += f'.{key}'
                    tmp_file_name += f'.{timestamp}'
                    tmp_file_name += f'.nc'
                    tmp_file_dict.update({key:tmp_file_name})
                    if verbose: print(self.verbose_indent+f'  {key:{max_key_len}}   {tmp_file_name}')

        return tmp_file_dict
    # --------------------------------------------------------------------------
    def get_multifile_dict_eam(self,var_name_dict,verbose=None,timestamp=None,dst_horz_grid=None):
        """
        Create dict of temporary file names associated with each variable
        """
        if self.do_timers: timer_start = perf_counter()
        if verbose is None: verbose = self.verbose
        if verbose: print(f'\n{self.verbose_indent}Creating list of temporary files...')

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
            if verbose: print(self.verbose_indent+f'  {key:10}   {tmp_file_name}')

        return tmp_file_dict
    # --------------------------------------------------------------------------
    def rename_vars(self,file_name,verbose=None):
        """ 
        Rename variables in file according to variable name dictionaries 
        """
        if self.do_timers: timer_start = perf_counter()
        if verbose is None: verbose = self.verbose
        if verbose: print(f'\n{self.verbose_indent}Renaming variables to match model variable names...')

        check_dependency('ncrename')

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
        if self.do_timers: self.print_timer(timer_start)

        # Do additional variable/attribute renaming specific to the input data
        self.rename_vars_special(file_name,verbose)

        return
    # --------------------------------------------------------------------------
    def rename_vars_multifile(self,file_dict,verbose=None):
        """ 
        Rename variables in file list according to variable name dictionaries 
        This approach was developed specifically for very fine grids like ne1024
        """
        if self.do_timers: timer_start = perf_counter()
        if verbose is None: verbose = self.verbose
        if verbose: print(f'\n{self.verbose_indent}Renaming variables to match model variable names...')

        check_dependency('ncrename')

        if 'lat' in self.atm_var_name_dict: lat_var = self.atm_var_name_dict['lat']
        if 'lon' in self.atm_var_name_dict: lon_var = self.atm_var_name_dict['lon']

        var_dict_all = self.atm_var_name_dict.copy()
        var_dict_all.update(self.sfc_var_name_dict)

        new_lev_name = None
        if self.lev_name=='level': new_lev_name = 'plev'

        for var, file_name in file_dict.items():

            with xr.open_dataset(file_name,decode_cf=False) as ds_data:
                ds_data.load()
                # Rename the variable
                if var_dict_all[var] in ds_data: 
                    ds_data = ds_data.rename({var_dict_all[var]:var})
                # Rename the vertical coordinate
                if new_lev_name is not None \
                and '_sfc_' not in file_name \
                and self.lev_name in ds_data:
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

        if self.do_timers: self.print_timer(timer_start)

        return
    # --------------------------------------------------------------------------
    def add_reference_pressure(self,file_name,verbose=None):
        """ 
        Add P0 variable 
        """
        if self.do_timers: timer_start = perf_counter()
        if verbose is None: verbose = self.verbose
        if verbose: print(f'\n{self.verbose_indent}Adding reference pressure (P0)...')

        check_dependency('ncap2')
        check_dependency('ncatted')

        # Add the variable
        run_cmd(f"ncap2 --hst -A -s 'P0=100000.' {file_name} {file_name}",
                verbose,prepend_line=False,shell=True)

        # add long_name and units attributes
        run_cmd(f"ncatted --hst -a long_name,P0,a,c,'reference pressure' -a units,P0,a,c,'Pa' {file_name}",
                verbose,prepend_line=False,shell=True)
        
        if self.do_timers: self.print_timer(timer_start)
        return
    # --------------------------------------------------------------------------
    def remap_horizontal(self,output_file_name,verbose=None):
        """  
        Horizontally remap data and combine into single file 
        """
        if self.do_timers: timer_start = perf_counter()
        if verbose is None: verbose = self.verbose
        if verbose: print(f'\n{self.verbose_indent}Horizontally remapping the data to temporary files...')

        if self.map_file is None: raise ValueError('map_file cannot be None!')
        if self.atm_file is None: raise ValueError('atm_file cannot be None!')
        if self.sfc_file is None: raise ValueError('sfc_file cannot be None!')

        # Define temporary files that will be deleted at the end
        atm_tmp_file_name = f'{self.tmp_dir}/tmp_atm_data.nc'
        sfc_tmp_file_name = f'{self.tmp_dir}/tmp_sfc_data.nc'

        # Remove temporary files if they exist
        if os.path.isfile(atm_tmp_file_name): run_cmd(f'rm {atm_tmp_file_name} ',verbose)
        if os.path.isfile(sfc_tmp_file_name): run_cmd(f'rm {sfc_tmp_file_name} ',verbose)

        check_dependency('ncremap')
        check_dependency('ncks')

        # Horzontally remap atmosphere data
        var_list = ','.join(self.atm_var_name_dict.values())
        cmd =  f'ncremap'
        cmd += f' --map_file={self.map_file} '
        cmd += f' --in_file={self.atm_file} '
        cmd += f' --out_file={atm_tmp_file_name} '
        cmd += f' --var_lst={var_list} '
        cmd += f' --fl_fmt={ncremap_file_fmt} '
        run_cmd(cmd,verbose)

        # Horzontally remap surface data
        var_list = ','.join(self.sfc_var_name_dict.values())
        cmd =  f'ncremap'
        cmd += f' --map_file={self.map_file} '
        cmd += f' --in_file={self.sfc_file} '
        cmd += f' --out_file={sfc_tmp_file_name} '
        cmd += f' --var_lst={var_list} '
        cmd += f' --fl_fmt={ncremap_file_fmt} '
        run_cmd(cmd,verbose)

        # Remove output file if it already exists
        if os.path.isfile(output_file_name): run_cmd(f'rm {output_file_name} ',verbose)

        if verbose: print(f'\n{self.verbose_indent}Combining temporary remapped files...')

        # Add atmosphere temporary file data into the final output file
        run_cmd(f'ncks -A --hdr_pad={hdr_pad} {atm_tmp_file_name} {output_file_name} ',
                verbose,prepend_line=False)

        # Add surface temporary file data into the final output file
        run_cmd(f'ncks -A --hdr_pad={hdr_pad} {sfc_tmp_file_name} {output_file_name} ',
                verbose,prepend_line=False)

        # delete the temporary files
        run_cmd(f'rm {sfc_tmp_file_name} {atm_tmp_file_name} ',verbose,prepend_line=False)

        if self.do_timers: self.print_timer(timer_start)
        return
    # --------------------------------------------------------------------------
    def remap_horizontal_multifile(self,file_dict,verbose=None):
        """  
        Horizontally remap data into seperate files for each variable
        This approach was developed specifically for very fine grids like ne1024
        """
        if self.do_timers: timer_start = perf_counter()
        if verbose is None: verbose = self.verbose
        if verbose: print(f'\n{self.verbose_indent}Horizontally remapping the multi-file data to temporary files...')

        if self.map_file is None: raise ValueError('map_file cannot be None!')
        if self.atm_file is None: raise ValueError('atm_file cannot be None!')
        if self.sfc_file is None: raise ValueError('sfc_file cannot be None!')

        check_dependency('ncremap')

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
            cmd  = f'ncremap'
            cmd += f" --nco_opt='-O --no_tmp_fl --hdr_pad={hdr_pad}' "
            cmd += f' --map_file={self.map_file}'
            cmd += f' --in_file={in_file}'
            cmd += f' --out_file={tmp_file_name}'
            # cmd += f' --var_lst={in_var},{lat_var},{lon_var}'
            cmd += f' --var_lst={in_var_list}'
            cmd += f' --fl_fmt={ncremap_file_fmt}'
            run_cmd(cmd,verbose,shell=True)

        if self.do_timers: self.print_timer(timer_start)

        return 
    # --------------------------------------------------------------------------
    def remap_horizontal_multifile_eam(self,file_dict,verbose=None):
        """  
        Horizontally remap data into seperate files for each variable
        This approach was developed specifically for very fine grids like ne1024
        """
        if self.do_timers: timer_start = perf_counter()
        if verbose is None: verbose = self.verbose
        if verbose: print(f'\n{self.verbose_indent}Horizontally remapping the multi-file data to temporary files...')

        if self.map_file is None: raise ValueError('map_file cannot be None!')
        if self.atm_file is None: raise ValueError('atm_file cannot be None!')
        if self.sfc_file is None: raise ValueError('sfc_file cannot be None!')

        check_dependency('ncremap')

        # Horzontally remap atmosphere and surface data to individual files
        for var,tmp_file_name in file_dict.items():
            in_var = var
            in_file = self.atm_file
            # Remove temporary files if they exist
            if os.path.isfile(tmp_file_name): run_cmd(f'rm {tmp_file_name}',verbose)
            # Remap the data
            in_var_list = f'{in_var}'
            cmd  = f'ncremap'
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

        if self.do_timers: self.print_timer(timer_start)

        return 
    # --------------------------------------------------------------------------
    def surface_adjustment_multifile(self,file_dict,verbose=None,
                                    adj_TS=True,adj_PS=True):
        """
        Perform surface temperature and pressure adjustments 
        using a multifile xarray dataset
        """
        if self.do_timers: timer_start = perf_counter()
        if verbose is None: verbose = self.verbose
        if verbose: print(f'\n{self.verbose_indent}Performing surface adjustments...')

        # update lev name in case it has not been updated previously
        self.lev_name = self.new_lev_name

        if print_memory_usage: print_mem_usage(msg='start surface_adjustment_multifile()')

        adj_TS_warning_msg = 'WARNING - surface_adjustment_multifile: '+\
        f'\n  adj_TS in  is not supported for {self.target_model}, disabling.'+\
        f'\n  Set adj_TS=False  to suppress this warning.'

        # build dict of variables needed for adjustment
        var_dict = {}
        if self.target_model=='EAM':
            if adj_TS:
                var_dict['TS'] = 'TS'
            if adj_PS:
                var_dict['PS']   = 'PS'
                var_dict['PHIS'] = 'PHIS'
                var_dict['T']    = 'T'
        if self.target_model=='EAMXX':
            if adj_TS: adj_TS = False ; print(adj_TS_warning_msg)
            if adj_PS:
                var_dict['PS']   = 'ps'
                var_dict['PHIS'] = 'phis'
                var_dict['T']    = 'T_mid'

        file_list = []
        for var,file_name in file_dict.items():
            if var in var_dict.values():
                file_list.append(file_name)

        # Load topo data for surface adjustment - use same chunking
        ds_topo = xr.open_dataset(self.topo_file,chunks=self.get_chunks())

        # Adjust surface temperature to match new surface height
        if adj_TS:
            if self.do_timers: timer_start_adj = perf_counter()
            with xr.open_mfdataset(file_list,combine='by_coords',chunks=self.get_chunks()) as ds_data:
                ds_data = ds_data.rename(dict((val,key) for key,val in var_dict.items()))
                hsa.adjust_surface_temperature( ds_data, ds_topo, verbose=verbose, 
                                                verbose_indent=self.verbose_indent )
            ds_data = ds_data.rename(var_dict)
            ds_data[var_dict['TS']].to_netcdf(file_dict[var_dict['TS']],format=hiccup_atm_nc_format,mode='a')
            ds_data.close()
            if self.do_timers: self.print_timer(timer_start_adj,caller='adjust_surface_temperature')

        if print_memory_usage: print_mem_usage(msg='after adj_TS')

        # Adjust surface pressure to match new surface height
        if adj_PS:
            if self.do_timers: timer_start_adj = perf_counter()
            with xr.open_mfdataset(file_list,combine='by_coords',chunks=self.get_chunks()) as ds_data:
                ds_data = ds_data.rename(dict((val,key) for key,val in var_dict.items()))
                hsa.adjust_surface_pressure( ds_data, ds_topo, pressure_var_name=self.lev_name,
                                             lev_coord_name=self.lev_name, verbose=verbose, 
                                             verbose_indent=self.verbose_indent )
            ds_data = ds_data.rename(var_dict)
            ds_data[var_dict['PS']].to_netcdf(file_dict[var_dict['PS']],format=hiccup_atm_nc_format,mode='a')
            ds_data.close()
            if self.do_timers: self.print_timer(timer_start_adj,caller='adjust_surface_pressure')

        if print_memory_usage: print_mem_usage(msg='after adj_PS')

        if self.do_timers: self.print_timer(timer_start)

        return

    # --------------------------------------------------------------------------
    def remap_vertical(self,input_file_name,output_file_name,
                       vert_file_name,ps_name='PS',vert_remap_var_list=None,
                       verbose=None):
        """  
        Vertically remap data and combine into single file 
        """
        if self.do_timers: timer_start = perf_counter()
        if verbose is None: verbose = self.verbose
        if verbose: print(f'\n{self.verbose_indent}Vertically remapping the data...')

        check_dependency('ncremap')

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
        cmd  = 'ncremap'
        cmd += f" --nco_opt='-O --no_tmp_fl --hdr_pad={hdr_pad}' " # doesn't work with vertical regridding?
        cmd += f' --vrt_fl={vert_file_name}'
        cmd += f' --ps_nm={ps_name}'
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

        if self.do_timers: self.print_timer(timer_start)
        return
    # --------------------------------------------------------------------------
    def remap_vertical_multifile(self,file_dict,vert_file_name,verbose=None):
        """
        wrapper around remap_vertical to support multi-file workflow
        specifically needed for very fine grids like ne1024
        """
        if verbose is None: verbose = self.verbose
        if verbose: print(f'\n{self.verbose_indent}Vertically remapping the multi-file data...')

        # temporarily disable timers and put a timer around the vertical remap loop
        global do_timers
        if self.do_timers: timer_start = perf_counter()
        prev_do_timers = self.do_timers
        self.do_timers = False

        ps_var_name = None
        if self.target_model=='EAM'  : ps_var_name = 'PS'
        if self.target_model=='EAMXX': ps_var_name = 'ps'

        if ps_var_name is None:
            raise ValueError(f'Error: ps_var_name is not specified for target_model: {self.target_model}')

        ps_file_name = file_dict[ps_var_name]

        for var,file_name in file_dict.items() :
            if '_sfc_' not in file_name :
                # Append surface pressure for vertical interpolation
                ds = xr.open_dataset(file_name)
                ds_ps = xr.open_dataset(ps_file_name)
                ds[ps_var_name] = ds_ps[ps_var_name]
                tmp_file_name = file_name.replace('.nc',f'.tmp.nc')
                ds.to_netcdf(tmp_file_name,mode='w')
                ds.close(); ds_ps.close()
                run_cmd(f'mv {tmp_file_name} {file_name}',verbose=False)

                # Do the vertical interpolation for this file
                self.remap_vertical(input_file_name=file_name,
                                   output_file_name=file_name,
                                   vert_file_name=vert_file_name,
                                   vert_remap_var_list=[var],
                                   ps_name=ps_var_name)
        
        # Re-set do_timers to previous value
        self.do_timers = prev_do_timers

        if self.do_timers: self.print_timer(timer_start)
        return
    # --------------------------------------------------------------------------
    def atmos_state_adjustment_multifile(self,file_dict,verbose=None,
                                        adjust_sat=True,adjust_wtr=True,
                                        convert_ozone=True):
        """
        Perform post-remapping atmospheric state adjustments 
        for the multifile workflow
        """
        if self.do_timers: timer_start = perf_counter()
        if verbose is None: verbose = self.verbose
        if verbose: print(f'\n{self.verbose_indent}Performing state adjustments...')

        # get list of file names for variables needed for adjustment
        def get_adj_file_list(var_list):
            file_list = []
            for var,file_name in file_dict.items():
                if var in var_list: file_list.append(file_name)
            return file_list

        if print_memory_usage: print_mem_usage(msg='start atmos_state_adjustment_multifile()')

        if adjust_sat:
            if self.target_model=='EAM'  : var_dict = {'Q':'Q', 'T':'T',    'PS':'PS'}
            if self.target_model=='EAMXX': var_dict = {'Q':'qv','T':'T_mid','PS':'ps'}
            file_list = get_adj_file_list(var_dict.values())
            if self.do_timers: timer_start_adj = perf_counter()

            with xr.open_mfdataset(file_list,combine='by_coords',chunks=self.get_chunks()) as ds_data:
                ds_data = ds_data.rename(dict((val,key) for key,val in var_dict.items()))
                # adjust water vapor to eliminate supersaturation
                if print_memory_usage: print_mem_usage(msg='before remove_supersaturation')
                hsa.remove_supersaturation( ds_data, hybrid_lev=True, verbose=verbose, verbose_indent=self.verbose_indent )
                if print_memory_usage: print_mem_usage(msg='after remove_supersaturation')
                # Write adjusted data back to data files
                ds_data = ds_data.rename(var_dict)
                ds_data[var_dict['Q']].to_netcdf(file_dict[var_dict['Q']],format=hiccup_atm_nc_format,mode='a')
                ds_data.close()
            if self.do_timers: self.print_timer(timer_start_adj,caller='remove_supersaturation')

        if adjust_wtr:
            if self.target_model=='EAM'  : var_dict = {'CLDLIQ':'CLDLIQ','CLDICE':'CLDICE'}
            if self.target_model=='EAMXX': var_dict = {'CLDLIQ':'qc',    'CLDICE':'qi'}
            file_list = get_adj_file_list(var_dict.values())
            if self.do_timers: timer_start_adj = perf_counter()
            with xr.open_mfdataset(file_list,combine='by_coords',chunks=self.get_chunks()) as ds_data:
                ds_data = ds_data.rename(dict((val,key) for key,val in var_dict.items()))
                # adjust cloud water to remove negative values
                if print_memory_usage: print_mem_usage(msg='before adjust_cld_wtr')
                hsa.adjust_cld_wtr( ds_data, verbose=verbose, verbose_indent=self.verbose_indent )
                if print_memory_usage: print_mem_usage(msg='after adjust_cld_wtr')
                # Write adjusted data back to data files
                ds_data = ds_data.rename(var_dict)
                for var in var_dict.values():
                    if var in self.atm_var_name_dict.keys():
                        ds_data[var].to_netcdf(file_dict[var],format=hiccup_atm_nc_format,mode='a')
                ds_data.close()
            if self.do_timers: self.print_timer(timer_start_adj,caller='adjust_cld_wtr')

        if convert_ozone:
            if verbose: print(f'\n{self.verbose_indent}Converting Ozone to molecular/volume mixing ratio...')
            if self.do_timers: timer_start_adj = perf_counter()
            if self.target_model=='EAM'  : O3_name = 'O3'
            if self.target_model=='EAMXX': O3_name = 'o3_volume_mix_ratio'
            ds_data = xr.open_mfdataset(file_dict[O3_name],combine='by_coords',chunks=self.get_chunks())
            # Convert mass mixing ratio to molecular/volume mixing ratio
            ds_data[O3_name] = ds_data[O3_name] * MW_dryair / MW_ozone
            ds_data[O3_name].attrs['units'] = 'mol/mol'
            ds_data.to_netcdf(file_dict[O3_name],format=hiccup_atm_nc_format,mode='a')
            ds_data.close()
            if self.do_timers: self.print_timer(timer_start_adj,caller='convert_ozone')

        if self.do_timers: self.print_timer(timer_start)
        return
    # --------------------------------------------------------------------------
    def atmos_state_apply_perturbations_multifile(self,file_dict,seed=None,verbose=None):
        """
        apply post-remapping atmospheric perturbations
        for the multifile workflow
        """
        if self.do_timers: timer_start = perf_counter()
        if verbose is None: verbose = self.verbose
        if verbose: print(f'\n{self.verbose_indent}Applying random perturbations...')

        # build list of file names for variables to be perturbed
        file_list = []
        var_list = ['T','PS','U','V']
        for var,file_name in file_dict.items():
            if var in var_list: file_list.append(file_name)

        with xr.open_mfdataset(file_list,combine='by_coords',chunks=self.get_chunks()) as ds_data:

            # adjust cloud water to remove negative values
            hsa.apply_random_perturbations( ds_data, var_list=var_list, seed=seed, 
                                            verbose=False, verbose_indent=self.verbose_indent )
            ds_data.compute()

        # Write perturbed data back to the individual data files
        for var in var_list:
            if var in self.atm_var_name_dict.keys():
                ds_data[var].to_netcdf(file_dict[var],format=hiccup_atm_nc_format,mode='a')

        ds_data.close()

        if self.do_timers: self.print_timer(timer_start)
        return
    # --------------------------------------------------------------------------
    def add_time_date_variables(self,ds,verbose=None,do_timers=None):
        """
        Check final output file and add necessary time and date information
        """
        if do_timers is None: do_timers = self.do_timers
        if do_timers: timer_start = perf_counter()
        if verbose is None: verbose = self.verbose
        if verbose: print(f'\n{self.verbose_indent}Editing time and date variables...')

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
        for t in range(len(ds['time'])) : time_bnds[t,:] = ds['time'][t].values
        ds['time_bnds'] = xr.DataArray( time_bnds, coords=time_coord, dims=['time','nbnd'] )
        ds['time_bnds'].attrs['long_name'] = 'time interval endpoints'

        # string representation of date
        # date_str = str(year).zfill(4)+str(month).zfill(2)+str(day).zfill(2)

        # Integer representation of date
        date_int = year*1e4 + month*1e2 + day

        # miscellaneous time/date variables
        if 'date' not in ds.variables :
            ds['date'] = xr.DataArray( np.array(date_int,dtype=int),
                                        coords=time_coord, dims=time_dim )
            ds['date'].attrs['long_name'] = 'current date (YYYYMMDD)'
        if 'datesec' not in ds.variables :
            ds['datesec'] = xr.DataArray( np.array(sec,dtype=int),
                                        coords=time_coord, dims=time_dim )
            ds['datesec'].attrs['long_name'] = 'current seconds of current date'
        if 'ndcur' not in ds.variables :
            ds['ndcur'] = xr.DataArray( np.full(time_shape,0.,dtype=int),
                                        coords=time_coord, dims=time_dim )
            ds['ndcur'].attrs['long_name'] = 'current day (from base day)'
        if 'nscur' not in ds.variables :
            ds['nscur'] = xr.DataArray( np.array(sec,dtype=int),
                                        coords=time_coord, dims=time_dim )
            ds['nscur'].attrs['long_name'] = 'current seconds of current day'
        if 'nsteph' not in ds.variables :
            ds['nsteph'] = xr.DataArray( np.full(time_shape,0.,dtype=int),
                                        coords=time_coord, dims=time_dim )
            ds['nsteph'].attrs['long_name'] = 'current timestep'

        # Base day time - not sure what this is for
        if 'nbdate' not in ds.variables :
            ds['nbdate'] = xr.DataArray( np.array(date_int,dtype=int),
                                        coords=time_coord, dims=time_dim )
            ds['nbdate'].attrs['long_name'] = 'base date (YYYYMMDD)'
        if 'ndbase' not in ds.variables :
            ds['ndbase'] = xr.DataArray( np.array(day,dtype=int),
                                        coords=time_coord, dims=time_dim )
            ds['ndbase'].attrs['long_name'] = 'base day'
        if 'nsbase' not in ds.variables :
            ds['nsbase'] = xr.DataArray( np.array(sec,dtype=int),
                                        coords=time_coord, dims=time_dim )
            ds['nsbase'].attrs['long_name'] = 'seconds of base day'
        if 'nbsec' not in ds.variables :
            ds['nbsec'] = xr.DataArray( np.array(sec,dtype=int),
                                        coords=time_coord, dims=time_dim )
            ds['nbsec'].attrs['long_name'] = 'seconds of base date'

        # (Not sure how to handle string/char variables, but are they needed?)

        # if 'date_written' not in ds.variables :
        #     ds['date_written'] = xr.DataArray( np.full(time_shape,0.), coords=time_coord, dims=time_dim )
        #     ds['date_written'].attrs['long_name'] = ''

        # if 'time_written' not in ds.variables :
        #     ds['time_written'] = xr.DataArray( np.full(time_shape,0.), coords=time_coord, dims=time_dim )
        #     ds['time_written'].attrs['long_name'] = ''

        if do_timers: self.print_timer(timer_start)
        return
    # --------------------------------------------------------------------------
    def add_time_date_variables_multifile(self,file_dict,verbose=None,ref_date='1850-01-01'):
        """
        """
        if self.do_timers: timer_start = perf_counter()
        if verbose is None: verbose = self.verbose
        if verbose: print(f'\n{self.verbose_indent}Editing time and date variables...')

        # xarray will automatically convert the time coordinate,
        # which can be problematic, especially when generating nudging data.
        # By specifying the encoding we can avoid this problem
        time_encoding_dict = {'time':{'units': f'hours since {ref_date} 00:00:00'}}

        for file_name in file_dict.values() :
            # run_cmd(f'ncdump {file_name} -v time | tail | grep "time =" ',verbose,shell=True)
            # run_cmd(f'ncdump {file_name} -h | grep "time:units" ',verbose,shell=True)
            with xr.open_dataset(file_name) as ds_data:
                ds_data.load()
                self.add_time_date_variables(ds_data,verbose=False,do_timers=False)
            ds_data.to_netcdf(file_name,format=hiccup_atm_nc_format,mode='w',encoding=time_encoding_dict)
            ds_data.close()
            # run_cmd(f'ncdump {file_name} -v time | tail | grep "time =" ',verbose,shell=True)
            # run_cmd(f'ncdump {file_name} -h | grep "time:units" ',verbose,shell=True)
        if self.do_timers: self.print_timer(timer_start)
        return
    # --------------------------------------------------------------------------
    def add_time_date_variables_multifile_eam(self,file_dict,verbose=None):
        """
        copy necessary time and date information
        """
        if self.do_timers: timer_start = perf_counter()
        if verbose is None: verbose = self.verbose
        if verbose: print(f'\n{self.verbose_indent}Editing time and date variables...')

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
        if self.do_timers: self.print_timer(timer_start)
        return
    # --------------------------------------------------------------------------
    def convert_to_single_precision_multifile(self,file_dict,verbose=None):
        """
        Convert all files in file_dict to single precision
        """
        if self.do_timers: timer_start = perf_counter()
        if verbose is None: verbose = self.verbose
        if verbose: print(f'\n{self.verbose_indent}Converting all data files to single precision...')

        check_dependency('ncpdq')

        for var,file_name in file_dict.items() :
            tmp_file_name = file_name.replace('.nc',f'.tmp.nc')
            cmd = f'ncpdq -O --pck_map dbl_flt {file_name} {tmp_file_name}'
            run_cmd(cmd,verbose)
            run_cmd(f'mv {tmp_file_name} {file_name}',verbose)

        if self.do_timers: self.print_timer(timer_start)
        return
    # --------------------------------------------------------------------------
    def combine_files(self,file_dict,output_file_name,delete_files=False,
                      method='xarray',use_single_precision=False,
                      permute_dimensions=True,permute_dim_list=['time','ncol','lev'],
                      combine_uv=False,verbose=None):
        """
        Combine files in file_dict into single output file
        """
        if self.do_timers: timer_start = perf_counter()
        if verbose is None: verbose = self.verbose
        if verbose: print(f'\n{self.verbose_indent}Combining temporary files into new file...')
        if method not in ['xarray']:
            raise ValueError(f'combine_files method argument "{method}" is invalid')

        if os.path.isfile(output_file_name): 
            run_cmd(f'rm {output_file_name} ',verbose)

        if self.target_model=='EAM'  : 
            u_name,v_name,uv_name = 'U','V','UV'
        if self.target_model=='EAMXX': 
            u_name,v_name,uv_name = 'horiz_winds_u','horiz_winds_v','horiz_winds'

        # Combine temporary files into the final output file
        if method=='xarray':
            ds_out = xr.Dataset()
            for var,file_name in file_dict.items() :
                ds_tmp = xr.open_dataset(file_name)
                ds_out[var] = ds_tmp[var]
                if use_single_precision:
                    ds_out[var] = ds_tmp[var].astype('float32')
                else:
                    ds_out[var] = ds_tmp[var].astype('float64')
                ds_tmp.close()
            if permute_dimensions: 
                ds_out.transpose(permute_dim_list[0],
                                 permute_dim_list[1],
                                 permute_dim_list[2])
            if combine_uv:
                ds_out[uv_name] = xr.concat([ds_out[u_name], ds_out[v_name]], dim='dim2')
                ds_out[uv_name] = ds_out[uv_name].transpose('time','ncol','dim2','lev')
                ds_out = ds_out.drop_vars([u_name,v_name])
            ds_out.to_netcdf(output_file_name)
            ds_out.close()

        # # NCO "append" method
        # # this method tends to be slower than xarray, but also I couldn't figure
        # # out how to support the combine_uv and use_single_precision options
        # # so I've disabled this option for now
        # if method=='nco':
        #     check_dependency('ncks')
        #     if use_single_precision:
        #         self.convert_to_single_precision_multifile(file_dict,verbose)
        #     # if combine_uv:
        #     #     run_cmd(f'ncap2 --overwrite -s \'{uv_name}={u_name}\' ')
        #     for var,file_name in file_dict.items() :
        #         cmd = f'ncks -A --hdr_pad={hdr_pad} --no_tmp_fl'
        #         cmd+= f' --fl_fmt={ncremap_file_fmt}'
        #         cmd+= f' {file_name} {output_file_name} '
        #         run_cmd(cmd,verbose,prepend_line=False)
        #     if permute_dimensions:
        #         dim_str = ','.join(permute_dim_list)
        #         tmp_output_file_name = output_file_name.replace('.nc',f'.tmp.nc')
        #         cmd = f'ncpdq -O --permute {dim_str}'
        #         cmd+= f' {output_file_name} {tmp_output_file_name} '
        #         run_cmd(cmd,verbose)
        #         run_cmd(f'mv {tmp_output_file_name} {output_file_name}',verbose)
            

        # make sure time dimension is "unlimited"
        cmd = f'ncks -O --fl_fmt={ncremap_file_fmt}'
        cmd+= f' --mk_rec_dmn time'
        cmd+= f' {output_file_name} {output_file_name} '
        run_cmd(cmd,verbose,prepend_line=False)

        # Delete temp files
        if delete_files:
            if verbose: print(f'\n{self.verbose_indent}Deleting temporary files...')
            for file_name in file_dict.values() :
                run_cmd(self.verbose_indent+f'rm {file_name}',verbose,prepend_line=False)

        if self.do_timers: self.print_timer(timer_start)
        return
    # --------------------------------------------------------------------------
    def clean_global_attributes(self,file_name,method='nco',verbose=None):
        """
        Remove messy global attributes of the file
        """
        if self.do_timers: timer_start = perf_counter()
        if verbose is None: verbose = self.verbose
        if verbose: print(f'\n{self.verbose_indent}Cleaning up excessive global attributes...')

        global_att_list = ['history_of_appended_files', 'nco_openmp_thread_number', 
                           'input_file', 'map_file', 'remap_version', 'remap_hostname', 
                           'remap_command', 'remap_script', 'NCO' ]

        # Remove the attributes listed in global_att_list using xarray
        if method=='xarray':
            ds = xr.open_dataset(file_name)
            for att in global_att_list:
                if att in ds.attrs: del ds.attrs[att]
            ds.attrs['history'] = ''
            ds.to_netcdf(file_name)
            ds.close()
            # reset file format since xarray will automatically convert to netcdf4
            cmd = f'ncks -O --fl_fmt={ncremap_file_fmt} {file_name} {file_name} '
            run_cmd(cmd,verbose,prepend_line=False)

        # Remove the attributes listed in global_att_list using ncatted
        if method=='nco':
            check_dependency('ncatted')
            cmd = 'ncatted -O '
            for att in global_att_list: cmd += f' -a {att},global,d,, '
            cmd += f' {file_name} {file_name} '
            run_cmd(cmd,verbose)

        # Also reset the history attribute
        run_cmd(f'ncatted -h -O -a history,global,o,c, {file_name} {file_name}',
                verbose,prepend_line=False)

        if self.do_timers: self.print_timer(timer_start,caller=f'clean_global_attributes_{method}')
        return
# ------------------------------------------------------------------------------
# HICCUP Subclasses (associated with the source data)
# ------------------------------------------------------------------------------
import hiccup.hiccup_data_subclass_ERA5
import hiccup.hiccup_data_subclass_CAMS
import hiccup.hiccup_data_subclass_NOAA
import hiccup.hiccup_data_subclass_EAM
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
