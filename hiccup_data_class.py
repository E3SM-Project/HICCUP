# Class structure for defining parameters specific to the 
# input data used for creating initial condition files
# 
# NOTE: Variable name dictionaries are defined with the key as the model's 
# variable name and the value as the reanalysis data variable name

import numpy as np
import xarray as xr
import subprocess as sp
import shutil 
import re
import glob

# use current path as default output location
default_output_dir  = './'

# algorithm flag for ncremap
ncremap_alg         = ' --alg_typ=tempest '        

# log file for Tempest output
tempest_log_file    = 'TempestRemap.log'

hiccup_verbose = False

# Set numpy to ignore overflow errors
np.seterr(over='ignore')

# Set up terminal colors
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

# ------------------------------------------------------------------------------
# Common method for printing and running commands
# ------------------------------------------------------------------------------
def run_cmd(cmd,verbose=None,prefix='\n  ',suffix='',use_color=True,shell=False):
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
    """ Check for required system commands"""
    if shutil.which(cmd) is None : raise OSError(f'{cmd} is not in system path')
    return
# ------------------------------------------------------------------------------
# Method for returning class object
# ------------------------------------------------------------------------------
def create_hiccup_data(name,atm_file,sfc_file,dst_horz_grid,dst_vert_grid,
                       output_dir=default_output_dir,lev_type='',verbose=False):
    """ Return HICCUP data class object """
    global hiccup_verbose
    hiccup_verbose = verbose
    for subclass in hiccup_data.__subclasses__():
        if subclass.is_name_for(name):
            # Create the object
            obj = subclass(name
                      ,atm_file=atm_file
                      ,sfc_file=sfc_file
                      ,dst_horz_grid=dst_horz_grid
                      ,dst_vert_grid=dst_vert_grid
                      ,output_dir=output_dir
                      ,lev_type=lev_type)
            # Check input files for for required variables
            obj.check_file_vars()
            # Return the object if everything checks out
            return obj
    raise ValueError(f'{name} is not a valid HICCUP dataset name')
# ------------------------------------------------------------------------------
# Base Class
# ------------------------------------------------------------------------------
class hiccup_data(object):
    """ Base class for HICCUP data object """
    def __init__(self,atm_file,sfc_file,dst_horz_grid,dst_vert_grid,
                 output_dir=default_output_dir,lev_type=''):
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

        # Make sure output directory is formatted correctly
        if output_dir=='' : output_dir = './'
        if not output_dir.endswith('/'): output_dir += '/'
        self.output_dir = output_dir

        # Load input files into xarray datasets
        self.ds_atm = xr.open_dataset(self.atm_file)
        self.ds_sfc = xr.open_dataset(self.sfc_file)
    # --------------------------------------------------------------------------
    def __str__(self):
        str_out = ''
        for key in self.__dict__.keys(): 
            attribute = getattr(self,key)
            if isinstance(attribute,dict):
                str_out = str_out+f'  {key:15}:\n'
                for k in attribute.keys(): 
                    str_out = str_out+f'      {k:8}  {attribute[k]}\n'
            else:
                if attribute!='' :
                    str_out = str_out+f'  {key:15}:  {attribute}\n'
        return str_out
    # --------------------------------------------------------------------------
    def check_file_vars(self):
        """ Check that required variables are in the input files """

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
    def create_dst_grid_file(self,verbose=None):
        """ Generate destination model grid file """
        if verbose is None : verbose = hiccup_verbose
        if verbose : print('\nGenerating dst grid file...')

        if 'ne' in self.dst_horz_grid and 'np' in self.dst_horz_grid : 
            
            # Spectral element grid with physics on GLL nodes
            ne = re.search('ne(.*)np', self.dst_horz_grid).group(1)
            self.dst_grid_file = self.output_dir+f'exodus_ne{ne}.g'
            
            check_dependency('GenerateCSMesh')
            cmd = f'GenerateCSMesh --res {ne} --file {self.dst_grid_file}'
            cmd += f' >> {tempest_log_file}'
            run_cmd(cmd,verbose,shell=True)

        elif 'ne' in self.dst_horz_grid and 'pg' in self.dst_horz_grid : 
            
            # Spectral element grid with FV physics grid (ex. ne30pg2)
            ne  = re.search('ne(.*)pg', self.dst_horz_grid).group(1)
            npg = re.search('pg(.*)', self.dst_horz_grid).group(1)
            exodus_file = self.output_dir+f'exodus_ne{ne}.g'

            # First create exodus file
            check_dependency('GenerateCSMesh')
            cmd = f'GenerateCSMesh --res {ne} --file {exodus_file}'
            cmd += f' >> {tempest_log_file}'
            run_cmd(cmd,verbose,shell=True)
            
            # Next switch to volumetric mesh that matches the physgrid
            self.dst_grid_file = self.output_dir+f'exodus_{self.dst_horz_grid}.nc'
            check_dependency('GenerateVolumetricMesh')
            cmd = 'GenerateVolumetricMesh'
            cmd += f' --in {exodus_file} '
            cmd += f' --out {self.dst_grid_file} '
            cmd += f' --np {npg} --uniform'
            cmd += f' >> {tempest_log_file}'
            run_cmd(cmd,verbose,shell=True)

            # # Create scrip file while we're at it (can be slow)
            check_dependency('ConvertExodusToSCRIP')
            scrip_file = self.output_dir+f'scrip_{self.dst_horz_grid}.nc'
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
        """ Generate mapping file after grid files have been created """
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
        else:
            self.map_opts = self.map_opts+' --out_type fv --out_np 1 --volumetric '
        
        cmd = f'ncremap {ncremap_alg} '
        cmd += f' --src_grd={self.src_grid_file}'
        cmd += f' --dst_grd={self.dst_grid_file}'
        cmd += f' --map_file={self.map_file}'
        cmd += f' --wgt_opt=\'{self.map_opts}\' '
        run_cmd(cmd,verbose,shell=True)

        return
    # --------------------------------------------------------------------------
    def rename_vars(self,file_name,verbose=None):
        """ rename variables in file according to variable name dictionaries """
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
            run_cmd(cmd,verbose,prefix='  ',suffix='')

        if verbose : print('\n Renaming ATM vars... \n')
        for key in self.atm_var_name_dict : rename_proc(key,self.atm_var_name_dict)

        if verbose : print('\n Renaming SFC vars... \n')
        for key in self.sfc_var_name_dict : rename_proc(key,self.sfc_var_name_dict)

        self.rename_vars_special(file_name,verbose)

        return
    # --------------------------------------------------------------------------
    def add_reference_pressure(self,file_name,verbose=None):
        """ add P0 variable """
        if verbose is None : verbose = hiccup_verbose
        if verbose : print('\nAdding reference pressure (P0)...')

        check_dependency('ncap2')
        check_dependency('ncatted')

        # Add the variable
        run_cmd(f"ncap2 --create_ram --no_tmp_fl --hst -A -s 'P0=100000.' {file_name} {file_name}",
                verbose,prefix='  ',suffix='',shell=True)

        # add long_name attribute
        run_cmd(f"ncatted --hst -A -a long_name,P0,a,c,'reference pressure' {file_name}",
                verbose,prefix='  ',suffix='',shell=True)
        
        # add units attribute
        run_cmd(f"ncatted --hst -A -a units,P0,a,c,'Pa' {file_name}",
                verbose,prefix='  ',suffix='',shell=True)

        return
    # --------------------------------------------------------------------------
    def remap_horizontal(self,output_file_name,verbose=None):
        """  Horizontally remap data and combine into single file """
        if verbose is None : verbose = hiccup_verbose
        if verbose : print('\nHorizontally remapping the data to temporary files...')

        if self.map_file is None: raise ValueError('map_file cannot be None!')
        if self.atm_file is None: raise ValueError('atm_file cannot be None!')
        if self.sfc_file is None: raise ValueError('sfc_file cannot be None!')

        # Define temporary files that will be deleted at the end
        atm_tmp_file_name = './tmp_atm_data.nc'
        sfc_tmp_file_name = './tmp_sfc_data.nc'

        check_dependency('ncremap')
        check_dependency('ncks')

        # Horzontally remap atmosphere data
        var_list = ','.join(self.atm_var_name_dict.values())
        cmd =  f'ncremap {ncremap_alg} '
        cmd += f' --map_file={self.map_file} '
        cmd += f' --in_file={self.atm_file} '
        cmd += f' --out_file={atm_tmp_file_name} '
        cmd += f' --var_lst={var_list} '
        run_cmd(cmd,verbose)

        # Horzontally remap surface data
        var_list = ','.join(self.sfc_var_name_dict.values())
        cmd =  f'ncremap {ncremap_alg} '
        cmd += f' --map_file={self.map_file} '
        cmd += f' --in_file={self.sfc_file} '
        cmd += f' --out_file={sfc_tmp_file_name} '
        cmd += f' --var_lst={var_list} '
        run_cmd(cmd,verbose)

        # Remove output file if it already exists
        if output_file_name in glob.glob('*') : run_cmd(f'rm {output_file_name} ',verbose)

        if verbose : print('\nCombining temporary remapped files...')

        # Add atmosphere temporary file data into the final output file
        run_cmd(f'ncks -A {atm_tmp_file_name} {output_file_name} ',
                verbose,prefix='  ',suffix='')

        # Add surface temporary file data into the final output file
        run_cmd(f'ncks -A {sfc_tmp_file_name} {output_file_name} ',
                verbose,prefix='  ',suffix='')

        # delete the temporary files
        run_cmd(f'rm {sfc_tmp_file_name} {atm_tmp_file_name} ',
                verbose,prefix='  ',suffix='')

        return
    # --------------------------------------------------------------------------
    def remap_vertical(self,input_file_name,output_file_name,
                       vert_file_name,vert_remap_var_list=None,
                       verbose=None):
        """  Vertically remap data and combine into single file """
        if verbose is None : verbose = hiccup_verbose
        if verbose : print('\nVertically remapping the data...')

        check_dependency('ncremap')

        # Build variable list from the input file if not supplied
        if vert_remap_var_list is None :
            vert_remap_var_list = []
            ds = xr.open_dataset(input_file_name)
            for key in ds.variables.keys(): 
                vert_remap_var_list.append(key)
                # only remap varaibles with lev coord - this ignores other variables like TS!
                # if self.lev_name in ds[key].dims and key!=self.lev_name :
                #     vert_remap_var_list.append(key)
        vert_remap_var_list = ','.join(vert_remap_var_list)

        # Perform the vertical remapping
        cmd = 'ncremap'
        cmd+= f' --vrt_fl={vert_file_name}'
        cmd+= f' --var_lst={vert_remap_var_list}'
        cmd+= f' --in_fl={input_file_name}'
        cmd+= f' --out_fl={output_file_name}'
        run_cmd(cmd,verbose,shell=True)

        return
# ------------------------------------------------------------------------------
# Subclasses
# ------------------------------------------------------------------------------
class ERA5(hiccup_data):
    @classmethod
    def is_name_for(cls,name) : return name == 'ERA5'
    def __init__(self,name,atm_file,sfc_file,dst_horz_grid,dst_vert_grid,
                 output_dir=default_output_dir,lev_type=''):
        super().__init__(atm_file=atm_file
                        ,sfc_file=sfc_file
                        ,dst_horz_grid=dst_horz_grid
                        ,dst_vert_grid=dst_vert_grid
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
        # self.sfc_var_name_dict.update({'TS1':'stl1'})      # Soil temperature level 1 
        # self.sfc_var_name_dict.update({'TS2':'stl2'})      # Soil temperature level 2 
        # self.sfc_var_name_dict.update({'TS3':'stl3'})      # Soil temperature level 3 
        # self.sfc_var_name_dict.update({'TS4':'stl4'})      # Soil temperature level 4 
        # self.sfc_var_name_dict.update({'':'lai_hv'})       # Leaf area index, high vegetation 
        # self.sfc_var_name_dict.update({'':'lai_lv'})       # Leaf area index, low vegetation 
        # self.sfc_var_name_dict.update({'':'src'})          # Skin reservoir content 
        # self.sfc_var_name_dict.update({'':'asn'})          # Snow albedo 
        # self.sfc_var_name_dict.update({'':'rsn'})          # Snow density 
        # self.sfc_var_name_dict.update({'':'sd'})           # Snow depth 
        # self.sfc_var_name_dict.update({'':'tsn'})          # Temperature of snow layer 
        # self.sfc_var_name_dict.update({'':'swvl1'})        # Volumetric soil water level 1 
        # self.sfc_var_name_dict.update({'':'swvl2'})        # Volumetric soil water level 2 
        # self.sfc_var_name_dict.update({'':'swvl3'})        # Volumetric soil water level 3 
        # self.sfc_var_name_dict.update({'':'swvl4'})        # Volumetric soil water level 4 

        self.nlat = len( self.ds_atm[ self.atm_var_name_dict['lat'] ].values )
        self.nlon = len( self.ds_atm[ self.atm_var_name_dict['lon'] ].values )

        self.src_grid_name = f'{self.nlat}x{self.nlon}'
        self.src_grid_file = self.output_dir+f'scrip_{self.name}_{self.src_grid_name}.nc'

        self.map_file = self.output_dir+f'map_{self.src_grid_name}_to_{self.dst_horz_grid}.nc'

    # --------------------------------------------------------------------------
    def create_src_grid_file(self,verbose=None):
        """ Generate source grid file """
        if verbose is None : verbose = hiccup_verbose
        if verbose : print('\nGenerating src grid file...')

        # Remove the file here to prevent the warning message when ncremap overwrites it
        if self.src_grid_file in glob.glob('*') : 
            run_cmd(f'rm {self.src_grid_file} ',verbose)

        check_dependency('ncremap')

        cmd  = f'ncremap {ncremap_alg} ' 
        cmd += f' --tmp_dir=./tmp'
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
        """ rename file vars specific to this subclass """
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
                verbose,prefix='  ',suffix='',shell=True)

        # change units attribute
        run_cmd(f"ncatted --hst -A -a units,{new_lev_name},a,c,'Pa' {file_name}",
                verbose,prefix='  ',suffix='',shell=True)

        # Remove lat/lon vertices variables since they are not needed
        run_cmd(f'ncks -C -O  -x -v lat_vertices,lon_vertices {file_name} {file_name}',
                verbose,prefix='  ',suffix='',shell=True)

        # also remove "bounds" attribute
        run_cmd(f'ncatted -O -a bounds,lat,d,, {file_name} {file_name}',
                verbose,prefix='  ',suffix='',shell=True)
        run_cmd(f'ncatted -O -a bounds,lon,d,, {file_name} {file_name}',
                verbose,shell=True)

        return
    # --------------------------------------------------------------------------
    def clean_global_attributes(self,file_name,verbose=None):
        """ remove messy global attributes of the file """
        if verbose is None : verbose = hiccup_verbose
        if verbose: print('\nCleaning up excessive global attributes...')
        
        global_att_list = ['history_of_appended_files', 'nco_openmp_thread_number', 
                           'input_file', 'map_file', 'remap_version', 'remap_hostname', 
                           'remap_command', 'remap_script', 'NCO' ]
        
        check_dependency('ncatted')

        # Remove the attributes listed in global_att_list
        for att in global_att_list:
            run_cmd(f'ncatted -O -a {att},global,d,, {file_name} {file_name}',
                    verbose,prefix='  ',suffix='')

        # Also reset the history attribute
        run_cmd(f'ncatted -h -O -a history,global,o,c, {file_name} {file_name}',
                verbose,prefix='  ',suffix='')
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
