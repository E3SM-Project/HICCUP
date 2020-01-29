# Class structure for defining parameters specific to the 
# input data used for creating initial condition files
# 
# NOTE: Variable name dictionaries are defined with the key as the model's 
# variable name and the value as the reanalysis data variable name

import xarray as xr
import subprocess as sp
import shutil 
import re

ncremap_alg = ' --alg_typ=tempest '        # algorithm flag for ncremap

#-------------------------------------------------------------------------------
# Method for checking if required software is installed
#-------------------------------------------------------------------------------
def check_dependency(cmd):
    """ Check for required system commands"""
    if shutil.which(cmd) is None : raise OSError(f'{cmd} is not in system path')
    return
#-------------------------------------------------------------------------------
# Method for returning class object
#-------------------------------------------------------------------------------
def create_hiccup_data(name,atm_file,sfc_file,dst_grid_name,lev_type=''):
    """ Return HICCUP data class object """
    check_requirements()
    for subclass in hiccup_data.__subclasses__():
        if subclass.is_name_for(name):
            return subclass(name
                           ,atm_file=atm_file
                           ,sfc_file=sfc_file
                           ,dst_grid_name=dst_grid_name
                           ,lev_type=lev_type)
    raise ValueError(f'{name} is not a valid HICCUP dataset name')
#-------------------------------------------------------------------------------
# Base Class
#-------------------------------------------------------------------------------
class hiccup_data(object):

    def __init__(self,name,lev_type='',atm_file='',sfc_file='',dst_grid_name=''):
        self.name = name
        self.lev_type = lev_type
        self.atm_file = atm_file
        self.sfc_file = sfc_file
        self.atm_var_name_dict = {}
        self.sfc_var_name_dict = {}
        self.lnd_var_name_dict = {}
        self.nlat = -1
        self.nlon = -1
        self.dst_grid_name = dst_grid_name
        self.src_grid_name = ''
        self.src_grid_file = None
        self.dst_grid_file = None
        self.map_file = None

        # Load input files into xarray datasets
        self.ds_atm = xr.open_dataset(self.atm_file)
        self.ds_sfc = xr.open_dataset(self.sfc_file)
    #---------------------------------------------------------------------------
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
    #---------------------------------------------------------------------------
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
    #---------------------------------------------------------------------------
    def create_dst_grid_file(self):
        """ Generate destination model grid file """
        
        if 'ne' in self.dst_grid_name and 'np' in self.dst_grid_name : 
            # Spectral element grid
            ne = re.search('ne(.*)np', self.dst_grid_name).group(1)
            self.dst_grid_file = f'exodus_ne{ne}.g'
            check_dependency('GenerateCSMesh')
            cmd  = f'GenerateCSMesh --res {ne} --file {self.dst_grid_file}'
            sp.call(cmd, shell=True)

        elif 'ne' in self.dst_grid_name and 'pg' in self.dst_grid_name : 
            # Spectral element grid with FV physics grid (ex. ne30pg2)
            ne  = re.search('ne(.*)pg', self.dst_grid_name).group(1)
            npg = re.search('pg(.*)', self.dst_grid_name).group(1)
            # First create exodus file
            exodus_file = f'exodus_ne{ne}.g'
            check_dependency('GenerateCSMesh')
            cmd  = f'GenerateCSMesh --res {ne} --file {exodus_file}'
            sp.call(cmd, shell=True)
            # Next create script file for FV physgrid
            self.dst_grid_file = f'scrip_{self.dst_grid_name}.nc'
            check_dependency('GenerateVolumetricMesh')
            cmd = f'GenerateVolumetricMesh --in {exodus_file} --out {self.dst_grid_file} --np {npg} --uniform'
            sp.call(cmd, shell=True)

        else:
            raise ValueError(f'grid_name={self.dst_grid_name} is not currently supported')

        return 
    #---------------------------------------------------------------------------
    def create_map_file(self):
        """ Generate mapping file aftergrid files have been created """
        if self.src_grid_file == None : 
            raise ValueError('src_grid_file is not defined for hiccup_data object')
        if self.dst_grid_file == None : 
            raise ValueError('dst_grid_file is not defined for hiccup_data object')
        self.map_file = f'map_{self.src_grid_name}_to_{self.dst_grid_name}.nc'
        # self.map_opts = '--in_type fv --in_np 1 --mono --out_format Classic '
        self.map_opts = '--in_type fv --in_np 1 --mono --out_double '
        if 'ne' in self.dst_grid_name and 'np' in self.dst_grid_name : 
            self.map_opts = self.map_opts+' --out_type cgll --out_np 4 ' # options for SE grid
        else:
            self.map_opts = self.map_opts+' --out_type fv --out_np 1 --volumetric '
        check_dependency('ncremap')
        cmd = f'ncremap {ncremap_alg} '
        cmd = cmd+f' --src_grd={self.src_grid_file}'
        cmd = cmd+f' --dst_grd={self.dst_grid_file}'
        cmd = cmd+f' --map_file={self.map_file}'
        cmd = cmd+f' --wgt_opt=\'{self.map_opts}\' '
        print(f'\n{cmd}\n')
        sp.call(cmd, shell=True)
        return
#-------------------------------------------------------------------------------
# Subclasses
#-------------------------------------------------------------------------------
class ERA5(hiccup_data):
    @classmethod
    def is_name_for(cls,name) : return name == 'ERA5'
    def __init__(self,name,atm_file,sfc_file,dst_grid_name,lev_type=''):
        super().__init__(name,atm_file=atm_file
                        ,sfc_file=sfc_file
                        ,dst_grid_name=dst_grid_name
                        ,lev_type=lev_type)
        
        self.name = 'ERA5'

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
        self.sfc_var_name_dict.update({'SST':'sst'})       # sea sfc temperature 
        self.sfc_var_name_dict.update({'TS1':'stl1'})      # Soil temperature level 1 
        self.sfc_var_name_dict.update({'TS2':'stl2'})      # Soil temperature level 2 
        self.sfc_var_name_dict.update({'TS3':'stl3'})      # Soil temperature level 3 
        self.sfc_var_name_dict.update({'TS4':'stl4'})      # Soil temperature level 4 
        self.sfc_var_name_dict.update({'':'lai_hv'})       # Leaf area index, high vegetation 
        self.sfc_var_name_dict.update({'':'lai_lv'})       # Leaf area index, low vegetation 
        self.sfc_var_name_dict.update({'':'src'})          # Skin reservoir content 
        self.sfc_var_name_dict.update({'':'asn'})          # Snow albedo 
        self.sfc_var_name_dict.update({'':'rsn'})          # Snow density 
        self.sfc_var_name_dict.update({'':'sd'})           # Snow depth 
        self.sfc_var_name_dict.update({'':'tsn'})          # Temperature of snow layer 
        self.sfc_var_name_dict.update({'':'swvl1'})        # Volumetric soil water level 1 
        self.sfc_var_name_dict.update({'':'swvl2'})        # Volumetric soil water level 2 
        self.sfc_var_name_dict.update({'':'swvl3'})        # Volumetric soil water level 3 
        self.sfc_var_name_dict.update({'':'swvl4'})        # Volumetric soil water level 4 

        # Land model variables
        # self.lnd_var_name_dict.update({'TS1':'stl1'})      # Soil temperature level 1 
        # self.lnd_var_name_dict.update({'TS2':'stl2'})      # Soil temperature level 2 
        # self.lnd_var_name_dict.update({'TS3':'stl3'})      # Soil temperature level 3 
        # self.lnd_var_name_dict.update({'TS4':'stl4'})      # Soil temperature level 4 
        # self.lnd_var_name_dict.update({'':'lai_hv'})       # Leaf area index, high vegetation 
        # self.lnd_var_name_dict.update({'':'lai_lv'})       # Leaf area index, low vegetation 
        # self.lnd_var_name_dict.update({'':'src'})          # Skin reservoir content 
        # self.lnd_var_name_dict.update({'':'asn'})          # Snow albedo 
        # self.lnd_var_name_dict.update({'':'snowc'})        # Snow cover 
        # self.lnd_var_name_dict.update({'':'rsn'})          # Snow density 
        # self.lnd_var_name_dict.update({'':'sde'})          # Snow depth (liq water equivalent) 
        # self.lnd_var_name_dict.update({'':'sd'})           # Snow depth 
        # self.lnd_var_name_dict.update({'':'tsn'})          # Temperature of snow layer 
        # self.lnd_var_name_dict.update({'':'swvl1'})        # Volumetric soil water level 1 
        # self.lnd_var_name_dict.update({'':'swvl2'})        # Volumetric soil water level 2 
        # self.lnd_var_name_dict.update({'':'swvl3'})        # Volumetric soil water level 3 
        # self.lnd_var_name_dict.update({'':'swvl4'})        # Volumetric soil water level 4 

        self.nlat = len( self.ds_atm[ self.atm_var_name_dict['lat'] ].values )
        self.nlon = len( self.ds_atm[ self.atm_var_name_dict['lon'] ].values )
    #---------------------------------------------------------------------------
    def create_src_grid_file(self):
        """ Generate source grid file """
        self.src_grid_name = f'{self.nlat}x{self.nlon}'
        self.src_grid_file = f'scrip_{self.name}_{self.src_grid_name}.nc'

        check_dependency('ncremap')
        cmd  = f'ncremap {ncremap_alg} ' \
              +f' -G ttl=\'Equi-Angular grid {self.src_grid_name}\''     \
              +f'#latlon={self.nlat},{self.nlon}'                           \
              +f'#lat_typ=uni'                                              \
              +f'#lon_typ=grn_ctr '                                         \
              +f' -g {self.src_grid_file} '
        sp.call(cmd, shell=True)
        return 
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


