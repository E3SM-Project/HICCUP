# Class structure for defining parameters specific to the 
# input data used for creating initial condition files
# 
# NOTE: Variable name dictionaries are defined with the key as the model's 
# variable name and the value as the reanalysis data variable name

import xarray as xr
import subprocess as sp
import shutil 

ncremap_alg = ' -a tempest '        # algorithm flag for ncremap

#-------------------------------------------------------------------------------
# Method for checking if required software is installed
#-------------------------------------------------------------------------------
def check_requirements():
    """ Check for required system commands"""
    cmd_list = []
    cmd_list.append('GenerateCSMesh')
    cmd_list.append('GenerateVolumetricMesh')
    for cmd in cmd_list:
        if shutil.which(cmd) is None :
            raise OSError(f'{cmd} is not in system path')
    return
#-------------------------------------------------------------------------------
# Method for returning class object
#-------------------------------------------------------------------------------
def create_hiccup_data(name,lev_type='',atm_file='',sfc_file=''):
    """ Return HICCUP data class object """
    check_requirements()
    for subclass in hiccup_data.__subclasses__():
        if subclass.is_name_for(name):
            return subclass(name,lev_type=lev_type,atm_file=atm_file,sfc_file=sfc_file)
    raise ValueError(f'{name} is not a valid HICCUP dataset name')
#-------------------------------------------------------------------------------
# Base Class
#-------------------------------------------------------------------------------
class hiccup_data(object):

    def __init__(self,name,lev_type='',atm_file='',sfc_file=''):
        self.name = name
        self.lev_type = lev_type
        self.atm_file = atm_file
        self.sfc_file = sfc_file
        self.atm_var_name_dict = {}
        self.sfc_var_name_dict = {}
        self.lnd_var_name_dict = {}
        self.nlat = -1
        self.nlon = -1

        self.ds_atm = xr.open_dataset(self.atm_file)
        self.ds_sfc = xr.open_dataset(self.sfc_file)

        

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

    def create_dst_grid_file(self,grid):
        # Generate source grid file:
        self.grid_file = f'scrip_{self.grid}.nc'
        cmd  = f'ncremap {ncremap_alg}' \
              +f' -G ttl=\'Equi-Angular grid {self.nlat}x{self.nlon}\''     \
              +f'#latlon={self.nlat},{self.nlon}'                           \
              +f'#lat_typ=uni'                                              \
              +f'#lon_typ=grn_ctr '                                         \
              +f' -g {self.grid_file} '
        sp.call(cmd, shell=True)
        return 
#-------------------------------------------------------------------------------
# Subclasses
#-------------------------------------------------------------------------------
class ERA5(hiccup_data):
    @classmethod
    def is_name_for(cls,name) : return name == 'ERA5'
    def __init__(self,name,lev_type='',atm_file='',sfc_file=''):
        super().__init__(name,lev_type=lev_type,atm_file=atm_file,sfc_file=sfc_file)
        
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

    def create_src_grid_file(self):
        # Generate source grid file:
        self.grid_file = f'scrip_{self.name}_{self.nlat}x{self.nlon}.nc'
        cmd  = f'ncremap {ncremap_alg}' \
              +f' -G ttl=\'Equi-Angular grid {self.nlat}x{self.nlon}\''     \
              +f'#latlon={self.nlat},{self.nlon}'                           \
              +f'#lat_typ=uni'                                              \
              +f'#lon_typ=grn_ctr '                                         \
              +f' -g {self.grid_file} '
        sp.call(cmd, shell=True)
        return 
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


