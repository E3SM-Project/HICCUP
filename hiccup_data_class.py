# Class structure for defining parameters specific to the 
# input data used for creating initial condition files
# 
# NOTE: Variable name dictionaries are defined with the key as the model's 
# variable name and the value as the reanalysis data variable name
#-------------------------------------------------------------------------------
# Method for returning class object
#-------------------------------------------------------------------------------
def create_hiccup_data(name):
  for subclass in hiccup_data.__subclasses__():
    if subclass.is_name_for(name):
      return subclass(name)
  raise ValueError(f'{name} is not a valid HICCUP dataset name')
#-------------------------------------------------------------------------------
# Base Class
#-------------------------------------------------------------------------------
class hiccup_data(object):

    def __init__(self,name,lev_type=''):
        self.name = name
        self.lev_type = lev_type
        self.atm_var_name_dict = {}
        self.sfc_var_name_dict = {}
        self.lnd_var_name_dict = {}

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

#-------------------------------------------------------------------------------
# Subclasses
#-------------------------------------------------------------------------------
class ERA5(hiccup_data):
    @classmethod
    def is_name_for(cls,name) : return name == 'ERA5'
    def __init__(self,name):
        super().__init__(name)
        
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
        self.lnd_var_name_dict.update({'TS1':'stl1'})      # Soil temperature level 1 
        self.lnd_var_name_dict.update({'TS2':'stl2'})      # Soil temperature level 2 
        self.lnd_var_name_dict.update({'TS3':'stl3'})      # Soil temperature level 3 
        self.lnd_var_name_dict.update({'TS4':'stl4'})      # Soil temperature level 4 
        self.lnd_var_name_dict.update({'':'lai_hv'})       # Leaf area index, high vegetation 
        self.lnd_var_name_dict.update({'':'lai_lv'})       # Leaf area index, low vegetation 
        self.lnd_var_name_dict.update({'':'src'})          # Skin reservoir content 
        self.lnd_var_name_dict.update({'':'asn'})          # Snow albedo 
        self.lnd_var_name_dict.update({'':'snowc'})        # Snow cover 
        self.lnd_var_name_dict.update({'':'rsn'})          # Snow density 
        self.lnd_var_name_dict.update({'':'sde'})          # Snow depth (liq water equivalent) 
        self.lnd_var_name_dict.update({'':'sd'})           # Snow depth 
        self.lnd_var_name_dict.update({'':'tsn'})          # Temperature of snow layer 
        self.lnd_var_name_dict.update({'':'swvl1'})        # Volumetric soil water level 1 
        self.lnd_var_name_dict.update({'':'swvl2'})        # Volumetric soil water level 2 
        self.lnd_var_name_dict.update({'':'swvl3'})        # Volumetric soil water level 3 
        self.lnd_var_name_dict.update({'':'swvl4'})        # Volumetric soil water level 4 


# class DYAMOND(hiccup_data):
#     @classmethod
#     def is_name_for(cls,name) : return name == 'DYAMOND'
#     def __init__(self,name):
#         super().__init__(name)
#         # DYAMOND - hybrid coords renamed because they use a different dimension name
#         self.var_name_dict.update({'hyam':'ohyam'})
#         self.var_name_dict.update({'hybm':'ohybm'})
#         self.var_name_dict.update({'hyai':'ohyai'})
#         self.var_name_dict.update({'hybi':'ohybi'})
#         self.var_name_dict.update({'t':'T'})
#         self.var_name_dict.update({'q':'Q'})
#         self.var_name_dict.update({'u':'U'})
#         self.var_name_dict.update({'v':'V'})
#         self.var_name_dict.update({'clwc':'CLDLIQ'})
#         self.var_name_dict.update({'ciwc':'CLDICE'})
#         self.var_name_dict.update({'z_2':'PHIS'})
#         self.var_name_dict.update({'o3':'O3'})
#         self.var_name_dict.update({'stl1':'TS1'})
#         self.var_name_dict.update({'stl2':'TS2'})
#         self.var_name_dict.update({'stl3':'TS3'})
#         self.var_name_dict.update({'stl4':'TS4'})

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


