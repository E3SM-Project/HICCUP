# Class structure for defining parameters specific to the 
# input data used for creating initial condition files
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
        self.var_name_dict = {}

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
        self.var_name_dict.update({'hyam':'ohyam'})
        self.var_name_dict.update({'hybm':'ohybm'})
        self.var_name_dict.update({'hyai':'ohyai'})
        self.var_name_dict.update({'hybi':'ohybi'})
        self.var_name_dict.update({'t':'T'})
        self.var_name_dict.update({'q':'Q'})
        self.var_name_dict.update({'u':'U'})
        self.var_name_dict.update({'v':'V'})
        self.var_name_dict.update({'clwc':'CLDLIQ'})
        self.var_name_dict.update({'ciwc':'CLDICE'})
        self.var_name_dict.update({'z_2':'PHIS'})
        self.var_name_dict.update({'o3':'O3'})
        self.var_name_dict.update({'stl1':'TS1'})
        self.var_name_dict.update({'stl2':'TS2'})
        self.var_name_dict.update({'stl3':'TS3'})
        self.var_name_dict.update({'stl4':'TS4'})


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


