import os
import numpy as np
import xarray as xr
from time import perf_counter
from hiccup.hiccup_data_class import hiccup_data
from hiccup.hiccup_utilities import check_dependency
from hiccup.hiccup_utilities import run_cmd
from hiccup.hiccup_utilities import tcolor
# ------------------------------------------------------------------------------
# HICCUP subclass for ERA5 source data
# ------------------------------------------------------------------------------
class ERA5(hiccup_data):
    @classmethod
    def is_name_for(cls,src_data_name) : return src_data_name == 'ERA5'
    def __init__( self, src_data_name,
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
        super().__init__(   
                          target_model=target_model,
                          atm_file=atm_file,
                          sfc_file=sfc_file,
                          dst_horz_grid=dst_horz_grid,
                          dst_vert_grid=dst_vert_grid,
                          output_dir=output_dir,
                          grid_dir=grid_dir,
                          map_dir=map_dir,
                          tmp_dir=tmp_dir,
                          sstice_combined_file=sstice_combined_file,
                          sstice_name=sstice_name,
                          topo_file=topo_file,
                          sst_file=sst_file,
                          ice_file=ice_file,
                          lev_type=lev_type,
                          check_input_files=check_input_files,
                          RRM_grid=RRM_grid,
                          do_timers=do_timers,
                          verbose=verbose,
                          verbose_indent=verbose_indent,
                        )

        self.src_data_name = 'ERA5'
        self.lev_name = 'level'
        self.new_lev_name = 'plev'

        if self.target_model=='EAM':
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

        if self.target_model=='EAMXX':
            self.atm_var_name_dict.update({'lat':'latitude'})
            self.atm_var_name_dict.update({'lon':'longitude'})
            self.atm_var_name_dict.update({'T_mid':'t'})                # temperature
            self.atm_var_name_dict.update({'qv':'q'})                   # specific humidity
            self.atm_var_name_dict.update({'horiz_winds_u':'u'})        # zonal wind
            self.atm_var_name_dict.update({'horiz_winds_v':'v'})        # meridional wind
            self.atm_var_name_dict.update({'qc':'clwc'})                # specific cloud liq water
            self.atm_var_name_dict.update({'qi':'ciwc'})                # specific cloud ice water
            self.atm_var_name_dict.update({'o3_volume_mix_ratio':'o3'}) # ozone mass mixing ratio
            self.sfc_var_name_dict.update({'ps':'sp'})                  # sfc pressure
            self.sfc_var_name_dict.update({'phis':'z'})                 # surface geopotential

        if self.target_model=='EAMXX-nudging':
            self.atm_var_name_dict.update({'lat':'latitude'})
            self.atm_var_name_dict.update({'lon':'longitude'})
            self.sfc_var_name_dict.update({'PS':'sp'})          # sfc pressure
            self.atm_var_name_dict.update({'U':'u'})            # zonal wind
            self.atm_var_name_dict.update({'V':'v'})            # meridional wind

            if not RRM_grid:
                dst_ne = self.get_dst_grid_ne()
                self.npg = 2
                self.dst_horz_grid_np = self.dst_horz_grid.replace(f'pg{self.npg}','np4')
                self.dst_horz_grid_pg = self.dst_horz_grid.replace('np4',f'pg{self.npg}')
                self.dst_grid_file_np = f'{self.grid_dir}/exodus_ne{dst_ne}.g'
                self.dst_grid_file_pg = f'{self.grid_dir}/scrip_{ self.dst_horz_grid_pg}.nc'

        self.src_nlat = len( self.ds_atm[ self.atm_var_name_dict['lat'] ].values )
        self.src_nlon = len( self.ds_atm[ self.atm_var_name_dict['lon'] ].values )

        self.src_horz_grid = f'{self.src_nlat}x{self.src_nlon}'
        self.src_grid_file = f'{self.grid_dir}/scrip_{self.src_data_name}_{self.src_horz_grid}.nc'

        self.map_file = f'{self.map_dir}/map_{self.src_horz_grid}_to_{self.dst_horz_grid}.nc'

    # --------------------------------------------------------------------------
    def create_src_grid_file(self,verbose=None):
        """ 
        Generate source grid file 
        """
        if self.do_timers: timer_start = perf_counter()
        if verbose is None : verbose = self.verbose
        if verbose : print(self.verbose_indent+'\nGenerating src grid file...')

        # Remove the file here to prevent the warning message when ncremap overwrites it
        if os.path.isfile(self.src_grid_file): run_cmd(f'rm {self.src_grid_file} ',verbose)

        check_dependency('ncremap')

        cmd  = f'ncremap'
        cmd += f' --tmp_dir={self.tmp_dir}'
        cmd += f' -G ttl=\'Equi-Angular grid {self.src_horz_grid}\'' 
        cmd += f'#latlon={self.src_nlat},{self.src_nlon}'                    
        cmd +=  '#lat_typ=uni'
        cmd +=  '#lat_drc=n2s'
        cmd +=  '#lon_typ=grn_ctr '
        cmd += f' -g {self.src_grid_file} '
        run_cmd(cmd,verbose,shell=True)

        if self.do_timers: self.print_timer(timer_start)
        return
    # --------------------------------------------------------------------------
    def rename_vars_special(self,ds,verbose=None,do_timers=None
                           ,new_lev_name=None,change_pressure_name=True
                           ,adjust_pressure_units=True):
        """ 
        Rename file vars specific to this subclass 
        """
        if do_timers is None: do_timers = self.do_timers
        if do_timers: timer_start = perf_counter()
        if verbose is None : verbose = self.verbose

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

        # delete lat/lon coordinates attributes, which can be problematic later
        if 'lat' in ds.coords: ds = ds.reset_coords(names='lat', drop=True)
        if 'lon' in ds.coords: ds = ds.reset_coords(names='lon', drop=True)

        # for v in ds.variables:
        #     if 'eulaVlliF_' in ds[v].attrs: del ds[v].attrs['eulaVlliF_']
        
        if do_timers: self.print_timer(timer_start)
        return
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
