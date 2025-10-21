import os
import numpy as np
import xarray as xr
from time import perf_counter
from hiccup.hiccup_data_class import hiccup_data
from hiccup.hiccup_utilities import check_dependency
from hiccup.hiccup_utilities import run_cmd
from hiccup.hiccup_utilities import tcolor
# ------------------------------------------------------------------------------
# HICCUP subclass for IFS/ECMWF source data
# ------------------------------------------------------------------------------
class IFS(hiccup_data):
    @classmethod
    def is_name_for(cls,src_data_name) : return src_data_name == 'IFS'
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

        self.src_data_name = 'IFS'
        self.new_lev_name = 'plev'
        # set vertical coordinate name to accommodate old and new names from CDS upgrade
        self.lev_name = None
        if self.atm_file is not None:
            if 'pressure_level' in self.ds_atm.coords: self.lev_name = 'pressure_level'
            if 'level'          in self.ds_atm.coords: self.lev_name = 'level'
            if self.lev_name is None:
                raise ValueError('IFS subclass: lev_name cannot be set from input data coordinates')

        # set flag to indicate whether source data uses hybrid vertical coordinate
        self.src_hybrid_lev = False
        
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
            self.atm_var_name_dict.update({'lnsp':'lnsp'})      # logarithm of sfc pressure

            # Surface variables
            self.sfc_var_name_dict.update({'PS':'sp'})         # Logarithm of sfc pressure
            self.sfc_var_name_dict.update({'TS':'skt'})        # skin temperature
            self.sfc_var_name_dict.update({'PHIS':'z'})        # surface geopotential

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
        # we're going to hijack this routine in order to inject the vertical
        # grid coordinate data that we need for vertical remapping IFS data
        
        ds_lev = xr.open_dataset('files_vert/vrt_hyb_ifs_L138.nc')

        print(); print(ds)
        print(); print(ds_lev)
        exit()

        return
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
