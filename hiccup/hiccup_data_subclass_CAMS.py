import os
import numpy as np
import xarray as xr
from time import perf_counter
from hiccup.hiccup_data_class import hiccup_data
from hiccup.hiccup_data_class import hiccup_verbose
from hiccup.hiccup_utilities import check_dependency
from hiccup.hiccup_utilities import run_cmd
from hiccup.hiccup_utilities import tcolor
# ------------------------------------------------------------------------------
# HICCUP subclass for CAMS atmospheric composition data (dust, aerosols, trace gases, etc.)
# ------------------------------------------------------------------------------
class CAMS(hiccup_data):
    @classmethod
    def is_name_for(cls,src_data_name) : return src_data_name == 'CAMS'
    def __init__( self, name,
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
                        )

        self.src_data_name = 'CAMS'
        self.lev_name = 'level'
        self.new_lev_name = 'plev'

        if self.target_model=='EAM':
            self.atm_var_name_dict.update({'lat':'latitude'})
            self.atm_var_name_dict.update({'lon':'longitude'})
            self.atm_var_name_dict.update({'dst_a1_1':'aermr04'})   # Dust Aerosol (0.03 - 0.55 um) Mixing Ratio
            self.atm_var_name_dict.update({'dst_a1_2':'aermr05'})   # Dust Aerosol (0.55 - 0.90 um) Mixing Ratio
            self.atm_var_name_dict.update({'dst_a3':'aermr06'})     # Dust Aerosol (0.90 - 20.0 um) Mixing Ratio
            self.sfc_var_name_dict.update({'PS':'sp'})              # sfc pressure

        # if self.target_model=='EAMXX':
        #     self.atm_var_name_dict.update({'lat':'latitude'})
        #     self.atm_var_name_dict.update({'lon':'longitude'})
        #     self.atm_var_name_dict.update({'dst_a1_1':'aermr04'})   # Dust Aerosol (0.03 - 0.55 um) Mixing Ratio
        #     self.atm_var_name_dict.update({'dst_a1_2':'aermr05'})   # Dust Aerosol (0.55 - 0.90 um) Mixing Ratio
        #     self.atm_var_name_dict.update({'dst_a3':'aermr06'})     # Dust Aerosol (0.90 - 20.0 um) Mixing Ratio
        #     self.sfc_var_name_dict.update({'PS':'sp'})              # sfc pressure

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
        if verbose is None : verbose = hiccup_verbose
        if verbose : print(verbose_indent+'\nGenerating src grid file...')

        # Remove the file here to prevent the warning message when ncremap overwrites it
        if os.path.isfile(self.src_grid_file): run_cmd(f'rm {self.src_grid_file} ',verbose)

        check_dependency('ncremap')

        cmd  = f'ncremap'
        cmd += f' --tmp_dir={self.tmp_dir}'
        cmd += f' -G ttl=\'Equi-Angular grid {self.src_horz_grid}\''
        cmd += f'#latlon={self.src_nlat},{self.src_nlon}'
        cmd +=  '#lat_typ=uni#lat_drc=n2s#lon_typ=grn_ctr '
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

        # delete lat/lon coordinates attributes, which can be problematic later
        if 'lat' in ds.coords: ds = ds.reset_coords(names='lat', drop=True)
        if 'lon' in ds.coords: ds = ds.reset_coords(names='lon', drop=True)

        if do_timers: self.print_timer(timer_start)
        return
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
