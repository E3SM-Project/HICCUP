from hiccup.hiccup_data_class_common import __all__
from hiccup.hiccup_data_class import hiccup_data
from hiccup.hiccup_utilities import check_dependency
from hiccup.hiccup_utilities import run_cmd
from hiccup.hiccup_utilities import tcolor
# ------------------------------------------------------------------------------
# HICCUP subclass for NOAA source data
# ------------------------------------------------------------------------------
class NOAA(hiccup_data):
    @classmethod
    def is_name_for(cls,src_data_name) : return src_data_name == 'NOAA'
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
        
        self.src_data_name = 'NOAA'
        self.sstice_name=='NOAA'
        self.sst_name = 'sst'
        self.ice_name = 'icec'

        # set flag to indicate whether source data uses hybrid vertical coordinate
        self.src_hybrid_lev = False

        self.ds_sst = xr.open_dataset(self.sst_file)
        self.ds_ice = xr.open_dataset(self.sst_file)

        self.src_nlat = len( self.ds_sst['latitude'].values )
        self.src_nlon = len( self.ds_sst['longitude'].values )

        self.src_horz_grid = f'{self.src_nlat}x{self.src_nlon}'
        self.src_grid_file = f'{self.grid_dir}/scrip_{self.src_data_name}_{self.src_horz_grid}.nc'

        self.map_file = f'{self.map_dir}/map_{self.src_horz_grid}_to_{self.dst_horz_grid}.nc'

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
