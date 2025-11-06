import os
import numpy as np
import xarray as xr
from time import perf_counter
from hiccup.hiccup_data_class import hiccup_data
from hiccup.hiccup_utilities import check_dependency
from hiccup.hiccup_utilities import run_cmd
from hiccup.hiccup_utilities import tcolor

# log file for Tempest output
tempest_log_file = 'TempestRemap.log'

# ------------------------------------------------------------------------------
# HICCUP subclass for EAM source data - i.e. repurpose existing model init data
# ------------------------------------------------------------------------------
class EAM(hiccup_data):
    @classmethod
    def is_name_for(cls,src_data_name) : return src_data_name == 'EAM'
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
        
        self.src_data_name = 'EAM'
        self.lev_name = 'lev'
        self.new_lev_name = 'lev'

        # set flag to indicate whether source data uses hybrid vertical coordinate
        self.src_hybrid_lev = True

        self.npg = 2

        self.dst_horz_grid_np = self.dst_horz_grid
        self.dst_horz_grid_pg = self.dst_horz_grid.replace('np4',f'pg{self.npg}')

        # Determine source grid from input atmosphere file
        ds = xr.open_dataset(atm_file)

        if check_input_files:
            if 'ne' in ds.attrs:
                ne = ds.attrs['ne']
            elif 'ncol' in ds.sizes:
                # use ncol formula to solve for # elements (ne): ncol_dyn = ne^2*6*9+2
                ne = int( np.sqrt( (ds.sizes['ncol']-2)/(6*9) ) )
            else:
                raise KeyError(f'Cannot determine source grid from atm_file: {atm_file}')
        else:
            ne = 0

        self.src_horz_grid    = f'ne{ne}np4'
        self.src_horz_grid_np = f'ne{ne}np4'
        self.src_horz_grid_pg = f'ne{ne}pg{self.npg}'

        ncol_size_np = np.square(ne)*6*9+2
        ncol_size_pg = np.square(ne)*6*np.square(self.npg)

        src_ne = self.get_src_grid_ne()
        dst_ne = self.get_dst_grid_ne()

        self.src_grid_file_np = f'{self.grid_dir}/exodus_ne{src_ne}.g'
        self.src_grid_file_pg = f'{self.grid_dir}/scrip_{ self.src_horz_grid_pg}.nc'

        self.dst_grid_file_np = f'{self.grid_dir}/exodus_ne{dst_ne}.g'
        self.dst_grid_file_pg = f'{self.grid_dir}/scrip_{ self.dst_horz_grid_pg}.nc'

        self.map_file_np = f'{self.map_dir}/map_{self.src_horz_grid_np}_to_{self.dst_horz_grid_np}.nc'
        self.map_file_pg = f'{self.map_dir}/map_{self.src_horz_grid_pg}_to_{self.dst_horz_grid_pg}.nc'


        # Atmospheric variables - need separate treatment for np4 and pgN data
        self.atm_var_name_dict_np = {}
        self.atm_var_name_dict_pg = {}
        ds = xr.open_dataset(self.atm_file)

        # check that dimension sizes are consistent
        if 'ncol' in ds.dims:
            if len(ds['ncol'])!=ncol_size_np and len(ds['ncol'])!=ncol_size_pg:
                raise ValueError('input data file does not have expected dimension sizes')
        if 'ncol_d' in ds.dims:
            if len(ds['ncol_d'])!=ncol_size_np and len(ds['ncol_d'])!=ncol_size_pg:
                raise ValueError('input data file does not have expected dimension sizes')

        if self.target_model=='EAM':
            for key in ds.variables.keys(): 
                if key in ['lat','lon','lat_d','lon_d','lat_vertices','lon_vertices']:
                    continue
                if 'ncol_d' in ds[key].dims: 
                    self.atm_var_name_dict_np.update({key:key})
                if 'ncol' in ds[key].dims:
                    if ds.sizes['ncol']==ncol_size_np:
                        self.atm_var_name_dict_np.update({key:key})
                    if ds.sizes['ncol']==ncol_size_pg:
                        self.atm_var_name_dict_pg.update({key:key})

        if self.target_model=='EAMXX':
            self.atm_var_name_dict = {}
            self.sfc_var_name_dict = {}
            self.atm_var_name_dict.update({'lat':'lat'})
            self.atm_var_name_dict.update({'lon':'lon'})
            self.atm_var_name_dict.update({'T_mid':'T'})                # temperature
            self.atm_var_name_dict.update({'qv':'Q'})                   # specific humidity
            self.atm_var_name_dict.update({'horiz_winds_u':'U'})        # zonal wind
            self.atm_var_name_dict.update({'horiz_winds_v':'V'})        # meridional wind
            self.atm_var_name_dict.update({'o3_volume_mix_ratio':'O3'}) # ozone mass mixing ratio
            self.atm_var_name_dict.update({'qc':'CLDLIQ'})              # specific cloud liq water
            self.atm_var_name_dict.update({'qi':'CLDICE'})              # specific cloud ice water
            self.sfc_var_name_dict.update({'ps':'PS'})                  # sfc pressure
            self.sfc_var_name_dict.update({'phis':'PHIS'})              # surface geopotential

    # --------------------------------------------------------------------------
    def create_src_grid_file(self,verbose=None):
        """ 
        Generate source grid file 
        """
        if self.do_timers: timer_start = perf_counter()
        if verbose is None : verbose = self.verbose
        if verbose : print(self.verbose_indent+'\nGenerating src grid files (np+pg)...')

        # Remove the file here to prevent the warning message when ncremap overwrites it
        if self.src_grid_file is not None:
            if os.path.isfile(self.src_grid_file): run_cmd(f'rm {self.src_grid_file} ',verbose)

        check_dependency('ncremap')

        # Spectral element grid
        ne = self.get_src_grid_ne()
        npg = self.get_src_grid_npg()

        check_dependency('GenerateCSMesh')
        cmd = f'GenerateCSMesh --alt --res {ne} --file {self.src_grid_file_np}'
        cmd += f' >> {tempest_log_file}'
        run_cmd(cmd,verbose,shell=True)

        # Next switch to volumetric mesh that matches the physgrid
        tmp_exodus_file = f'{self.grid_dir}/exodus_{self.src_horz_grid_pg}.g'
        check_dependency('GenerateVolumetricMesh')
        cmd = 'GenerateVolumetricMesh'
        cmd += f' --in {self.src_grid_file_np} '
        cmd += f' --out {tmp_exodus_file} '
        cmd += f' --np {npg} --uniform'
        cmd += f' >> {tempest_log_file}'
        run_cmd(cmd,verbose,shell=True)

        # Create pgN scrip file
        check_dependency('ConvertMeshToSCRIP')
        scrip_file = f'{self.grid_dir}/scrip_{self.dst_horz_grid_pg}.nc'
        cmd = 'ConvertMeshToSCRIP'
        cmd += f' --in {tmp_exodus_file} '
        cmd += f' --out {self.src_grid_file_pg} '
        cmd += f' >> {tempest_log_file}'
        run_cmd(cmd,verbose,shell=True)

        # fix grid_imask type
        run_cmd(f'ncap2 --overwrite -s \'grid_imask=int(grid_imask)\' '
                +f'{self.src_grid_file_pg} {self.src_grid_file_pg}',verbose,shell=True)

        # delete temporary exodus file
        run_cmd(f'rm {tmp_exodus_file} ',verbose)

        if self.do_timers: self.print_timer(timer_start)
        return 
    # --------------------------------------------------------------------------
    def create_dst_grid_file(self,verbose=None):
        """ 
        Generate destination model grid file. For the case where EAM data 
        is the source we need both np4 and pgN grids
        """
        if self.do_timers: timer_start = perf_counter()
        if verbose is None : verbose = self.verbose
        if verbose : print(self.verbose_indent+'\nGenerating dst grid files (np+pg)...')
        
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

        # fix grid_imask type
        run_cmd('ncap2 --overwrite -s \'grid_imask=int(grid_imask)\' '
                +f'{self.dst_grid_file_pg} {self.dst_grid_file_pg}',verbose,shell=True)

        # delete temporary exodus file
        run_cmd(f'rm {tmp_exodus_file} ',verbose)

        if self.do_timers: self.print_timer(timer_start)
        return 
    # --------------------------------------------------------------------------
    def create_map_file(self,verbose=None):
        """ 
        Generate mapping files for EAM after grid files have been created 
        (overloads default routine that assumes only one map file is needed)
        """
        if self.do_timers: timer_start = perf_counter()
        if verbose is None : verbose = self.verbose
        if verbose : print(self.verbose_indent+'\nGenerating mapping files (np+pg)...')

        check_dependency('ncremap')

        dst_ne = self.get_dst_grid_ne()
        src_ne = self.get_src_grid_ne()

        # Check that grid file fields are not empty
        if self.src_grid_file_np is None : raise ValueError('src_grid_file_np is not defined!')
        if self.src_grid_file_pg is None : raise ValueError('src_grid_file_pg is not defined!')
        if self.dst_grid_file_np is None : raise ValueError('dst_grid_file_np is not defined!')
        if self.dst_grid_file_pg is None : raise ValueError('dst_grid_file_pg is not defined!')

        # Create the np4 map file
        cmd = f'ncremap -a se2se '
        cmd += f' --src_grd={self.src_grid_file_np}'
        cmd += f' --dst_grd={self.dst_grid_file_np}'
        cmd += f' --map_file={self.map_file_np}'
        if dst_ne>src_ne : cmd += ' --lrg2sml '
        run_cmd(cmd,verbose,shell=True)

        # Create the pgN map file
        cmd = f'ncremap -a fv2fv_flx '
        cmd += f' --src_grd={self.src_grid_file_pg}'
        cmd += f' --dst_grd={self.dst_grid_file_pg}'
        cmd += f' --map_file={self.map_file_pg}'
        if dst_ne>src_ne : cmd += ' --lrg2sml '
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
        return
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
