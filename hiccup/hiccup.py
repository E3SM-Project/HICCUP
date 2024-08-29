import os
from time import perf_counter
import numpy as np
# ------------------------------------------------------------------------------
import hiccup.hiccup_data_class as hdc
import hiccup.hiccup_utilities as hu
from hiccup.hiccup_utilities import check_dependency
from hiccup.hiccup_utilities import run_cmd
from hiccup.hiccup_utilities import tcolor
# from hiccup.hiccup_utilities import print_timer_summary
from hiccup.hiccup_data_class_timer_methods import timer_start_total
# ------------------------------------------------------------------------------

default_target_model = 'EAM'

# default output paths
default_output_dir  = './data'
default_grid_dir    = './files_grid'
default_map_dir     = './files_mapping'
default_tmp_dir     = './files_tmp'

# Global verbosity default
hiccup_verbose = False
verbose_indent = ''

hdc.hiccup_verbose = hiccup_verbose
hdc.verbose_indent = verbose_indent

hu.hiccup_verbose = hiccup_verbose
hu.verbose_indent = verbose_indent

# Set numpy to ignore overflow errors
np.seterr(over='ignore')

# Disable HDF file locking to prevent permission 
# errors when writing data to files
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

# ------------------------------------------------------------------------------
# Get default topography file name
# ------------------------------------------------------------------------------
def get_default_topo_file_name(grid,topo_file_root=None):
    """
    return default topo file associated with input grid name
    """
    root_err_msg = 'No default for topo_file_root! Topo file root must be manually specified.'
    file_err_msg = 'No default for topo_file_name! Topo file path must be manually specified.'

    # Set root directory path if not provided
    if topo_file_root is None:
        nersc_inputdata_path = '/global/cfs/projectdirs/e3sm/inputdata'
        olcf_inputdata_path  = '/gpfs/alpine/cli115/world-shared/e3sm/inputdata'
        if os.path.exists(nersc_inputdata_path): topo_file_root = f'{nersc_inputdata_path}/atm/cam/topo'
        if os.path.exists(olcf_inputdata_path) : topo_file_root = f'{ olcf_inputdata_path}/atm/cam/topo'
        if topo_file_root is None: raise ValueError(root_err_msg)
    
    # Set default topo file
    topo_file_name = None
    if grid=='ne1024np4': topo_file_name = f'{topo_file_root}/USGS-gtopo30_ne1024np4_16xconsistentSGH_20190528.nc'
    if grid=='ne256np4' : topo_file_name = f'{topo_file_root}/USGS-gtopo30_ne256np4pg2_16xdel2_20200213.nc'
    if grid=='ne120np4' : topo_file_name = f'{topo_file_root}/USGS-gtopo30_ne120np4_16xdel2-PFC-consistentSGH.nc'
    if grid=='ne45np4'  : topo_file_name = f'{topo_file_root}/USGS-gtopo30_ne45np4pg2_16xdel2.c20200615.nc'
    if grid=='ne30np4'  : topo_file_name = f'{topo_file_root}/USGS-gtopo30_ne30np4pg2_16xdel2.c20200108.nc'
    if grid=='ne16np4'  : topo_file_name = f'{topo_file_root}/USGS-gtopo30_ne16np4pg2_16xdel2_20200527.nc'
    if grid=='ne11np4'  : topo_file_name = f'{topo_file_root}/USGS-gtopo30_ne11np4_16xconsistentSGH.c20160612.nc'
    if grid=='ne4np4'   : topo_file_name = f'{topo_file_root}/USGS-gtopo30_ne4pg2_16xdel2-PFC-consistentSGH.c20190618.nc'
    
    if topo_file_name is None: raise ValueError(file_err_msg)
    return topo_file_name
# ------------------------------------------------------------------------------
# Method for returning class object
# ------------------------------------------------------------------------------
def create_hiccup_data( src_data_name,
                        target_model=default_target_model,
                        dst_horz_grid=None,
                        dst_vert_grid=None,
                        atm_file=None,
                        sfc_file=None,
                        output_dir=default_output_dir,
                        grid_dir=default_grid_dir,
                        map_dir=default_map_dir,
                        tmp_dir=default_tmp_dir,
                        sstice_combined_file=None,
                        sstice_name=None,
                        sst_file=None,
                        ice_file=None,
                        topo_file=None,
                        lev_type=None,
                        verbose=False,
                        check_input_files=True,
                        RRM_grid=False,
                      ):
    """ 
    Create HICCUP data class object, check for required input variables and 
    create specified output directories if they do not exist
    """
    if verbose is not None:
        global hiccup_verbose
        hiccup_verbose = verbose
        hdc.hiccup_verbose = hiccup_verbose
        hu.hiccup_verbose = hiccup_verbose
    hu.check_nco_version()
    for subclass in hdc.hiccup_data.__subclasses__():
        if subclass.is_name_for(src_data_name):
            # Create the object
            obj = subclass( src_data_name,
                            target_model=target_model,
                            atm_file=atm_file,
                            sfc_file=sfc_file,
                            sstice_name=sstice_name,
                            sst_file=sst_file,
                            ice_file=ice_file,
                            sstice_combined_file=sstice_combined_file,
                            topo_file=topo_file,
                            dst_horz_grid=dst_horz_grid,
                            dst_vert_grid=dst_vert_grid,
                            output_dir=output_dir,
                            grid_dir=grid_dir,
                            map_dir=map_dir,
                            tmp_dir=tmp_dir,
                            lev_type=lev_type,
                            check_input_files=check_input_files,
                            RRM_grid=RRM_grid,
                          )

            # Check input files for for required variables
            if check_input_files: obj.check_file_vars()

            # Create the output, grid, and map folders if they do not exist
            for d in [output_dir,grid_dir,map_dir]:
                if d is not None: 
                    if not os.path.exists(d): os.makedirs(d)

            # global timer_start_total
            obj.timer_start_total = perf_counter()

            # Return the object if everything checks out
            return obj
    raise ValueError(f'{src_data_name} is not a supported HICCUP source dataset')
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------