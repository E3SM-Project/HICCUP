# ------------------------------------------------------------------------------
# Methods for populating grid information for a hiccup_data object
# ------------------------------------------------------------------------------
from hiccup.hiccup_data_class_common import *
# ------------------------------------------------------------------------------
def get_src_grid_ne(self):
    """
    Return number of elements of source grid (if starting from model data)
    """
    if hasattr(self, 'src_horz_grid_np'):
        if self.src_horz_grid_np is None: return
        result = re.search('ne(.*)np', self.src_horz_grid_np)
        return result.group(1) if result else 0
    else:
        raise AttributeError('src_horz_grid_np attribute not found!')
# ------------------------------------------------------------------------------
def get_src_grid_npg(self):
    """
    Return number of FV physgrid cells (npg) of source grid (if starting from model data)
    """
    if hasattr(self, 'src_horz_grid_pg'):
        if self.src_horz_grid_pg is None: return
        result = re.search('pg(.*)', self.src_horz_grid_pg)
        return result.group(1) if result else 0
    else:
        raise AttributeError('src_horz_grid_pg attribute not found!')
# ------------------------------------------------------------------------------
def get_dst_grid_ne(self):
    """
    Return number of elements of target model grid
    """
    if hasattr(self, 'dst_horz_grid'):
        if self.dst_horz_grid is None: return
        if 'np4' in self.dst_horz_grid:
            result = re.search('ne(.*)np', self.dst_horz_grid)
        if 'pg' in self.dst_horz_grid:
            result = re.search('ne(.*)pg', self.dst_horz_grid)
        return result.group(1) if result else 0
    else:
        raise AttributeError('dst_horz_grid attribute not found!')
# ------------------------------------------------------------------------------
def get_dst_grid_npg(self):
    """
    Return number of FV physgrid cells (npg) of target model grid
    """
    if hasattr(self, 'dst_horz_grid_pg'):
        if self.dst_horz_grid_pg is None: return
        result = re.search('pg(.*)', self.dst_horz_grid_pg)
        return result.group(1) if result else 0
    elif hasattr(self, 'dst_horz_grid'):
        result = re.search('pg(.*)', self.dst_horz_grid)
        return result.group(1) if result else 0
    else:
        raise AttributeError('dst_horz_grid_pg and dst_horz_grid attributes not found!')
# ------------------------------------------------------------------------------
def get_dst_grid_ncol(self):
    """
    Return ncol for destination grid
    """
    if self.RRM_grid:
        # use map file to determine ncol
        if self.map_file is None : raise ValueError('get_dst_grid_ncol: ncol cannot be determined')
        with xr.open_dataset(self.map_file) as ds_grid:
            ncol = int(ds_grid['n_b'].values)
    else:
        ne  = int(self.get_dst_grid_ne())
        npg = int(self.get_dst_grid_npg())
        if npg==0: ncol = int(ne*ne*6*9+2)
        if npg>0 : ncol = int(ne*ne*6*npg)
    return ncol
# ------------------------------------------------------------------------------
def get_grid_cell_count(grid_file):
    """
    Return the number of cells/elements in a SCRIP or exodus grid file
    """
    import netCDF4
    # NOTE: use netCDF4 directly instead of xarray - some grid files (e.g. exodus)
    # declare dimensions (like num_elem) that are not attached to any variable,
    # so xarray's ds.sizes silently omits them even though ncdump -h shows them
    direct_dim_list   = ['grid_size','ncol','elementCount']
    spectral_dim_list = ['num_elem','num_el_in_blk1']
    with netCDF4.Dataset(grid_file) as ds_grid:
        for dim in direct_dim_list:
            if dim in ds_grid.dimensions: return int(ds_grid.dimensions[dim].size)
        for dim in spectral_dim_list:
            if dim in ds_grid.dimensions: return int(ds_grid.dimensions[dim].size*9+2)
    raise ValueError(f'get_grid_cell_count: could not determine the number of cells'
                     f' in {grid_file} - expected one of these dimensions: '
                     f'{direct_dim_list+spectral_dim_list}')
# ------------------------------------------------------------------------------
def check_lrg2sml(self,src_grid_file=None,dst_grid_file=None,verbose=None):
    """
    Determine whether the ncremap --lrg2sml flag is needed. This flag swaps the
    order of the grid arguments given to GenerateOverlapMesh, which requires the
    grid with the smaller cells to come first, so the flag is needed whenever the
    destination grid is finer than the source grid (i.e. "very fine" RRM grids).
    """
    if verbose is None: verbose = self.verbose
    if src_grid_file is None: src_grid_file = self.src_grid_file
    if dst_grid_file is None: dst_grid_file = self.dst_grid_file
    src_cell_cnt = get_grid_cell_count(src_grid_file)
    dst_cell_cnt = get_grid_cell_count(dst_grid_file)
    lrg2sml = dst_cell_cnt > src_cell_cnt
    if verbose:
        print(f'{self.verbose_indent}  src grid cells: {src_cell_cnt}'
              f' / dst grid cells: {dst_cell_cnt} => lrg2sml = {lrg2sml}')
    return lrg2sml
# ------------------------------------------------------------------------------
