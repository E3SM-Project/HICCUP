# ------------------------------------------------------------------------------
# Methods for populating grid information for a hiccup_data object
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
        ds_grid = xr.open_dataset(self.map_file)
        ncol = int(ds_grid['n_b'].values)
    else:
        ne  = int(self.get_dst_grid_ne())
        npg = int(self.get_dst_grid_npg())
        if npg==0: ncol = int(ne*ne*6*9+2)
        if npg>0 : ncol = int(ne*ne*6*npg)
    return ncol
# ------------------------------------------------------------------------------
