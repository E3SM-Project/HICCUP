#!/usr/bin/env python
#===================================================================================================
# Unit testing for hiccup_data_class module
#===================================================================================================
import os, tempfile, unittest
import numpy as np, xarray as xr, pandas as pd
from time import perf_counter
from unittest.mock import patch
from hiccup.hiccup_data_class_timer_methods import print_timer
from hiccup.hiccup_data_class import _drop_ps, get_adj_file_list
from hiccup.hiccup_utilities import compare_version, parse_version
from hiccup import hiccup

verbose_default = False # local verbosity default

# Resolve test_data path relative to this file so tests run from any directory
TEST_DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'test_data')

#===============================================================================
class hiccup_data_class_test_case(unittest.TestCase):
  """ 
  Tests for hiccup_data_class_sstice_methods.py 
  """
  # ----------------------------------------------------------------------------
  def setUp(self):
    self.hiccup_data_ERA5 = hiccup.create_hiccup_data( src_data_name='ERA5',
                                                       input_file_list=[f'{TEST_DATA}/HICCUP_TEST.ERA5.atm.low-res.nc',
                                                                        f'{TEST_DATA}/HICCUP_TEST.ERA5.sfc.low-res.nc'])
    self.hiccup_data_NOAA = hiccup.create_hiccup_data( src_data_name='NOAA',
                                                       sst_file=f'{TEST_DATA}/HICCUP_TEST.NOAA.sst.nc',
                                                       ice_file=f'{TEST_DATA}/HICCUP_TEST.NOAA.ice.nc')
    self.hiccup_data_CAMS = hiccup.create_hiccup_data( src_data_name='CAMS',
                                                       input_file_list=[f'{TEST_DATA}/HICCUP_TEST.CAMS.atm.low-res.nc',
                                                                        f'{TEST_DATA}/HICCUP_TEST.CAMS.sfc.low-res.nc'])
    self.obj_list = []
    self.obj_list.append(self.hiccup_data_ERA5)
    self.obj_list.append(self.hiccup_data_NOAA)
    self.obj_list.append(self.hiccup_data_CAMS)
  # ----------------------------------------------------------------------------
  def test_hiccup_data_class_print(self):
    """ 
    test print methods of a hiccup data class for ERA5
    """
    timer_start = perf_counter()
    print_msg = self.hiccup_data_ERA5.__str__()
    print_msg = self.hiccup_data_NOAA.__str__()
    print_timer(timer_start,caller='test_hiccup_data_class_print')
  # ----------------------------------------------------------------------------
  def test_hiccup_data_class_atrributes(self):
    """ 
    check that various attributes are defined
    """
    timer_start = perf_counter()
    attr_list = []
    attr_list.append('src_data_name')
    attr_list.append('target_model')
    attr_list.append('input_file_list')
    attr_list.append('topo_file')
    attr_list.append('atm_var_name_dict')
    attr_list.append('sfc_var_name_dict')
    attr_list.append('src_horz_grid')
    attr_list.append('dst_horz_grid')
    attr_list.append('dst_vert_grid')
    attr_list.append('RRM_grid')
    attr_list.append('do_timers')
    attr_list.append('verbose')
    attr_list.append('verbose_indent')
    attr_list.append('grid_dir')
    attr_list.append('map_dir')
    attr_list.append('tmp_dir')
    attr_list.append('sstice_name')
    for obj in self.obj_list:
      for attr in attr_list:
        if not hasattr(obj, attr):
          raise AttributeError(f'{attr} attribute not found for hiccup data class with src_data_name: {obj.src_data_name}')
    print_timer(timer_start,caller='test_hiccup_data_class_atrributes')
  # ----------------------------------------------------------------------------
  def test_hiccup_data_class_get(self):
    """ 
    test "get" methods of the hiccup data class
    """
    timer_start = perf_counter()
    for obj in self.obj_list:
      self.assertIsNone( obj.get_src_grid_ne() )
      self.assertIsNone( obj.get_src_grid_npg() )
      self.assertIsNone( obj.get_dst_grid_ne() )
      self.assertIsNone( obj.get_dst_grid_npg() )
    print_timer(timer_start,caller='test_hiccup_data_class_get')
  # ----------------------------------------------------------------------------
  def test_hiccup_data_class_timer(self):
    """ 
    test that class timer attributes are properly adjusted
    """
    timer_start = perf_counter()
    for obj in self.obj_list:
      self.assertTrue(isinstance(obj.timer_start_total, (int, float, complex)))
      self.assertTrue(len(obj.timer_msg_all)==0)
      timer_start = perf_counter(); obj.print_timer(timer_start, use_color=False, print_msg=False)
      timer_start = perf_counter(); obj.print_timer(timer_start, use_color=False, print_msg=False)
      self.assertTrue(len(obj.timer_msg_all)==2)
    print_timer(timer_start,caller='test_hiccup_data_class_timer')
  # ----------------------------------------------------------------------------
  def test_hiccup_data_class_check_file_vars(self):
    """ 
    test that check_file_vars method doesn't cause any problems
    """
    timer_start = perf_counter()
    for obj in self.obj_list:
      obj.check_file_vars()
    print_timer(timer_start,caller='test_hiccup_data_class_check_file_vars')
  # ----------------------------------------------------------------------------
  def test_hiccup_data_class_get_multifile_dict(self):
    """ 
    test that get_multifile_dict method doesn't cause any problems
    """
    timer_start = perf_counter()
    for obj in self.obj_list:
      file_dict = obj.get_multifile_dict(timestamp=999)
      num_vars = len(obj.atm_var_name_dict) + len(obj.sfc_var_name_dict)
      if num_vars>0: num_vars = num_vars-2 # subtract 2 for lat/lon coords
      self.assertTrue(len(file_dict)==num_vars)
    print_timer(timer_start,caller='test_hiccup_data_class_get_multifile_dict')
  # ----------------------------------------------------------------------------
  def test_hiccup_data_class_get_dst_grid_ne(self):
    """
    test that get_dst_grid_ne correctly parses ne value from grid name strings
    """
    timer_start = perf_counter()
    obj = self.hiccup_data_ERA5
    for grid_name, expected_ne in [('ne30np4','30'), ('ne120np4','120'), ('ne30pg2','30'), ('ne256pg2','256')]:
      obj.dst_horz_grid = grid_name
      self.assertEqual(obj.get_dst_grid_ne(), expected_ne,
                       msg=f'get_dst_grid_ne failed for dst_horz_grid={grid_name!r}')
    print_timer(timer_start,caller='test_hiccup_data_class_get_dst_grid_ne')
  # ----------------------------------------------------------------------------
  def test_hiccup_data_class_get_dst_grid_npg(self):
    """
    test that get_dst_grid_npg correctly parses npg value from grid name strings
    """
    timer_start = perf_counter()
    obj = self.hiccup_data_ERA5
    # When dst_horz_grid_pg is explicitly set, it should take precedence
    for pg_str, expected_npg in [('pg2','2'), ('pg3','3')]:
      obj.dst_horz_grid_pg = pg_str
      self.assertEqual(obj.get_dst_grid_npg(), expected_npg,
                       msg=f'get_dst_grid_npg failed for dst_horz_grid_pg={pg_str!r}')
    # When dst_horz_grid_pg is absent, fall back to dst_horz_grid parsing
    del obj.dst_horz_grid_pg
    obj.dst_horz_grid = 'ne30pg2'
    self.assertEqual(obj.get_dst_grid_npg(), '2')
    obj.dst_horz_grid = 'ne30np4'
    self.assertEqual(obj.get_dst_grid_npg(), 0)  # no 'pg' in string -> returns 0
    print_timer(timer_start,caller='test_hiccup_data_class_get_dst_grid_npg')
  # ----------------------------------------------------------------------------
  def test_hiccup_data_class_get_dst_grid_ncol(self):
    """
    test that get_dst_grid_ncol returns correct column count for np4 and pg grids
    """
    timer_start = perf_counter()
    obj = self.hiccup_data_ERA5
    del obj.dst_horz_grid_pg  # use dst_horz_grid as fallback for npg parsing
    # np4 grid: ncol = ne^2 * 6 * 9 + 2
    obj.dst_horz_grid = 'ne30np4'
    self.assertEqual(obj.get_dst_grid_ncol(), 30*30*6*9+2)
    # pg2 grid: ncol = ne^2 * 6 * npg
    obj.dst_horz_grid = 'ne30pg2'
    self.assertEqual(obj.get_dst_grid_ncol(), 30*30*6*2)
    print_timer(timer_start,caller='test_hiccup_data_class_get_dst_grid_ncol')
  # ----------------------------------------------------------------------------
  def test_hiccup_data_class_get_chunks(self):
    """
    test that get_chunks returns a dict with the expected dimension keys
    """
    timer_start = perf_counter()
    chunks_ncol = self.hiccup_data_ERA5.get_chunks(ncol_only=True)
    self.assertIsInstance(chunks_ncol, dict)
    self.assertIn('ncol', chunks_ncol)
    chunks_all = self.hiccup_data_ERA5.get_chunks(ncol_only=False)
    self.assertIsInstance(chunks_all, dict)
    self.assertIn('ncol', chunks_all)
    self.assertIn('lev', chunks_all)
    print_timer(timer_start,caller='test_hiccup_data_class_get_chunks')
  # ----------------------------------------------------------------------------
  def test_hiccup_data_class_check_file_vars_missing(self):
    """
    test that check_file_vars_impl raises ValueError when a required variable is absent
    """
    timer_start = perf_counter()
    ds = xr.Dataset({'existing_var': xr.DataArray([1, 2, 3])})
    var_name_dict = {'MODEL_VAR': 'missing_src_var'}
    self.assertRaises(ValueError, self.hiccup_data_ERA5.check_file_vars_impl,
                      ds, var_name_dict, 'fake_file.nc')
    print_timer(timer_start,caller='test_hiccup_data_class_check_file_vars_missing')
  # ----------------------------------------------------------------------------
  def test_hiccup_data_class_add_time_date_variables(self):
    """
    test that add_time_date_variables adds the expected time/date variables
    """
    timer_start = perf_counter()
    time_values = pd.date_range('2020-06-15 12:00', periods=1)
    ds = xr.Dataset({'T': xr.DataArray(np.zeros((1,4)), dims=['time','ncol'])},
                    coords={'time': time_values})
    self.hiccup_data_ERA5.add_time_date_variables(ds)
    expected_vars = ['date','datesec','ndcur','nscur','nsteph',
                     'nbdate','ndbase','nsbase','nbsec','time_bnds']
    for var in expected_vars:
      self.assertIn(var, ds.variables, msg=f'{var!r} not found in dataset after add_time_date_variables')
    # date should encode YYYYMMDD
    self.assertEqual(int(ds['date'].values[0]), 20200615)
    # datesec should reflect the seconds-of-day (12:00 = 43200 seconds)
    self.assertEqual(int(ds['datesec'].values[0]), 0)  # seconds component only
    print_timer(timer_start,caller='test_hiccup_data_class_add_time_date_variables')
  # ----------------------------------------------------------------------------
  def test_hiccup_data_class_get_multifile_dict_eam(self):
    """
    test that get_multifile_dict_eam returns one entry per variable with .nc paths
    """
    timer_start = perf_counter()
    var_dict = {'T': 't', 'Q': 'q', 'U': 'u', 'V': 'v'}
    file_dict = self.hiccup_data_ERA5.get_multifile_dict_eam(var_dict, timestamp=999)
    self.assertEqual(len(file_dict), len(var_dict))
    for key in var_dict:
      self.assertIn(key, file_dict)
      self.assertTrue(file_dict[key].endswith('.nc'),
                      msg=f'Expected .nc path for key {key!r}, got {file_dict[key]!r}')
    print_timer(timer_start,caller='test_hiccup_data_class_get_multifile_dict_eam')
  # ----------------------------------------------------------------------------
  def test_hiccup_data_class_stage_multifile(self):
    """
    verify stage_multifile() copies each variable into the multifile layout
    (same file_dict paths as remap_horizontal_multifile()) without performing
    any horizontal remap
    """
    timer_start = perf_counter()
    obj = self.hiccup_data_ERA5
    file_dict = obj.get_multifile_dict(timestamp='stagetest')
    try:
      obj.stage_multifile(file_dict=file_dict)
      all_var_dicts = {**obj.atm_var_name_dict, **obj.sfc_var_name_dict}
      lat_var = obj.atm_var_name_dict['lat']
      lon_var = obj.atm_var_name_dict['lon']
      for key, path in file_dict.items():
        self.assertTrue(os.path.isfile(path), msg=f'stage_multifile did not create {path}')
        src_var = all_var_dicts[key]
        with xr.open_dataset(path) as ds_out:
          self.assertIn(src_var, ds_out.variables,
                        msg=f'{src_var} missing from staged file for {key!r}')
          self.assertIn(lat_var, ds_out.variables, msg='lat missing from staged file')
          self.assertIn(lon_var, ds_out.variables, msg='lon missing from staged file')
      # spot check that the data itself was copied through unchanged
      src_path = obj._var_to_file_map['T']
      with xr.open_dataset(src_path) as ds_src, xr.open_dataset(file_dict['T']) as ds_T:
        np.testing.assert_allclose(ds_T['t'].values, ds_src['t'].values)
    finally:
      for path in file_dict.values():
        if os.path.isfile(path): os.remove(path)
    print_timer(timer_start,caller='test_hiccup_data_class_stage_multifile')
  # ----------------------------------------------------------------------------
  def test_hiccup_data_class_stage_multifile_eam(self):
    """
    verify stage_multifile_eam() copies each variable into the multifile layout
    (same file_dict paths as remap_horizontal_multifile_eam()) without performing
    any horizontal remap, using identity-named ("EAM-style") source variables
    """
    timer_start = perf_counter()
    obj = self.hiccup_data_ERA5   # reuse a real object; stage_multifile_eam is generic

    # build a small synthetic "EAM-style" source file (identity-named variables,
    # unlike ERA5 test data where e.g. 'T' is stored as 't')
    ncol = 4
    ds_in = xr.Dataset({
      'T':   (('ncol',), np.arange(ncol, dtype='float64')),
      'PS':  (('ncol',), np.arange(ncol, dtype='float64')*10.),
      'lat': (('ncol',), np.linspace(-1,1,ncol)),
      'lon': (('ncol',), np.linspace(0,2,ncol)),
    })
    tmp_in = tempfile.NamedTemporaryFile(suffix='.nc', delete=False)
    tmp_in.close()
    # disable xarray's default NaN _FillValue so check_file_FillValue() doesn't
    # decide to write a "modified" copy and rebuild _var_to_file_map on us
    no_fill_encoding = {v: {'_FillValue': None} for v in ds_in.data_vars}
    ds_in.to_netcdf(tmp_in.name, encoding=no_fill_encoding)

    file_dict = {'T':  os.path.join(obj.tmp_dir,'stage_eam_test.T.nc'),
                 'PS': os.path.join(obj.tmp_dir,'stage_eam_test.PS.nc')}

    orig_input_file_list = obj.input_file_list
    orig_var_to_file_map = obj._var_to_file_map
    try:
      obj.input_file_list  = [tmp_in.name]
      obj._var_to_file_map = {'T': tmp_in.name, 'PS': tmp_in.name}

      obj.stage_multifile_eam(file_dict=file_dict)

      for key, path in file_dict.items():
        self.assertTrue(os.path.isfile(path), msg=f'stage_multifile_eam did not create {path}')
        with xr.open_dataset(path) as ds_out:
          self.assertIn(key, ds_out.variables, msg=f'{key} missing from staged file')
          self.assertIn('lat', ds_out.variables, msg='lat missing from staged file')
          self.assertIn('lon', ds_out.variables, msg='lon missing from staged file')
          np.testing.assert_allclose(ds_out[key].values, ds_in[key].values)
    finally:
      obj.input_file_list  = orig_input_file_list
      obj._var_to_file_map = orig_var_to_file_map
      os.remove(tmp_in.name)
      for path in file_dict.values():
        if os.path.isfile(path): os.remove(path)

    print_timer(timer_start,caller='test_hiccup_data_class_stage_multifile_eam')
  # ----------------------------------------------------------------------------
  def test_parse_version(self):
    """
    test that parse_version correctly parses version strings with and without suffixes
    """
    timer_start = perf_counter()
    main, suffix = parse_version('5.3.1')
    self.assertEqual(main, (5, 3, 1))
    main, suffix = parse_version('5.3.1-alpha09')
    self.assertEqual(main, (5, 3, 1))
    self.assertEqual(suffix[1], 9)  # alpha suffix with numeric part 9
    print_timer(timer_start,caller='test_parse_version')
  # ----------------------------------------------------------------------------
  def test_compare_version(self):
    """
    test that compare_version correctly identifies whether version >= required_version
    """
    timer_start = perf_counter()
    self.assertTrue( compare_version('5.3.1', required_version='5.3.1'))   # equal
    self.assertTrue( compare_version('5.4.0', required_version='5.3.1'))   # newer minor
    self.assertTrue( compare_version('6.0.0', required_version='5.3.1'))   # newer major
    self.assertFalse(compare_version('5.2.0', required_version='5.3.1'))   # older minor
    self.assertFalse(compare_version('4.9.9', required_version='5.3.1'))   # older major
    self.assertFalse(compare_version('5.3.1-alpha09', required_version='5.3.1'))  # alpha < release
    print_timer(timer_start,caller='test_compare_version')
  # ----------------------------------------------------------------------------
  def test_build_var_to_file_map_prefers_2d_for_sfc_vars(self):
    """
    Verify that _build_var_to_file_map routes sfc variables to a file where the
    source variable has no vertical dimension, even when a same-named 3-D variable
    exists in an earlier file.

    ERA5 test data: both atm and sfc files contain 'z', but the atm copy has a
    pressure_level dimension while the sfc copy is 2-D.  PHIS/phis must map to
    the sfc file so that the surface-adjustment step does not pick up 3-D data.
    """
    timer_start = perf_counter()
    import xarray as xr, os

    obj = self.hiccup_data_ERA5   # ERA5/EAM object from setUp

    # Confirm 'z' is present in both files
    atm_file = obj.input_file_list[0]
    sfc_file = obj.input_file_list[1]
    ds_atm = xr.open_dataset(atm_file, decode_times=False)
    ds_sfc = xr.open_dataset(sfc_file, decode_times=False)
    self.assertIn('z', ds_atm, msg='test precondition: z must exist in atm file')
    self.assertIn('z', ds_sfc, msg='test precondition: z must exist in sfc file')
    self.assertTrue(
      any(d in {'pressure_level','level','plev','lev'} for d in ds_atm['z'].dims),
      msg='test precondition: z in atm file must have a vertical dimension')
    self.assertFalse(
      any(d in {'pressure_level','level','plev','lev'} for d in ds_sfc['z'].dims),
      msg='test precondition: z in sfc file must NOT have a vertical dimension')

    # The map must route PHIS (src='z') to the sfc file, not the atm file
    self.assertIn('PHIS', obj._var_to_file_map,
                  msg='PHIS must be present in _var_to_file_map')
    self.assertEqual(
      os.path.abspath(obj._var_to_file_map['PHIS']),
      os.path.abspath(sfc_file),
      msg='PHIS must map to sfc file (2-D z), not atm file (3-D z on pressure levels)')

    # Atmospheric variable T must still map to the atm file
    self.assertIn('T', obj._var_to_file_map,
                  msg='T must be present in _var_to_file_map')
    self.assertEqual(
      os.path.abspath(obj._var_to_file_map['T']),
      os.path.abspath(atm_file),
      msg='T must map to atm file')

    print_timer(timer_start,caller='test_build_var_to_file_map_prefers_2d_for_sfc_vars')
  # ----------------------------------------------------------------------------
  def test_hiccup_data_class_str_content(self):
    """
    verify __str__ output contains key attribute names and object header
    """
    timer_start = perf_counter()
    s = self.hiccup_data_ERA5.__str__()
    self.assertIn('HICCUP data object', s)
    self.assertIn('src_data_name', s)
    self.assertIn('ERA5', s)
    self.assertIn('atm_var_name_dict', s)
    self.assertIn('sfc_var_name_dict', s)
    print_timer(timer_start, caller='test_hiccup_data_class_str_content')
  # ----------------------------------------------------------------------------
  def test_hiccup_data_class_init_error_target_model_none(self):
    """
    verify that target_model=None raises ValueError
    """
    timer_start = perf_counter()
    self.assertRaises(
      ValueError,
      hiccup.create_hiccup_data,
      'ERA5',
      target_model=None,
      input_file_list=[f'{TEST_DATA}/HICCUP_TEST.ERA5.atm.low-res.nc',
                       f'{TEST_DATA}/HICCUP_TEST.ERA5.sfc.low-res.nc'])
    print_timer(timer_start, caller='test_hiccup_data_class_init_error_target_model_none')
  # ----------------------------------------------------------------------------
  def test_hiccup_data_class_init_error_missing_file(self):
    """
    verify that a non-existent input file raises ValueError when check_input_files=True
    """
    timer_start = perf_counter()
    self.assertRaises(
      ValueError,
      hiccup.create_hiccup_data,
      'ERA5',
      input_file_list=[f'{TEST_DATA}/HICCUP_TEST.ERA5.atm.low-res.nc',
                       '/nonexistent/file.nc'],
      check_input_files=True)
    print_timer(timer_start, caller='test_hiccup_data_class_init_error_missing_file')
  # ----------------------------------------------------------------------------
  def test_hiccup_data_class_init_error_no_input_files(self):
    """
    verify that an empty input_file_list raises ValueError for ERA5 (which has var dicts)
    """
    timer_start = perf_counter()
    self.assertRaises(
      ValueError,
      hiccup.create_hiccup_data,
      'ERA5',
      input_file_list=[])
    print_timer(timer_start, caller='test_hiccup_data_class_init_error_no_input_files')
  # ----------------------------------------------------------------------------
  def test_hiccup_data_class_check_file_vars_populates_map(self):
    """
    verify _var_to_file_map is populated after check_file_vars()
    """
    timer_start = perf_counter()
    obj = self.hiccup_data_ERA5
    obj.check_file_vars()
    self.assertIsInstance(obj._var_to_file_map, dict)
    self.assertGreater(len(obj._var_to_file_map), 0)
    # every key in the combined var dicts (except lat/lon) should be mapped
    all_vars = {**obj.atm_var_name_dict, **obj.sfc_var_name_dict}
    for var in all_vars:
      if var not in ('lat', 'lon'):
        self.assertIn(var, obj._var_to_file_map,
                      msg=f'Expected {var!r} in _var_to_file_map')
    print_timer(timer_start, caller='test_hiccup_data_class_check_file_vars_populates_map')
  # ----------------------------------------------------------------------------
  def test_hiccup_data_class_get_multifile_dict_prefix_logic(self):
    """
    verify atm vars get tmp_atm_data prefix, sfc vars get tmp_sfc_data prefix,
    and lat/lon coords are excluded from the dict
    """
    timer_start = perf_counter()
    obj = self.hiccup_data_ERA5
    file_dict = obj.get_multifile_dict(timestamp=999)
    lat_var = obj.atm_var_name_dict.get('lat')
    lon_var = obj.atm_var_name_dict.get('lon')
    # lat and lon should never appear as keys
    self.assertNotIn('lat', file_dict)
    self.assertNotIn('lon', file_dict)
    for key, path in file_dict.items():
      if key in obj.sfc_var_name_dict and key not in obj.atm_var_name_dict:
        self.assertIn('tmp_sfc_data', path,
                      msg=f'sfc var {key!r} expected tmp_sfc_data prefix, got {path!r}')
      elif key in obj.atm_var_name_dict and obj.atm_var_name_dict[key] not in [lat_var, lon_var]:
        self.assertIn('tmp_atm_data', path,
                      msg=f'atm var {key!r} expected tmp_atm_data prefix, got {path!r}')
    print_timer(timer_start, caller='test_hiccup_data_class_get_multifile_dict_prefix_logic')
  # ----------------------------------------------------------------------------
  def test_drop_ps_helper(self):
    """
    verify _drop_ps removes PS from non-PS files and keeps it in the PS file
    """
    timer_start = perf_counter()
    ps_path = '/fake/ps.nc'
    t_path  = '/fake/t.nc'
    file_dict = {'PS': ps_path, 'T': t_path}
    # Dataset sourced from PS file — PS must be kept
    ds_ps = xr.Dataset({'PS': xr.DataArray([1.0]), 'T': xr.DataArray([2.0])})
    ds_ps.encoding['source'] = ps_path
    # Dataset sourced from T file — PS must be dropped
    ds_t = xr.Dataset({'PS': xr.DataArray([1.0]), 'T': xr.DataArray([2.0])})
    ds_t.encoding['source'] = t_path
    with patch('hiccup.hiccup_data_class.os.path.samefile', side_effect=lambda a, b: a == b):
      result_ps = _drop_ps(ds_ps, file_dict)
      result_t  = _drop_ps(ds_t,  file_dict)
    self.assertIn('PS', result_ps, msg='PS must be kept in the PS file dataset')
    self.assertNotIn('PS', result_t, msg='PS must be dropped from the T file dataset')
    print_timer(timer_start, caller='test_drop_ps_helper')
  # ----------------------------------------------------------------------------
  def test_get_adj_file_list_helper(self):
    """
    verify get_adj_file_list returns only the files whose keys appear in var_list
    """
    timer_start = perf_counter()
    file_dict = {'T': '/tmp/t.nc', 'Q': '/tmp/q.nc', 'PS': '/tmp/ps.nc', 'U': '/tmp/u.nc'}
    result = get_adj_file_list(['Q', 'PS'], file_dict)
    self.assertEqual(result, ['/tmp/q.nc', '/tmp/ps.nc'])
    # empty var_list returns empty list
    self.assertEqual(get_adj_file_list([], file_dict), [])
    print_timer(timer_start, caller='test_get_adj_file_list_helper')
  # ----------------------------------------------------------------------------
  def test_surface_adjustment_multifile_error_guard(self):
    """
    verify adj_T_eam=True with adj_PS=False raises ValueError
    """
    timer_start = perf_counter()
    obj = self.hiccup_data_ERA5
    self.assertRaises(ValueError, obj.surface_adjustment_multifile,
                      {}, adj_T_eam=True, adj_PS=False)
    print_timer(timer_start, caller='test_surface_adjustment_multifile_error_guard')
  # ----------------------------------------------------------------------------
  def test_atmos_state_adjustment_multifile_adjust_wtr_empty(self):
    """
    verify early return when no cloud water files are in file_dict
    """
    timer_start = perf_counter()
    # file_dict with no CLDLIQ/CLDICE/qc/qi keys -> get_adj_file_list returns []
    obj = self.hiccup_data_ERA5
    obj.atmos_state_adjustment_multifile_adjust_wtr(file_dict={})  # should not raise
    print_timer(timer_start, caller='test_atmos_state_adjustment_multifile_adjust_wtr_empty')
  # ----------------------------------------------------------------------------
  def test_check_input_times_and_get_time_index(self):
    """
    verify check_input_times passes for a present time, raises for an absent one,
    and _get_time_index returns the correct integer index
    """
    timer_start = perf_counter()
    # build a tiny temp netCDF with a proper 'time' dimension
    tmp = tempfile.NamedTemporaryFile(suffix='.nc', delete=False)
    tmp.close()
    try:
      time_vals = pd.date_range('2008-10-01', periods=3, freq='6h')
      ds = xr.Dataset({'T': xr.DataArray(np.zeros((3, 4)), dims=['time', 'ncol'])},
                      coords={'time': time_vals})
      ds.to_netcdf(tmp.name)
      ds.close()
      obj = self.hiccup_data_ERA5
      obj.input_file_list = [tmp.name]
      # happy path — time is in the file
      obj.check_input_times('2008-10-01')
      # error path — time is absent
      self.assertRaises(ValueError, obj.check_input_times, '2020-01-01')
      # _get_time_index returns correct indices
      self.assertEqual(obj._get_time_index(tmp.name, '2008-10-01 00:00'), 0)
      self.assertEqual(obj._get_time_index(tmp.name, '2008-10-01 06:00'), 1)
      self.assertEqual(obj._get_time_index(tmp.name, '2008-10-01 12:00'), 2)
    finally:
      os.remove(tmp.name)
    print_timer(timer_start, caller='test_check_input_times_and_get_time_index')
  # ----------------------------------------------------------------------------
  def _make_dim_order_ds(self, dim_order=('time','lev','ncol')):
    """helper - tiny dataset whose 3D fields use the requested dim order"""
    sizes = {'time':1,'lev':4,'ncol':6}
    shape = tuple(sizes[d] for d in dim_order)
    return xr.Dataset({
      'T':  xr.DataArray(np.zeros(shape), dims=list(dim_order)),
      'PS': xr.DataArray(np.zeros((1,6)), dims=['time','ncol']),
    })
  # ----------------------------------------------------------------------------
  def test_check_dimension_order(self):
    """
    check_dimension_order should accept the expected layout and reject a transposed one
    """
    timer_start = perf_counter()
    obj = self.hiccup_data_ERA5

    # happy path - EAM layout
    ds_good = self._make_dim_order_ds(('time','lev','ncol'))
    obj.check_dimension_order(ds_good,['time','lev','ncol'],verbose=False)

    # error path - transposed field is the regression this check exists to catch
    ds_bad = self._make_dim_order_ds(('time','ncol','lev'))
    with self.assertRaises(ValueError) as ctx:
      obj.check_dimension_order(ds_bad,['time','lev','ncol'],verbose=False)
    self.assertIn('T',str(ctx.exception))

    # the same dataset is valid under the EAMxx expected order
    obj.check_dimension_order(ds_bad,['time','ncol','lev'],verbose=False)
    print_timer(timer_start, caller='test_check_dimension_order')
  # ----------------------------------------------------------------------------
  def test_check_dimension_order_ignores_extra_dims(self):
    """
    only the relative order of recognized dims matters - extra dims (dim2, nv, nbnd)
    may appear anywhere, so EAMxx horiz_winds and grid metadata must not trip the check
    """
    timer_start = perf_counter()
    obj = self.hiccup_data_ERA5
    ds = xr.Dataset({
      'horiz_winds':   xr.DataArray(np.zeros((1,6,2,4)), dims=['time','ncol','dim2','lev']),
      'lat_vertices':  xr.DataArray(np.zeros((6,4)),     dims=['ncol','nv']),
      'time_bnds':     xr.DataArray(np.zeros((1,2)),     dims=['time','nbnd']),
      'hyai':          xr.DataArray(np.zeros(5),         dims=['ilev']),
    })
    obj.check_dimension_order(ds,['time','ncol','lev'],verbose=False)
    print_timer(timer_start, caller='test_check_dimension_order_ignores_extra_dims')
  # ----------------------------------------------------------------------------
  def test_check_dimension_order_handles_ncol_d(self):
    """
    a dataset using ncol_d should be checked against ncol_d in the expected slot
    """
    timer_start = perf_counter()
    obj = self.hiccup_data_ERA5
    ds_good = xr.Dataset({'T': xr.DataArray(np.zeros((1,4,6)), dims=['time','lev','ncol_d'])})
    obj.check_dimension_order(ds_good,['time','lev','ncol'],verbose=False)
    ds_bad = xr.Dataset({'T': xr.DataArray(np.zeros((1,6,4)), dims=['time','ncol_d','lev'])})
    self.assertRaises(ValueError, obj.check_dimension_order, ds_bad, ['time','lev','ncol'], None, False)
    print_timer(timer_start, caller='test_check_dimension_order_handles_ncol_d')
  # ----------------------------------------------------------------------------
  def test_verify_target_model_valid_names(self):
    """
    all supported target models should be accepted and returned in canonical spelling
    """
    timer_start = perf_counter()
    for name in ['EAM','EAMXX','EAMXX-nudging']:
      self.assertEqual(hiccup.verify_target_model(name), name)
    print_timer(timer_start, caller='test_verify_target_model_valid_names')
  # ----------------------------------------------------------------------------
  def test_verify_target_model_is_case_insensitive(self):
    """
    input is matched case-insensitively but normalized back to the canonical spelling,
    since the rest of HICCUP compares against these exact strings - "EAMXX-nudging" is
    mixed case, so returning an upper-cased name would silently disable its branches
    """
    timer_start = perf_counter()
    self.assertEqual(hiccup.verify_target_model('eam'), 'EAM')
    self.assertEqual(hiccup.verify_target_model('eamxx'), 'EAMXX')
    for name in ['eamxx-nudging','EAMXX-NUDGING','EAMxx-Nudging']:
      self.assertEqual(hiccup.verify_target_model(name), 'EAMXX-nudging')
    print_timer(timer_start, caller='test_verify_target_model_is_case_insensitive')
  # ----------------------------------------------------------------------------
  def test_verify_target_model_invalid(self):
    """
    unrecognized names and non-string input should raise ValueError
    """
    timer_start = perf_counter()
    for bad in ['EAMXX-foo','CAM','',None,5,['EAM']]:
      self.assertRaises(ValueError, hiccup.verify_target_model, bad)
    print_timer(timer_start, caller='test_verify_target_model_invalid')
  # ----------------------------------------------------------------------------
  def test_create_hiccup_data_accepts_eamxx_nudging(self):
    """
    EAMXX-nudging must survive create_hiccup_data() - it is a supported target model
    with dedicated branches in combine_files() and remap_vertical_multifile()
    """
    timer_start = perf_counter()
    obj = hiccup.create_hiccup_data( src_data_name='ERA5',
                                     target_model='EAMXX-nudging',
                                     dst_horz_grid='ne30np4',
                                     dst_vert_grid='L128',
                                     input_file_list=[f'{TEST_DATA}/HICCUP_TEST.ERA5.atm.low-res.nc',
                                                      f'{TEST_DATA}/HICCUP_TEST.ERA5.sfc.low-res.nc'])
    self.assertEqual(obj.target_model, 'EAMXX-nudging')
    print_timer(timer_start, caller='test_create_hiccup_data_accepts_eamxx_nudging')
#===============================================================================
if __name__ == '__main__':
    unittest.main()