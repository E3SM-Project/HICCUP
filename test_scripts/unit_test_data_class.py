#!/usr/bin/env python
#===================================================================================================
# Unit testing for hiccup_data_class module
#===================================================================================================
import unittest, numpy as np, xarray as xr, pandas as pd
from time import perf_counter
from hiccup.hiccup_data_class_timer_methods import print_timer
from hiccup.hiccup_utilities import compare_version, parse_version
from hiccup import hiccup

verbose_default = False # local verbosity default

#===============================================================================
class hiccup_data_class_test_case(unittest.TestCase):
  """ 
  Tests for hiccup_data_class_sstice_methods.py 
  """
  # ----------------------------------------------------------------------------
  def setUp(self):
    self.hiccup_data_ERA5 = hiccup.create_hiccup_data( src_data_name='ERA5',
                                                       input_file_list=['test_data/HICCUP_TEST.ERA5.atm.low-res.nc',
                                                                        'test_data/HICCUP_TEST.ERA5.sfc.low-res.nc'])
    self.hiccup_data_NOAA = hiccup.create_hiccup_data( src_data_name='NOAA',
                                                       sst_file='test_data/HICCUP_TEST.NOAA.sst.nc',
                                                       ice_file='test_data/HICCUP_TEST.NOAA.ice.nc')
    self.hiccup_data_CAMS = hiccup.create_hiccup_data( src_data_name='CAMS',
                                                       input_file_list=['test_data/HICCUP_TEST.CAMS.atm.low-res.nc',
                                                                        'test_data/HICCUP_TEST.CAMS.sfc.low-res.nc'])
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
    attr_list.append('output_dir')
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
#===============================================================================
if __name__ == '__main__':
    unittest.main()