#!/usr/bin/env python
#===================================================================================================
# Unit testing for hiccup_data_class module
#===================================================================================================
import unittest, numpy as np, xarray as xr
from time import perf_counter
from hiccup.hiccup_data_class_timer_methods import print_timer
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
                                                       atm_file='test_data/HICCUP_TEST.ERA5.atm.low-res.nc',
                                                       sfc_file='test_data/HICCUP_TEST.ERA5.sfc.low-res.nc')
    self.hiccup_data_NOAA = hiccup.create_hiccup_data( src_data_name='NOAA',
                                                       sst_file='test_data/HICCUP_TEST.NOAA.sst.nc',
                                                       ice_file='test_data/HICCUP_TEST.NOAA.ice.nc')
    self.obj_list = []
    self.obj_list.append(self.hiccup_data_ERA5)
    self.obj_list.append(self.hiccup_data_NOAA)
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
    attr_list.append('atm_file')
    attr_list.append('sfc_file')
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
#===============================================================================
if __name__ == '__main__':
    unittest.main()