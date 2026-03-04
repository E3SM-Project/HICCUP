#!/usr/bin/env python
#===================================================================================================
# Unit testing for hiccup_utilities module
#===================================================================================================
import unittest
from unittest.mock import patch
import numpy as np
import xarray as xr
from time import perf_counter
from hiccup.hiccup_data_class_timer_methods import print_timer
from hiccup.hiccup_utilities import (tcolor, run_cmd, check_dependency,
                                      suffix_as_tuple, parse_version,
                                      compare_version, check_nco_version,
                                      print_stat, chk_finite)

verbose_default = False # local verbosity default

#===============================================================================
class tcolor_test_case(unittest.TestCase):
  """
  Tests for tcolor class
  """
  # ----------------------------------------------------------------------------
  def test_tcolor_attributes_are_strings(self):
    """
    test that all tcolor attributes are strings
    """
    timer_start = perf_counter()
    for attr in ['ENDC','BLACK','RED','GREEN','YELLOW','BLUE','MAGENTA','CYAN','WHITE']:
      self.assertIsInstance(getattr(tcolor, attr), str)
    print_timer(timer_start, caller='test_tcolor_attributes_are_strings')
  # ----------------------------------------------------------------------------
  def test_tcolor_attributes_are_escape_codes(self):
    """
    test that all tcolor attributes start with ESC character
    """
    timer_start = perf_counter()
    for attr in ['ENDC','BLACK','RED','GREEN','YELLOW','BLUE','MAGENTA','CYAN','WHITE']:
      self.assertTrue(getattr(tcolor, attr).startswith('\033['))
    print_timer(timer_start, caller='test_tcolor_attributes_are_escape_codes')
  # ----------------------------------------------------------------------------
  def test_tcolor_endc_resets(self):
    """
    test that ENDC is the standard reset code
    """
    timer_start = perf_counter()
    self.assertEqual(tcolor.ENDC, '\033[0m')
    print_timer(timer_start, caller='test_tcolor_endc_resets')
  # ----------------------------------------------------------------------------

#===============================================================================
class run_cmd_test_case(unittest.TestCase):
  """
  Tests for run_cmd()
  """
  # ----------------------------------------------------------------------------
  @patch('hiccup.hiccup_utilities.sp.check_call')
  def test_run_cmd_calls_check_call(self, mock_check_call):
    """
    test that run_cmd invokes subprocess check_call
    """
    timer_start = perf_counter()
    run_cmd('echo hello', use_color=False)
    mock_check_call.assert_called_once()
    print_timer(timer_start, caller='test_run_cmd_calls_check_call')
  # ----------------------------------------------------------------------------
  @patch('hiccup.hiccup_utilities.sp.check_call')
  def test_run_cmd_splits_command_by_default(self, mock_check_call):
    """
    test that run_cmd splits command string when shell=False
    """
    timer_start = perf_counter()
    run_cmd('echo hello world', shell=False, use_color=False)
    args, kwargs = mock_check_call.call_args
    self.assertEqual(args[0], ['echo', 'hello', 'world'])
    print_timer(timer_start, caller='test_run_cmd_splits_command_by_default')
  # ----------------------------------------------------------------------------
  @patch('hiccup.hiccup_utilities.sp.check_call')
  def test_run_cmd_shell_mode(self, mock_check_call):
    """
    test that run_cmd passes shell=True when requested
    """
    timer_start = perf_counter()
    run_cmd('echo hello', shell=True, use_color=False)
    args, kwargs = mock_check_call.call_args
    self.assertTrue(kwargs.get('shell', False))
    print_timer(timer_start, caller='test_run_cmd_shell_mode')
  # ----------------------------------------------------------------------------
  @patch('builtins.print')
  @patch('hiccup.hiccup_utilities.sp.check_call')
  def test_run_cmd_verbose_prints(self, mock_check_call, mock_print):
    """
    test that run_cmd prints when verbose=True
    """
    timer_start = perf_counter()
    run_cmd('echo hello', verbose=True, use_color=False)
    mock_print.assert_called_once()
    print_timer(timer_start, caller='test_run_cmd_verbose_prints')
  # ----------------------------------------------------------------------------
  @patch('builtins.print')
  @patch('hiccup.hiccup_utilities.sp.check_call')
  def test_run_cmd_silent_by_default(self, mock_check_call, mock_print):
    """
    test that run_cmd does not print when verbose=False
    """
    timer_start = perf_counter()
    run_cmd('echo hello', verbose=False, use_color=False)
    mock_print.assert_not_called()
    print_timer(timer_start, caller='test_run_cmd_silent_by_default')
  # ----------------------------------------------------------------------------

#===============================================================================
class check_dependency_test_case(unittest.TestCase):
  """
  Tests for check_dependency()
  """
  # ----------------------------------------------------------------------------
  @patch('hiccup.hiccup_utilities.shutil.which', return_value='/usr/bin/ls')
  def test_check_dependency_found(self, _):
    """
    test that check_dependency does not raise when command is found
    """
    timer_start = perf_counter()
    check_dependency('ls')  # should not raise
    print_timer(timer_start, caller='test_check_dependency_found')
  # ----------------------------------------------------------------------------
  @patch('hiccup.hiccup_utilities.shutil.which', return_value=None)
  def test_check_dependency_not_found_raises(self, _):
    """
    test that check_dependency raises OSError when command is not found
    """
    timer_start = perf_counter()
    with self.assertRaises(OSError):
      check_dependency('nonexistent_tool_xyz')
    print_timer(timer_start, caller='test_check_dependency_not_found_raises')
  # ----------------------------------------------------------------------------
  @patch('hiccup.hiccup_utilities.shutil.which', return_value=None)
  def test_check_dependency_error_message_contains_cmd(self, _):
    """
    test that OSError message contains the missing command name
    """
    timer_start = perf_counter()
    with self.assertRaises(OSError) as ctx:
      check_dependency('my_missing_tool')
    self.assertIn('my_missing_tool', str(ctx.exception))
    print_timer(timer_start, caller='test_check_dependency_error_message_contains_cmd')
  # ----------------------------------------------------------------------------

#===============================================================================
class version_parsing_test_case(unittest.TestCase):
  """
  Tests for suffix_as_tuple(), parse_version(), and compare_version()
  """
  # ----------------------------------------------------------------------------
  def test_suffix_as_tuple_empty(self):
    """
    test that empty suffix returns tuple with highest order index
    """
    timer_start = perf_counter()
    result = suffix_as_tuple('')
    self.assertIsInstance(result, tuple)
    self.assertEqual(len(result), 2)
    print_timer(timer_start, caller='test_suffix_as_tuple_empty')
  # ----------------------------------------------------------------------------
  def test_suffix_as_tuple_alpha(self):
    """
    test that alpha suffix has lower order than beta and release
    """
    timer_start = perf_counter()
    alpha = suffix_as_tuple('alpha09')
    beta  = suffix_as_tuple('beta01')
    rel   = suffix_as_tuple('')
    self.assertLess(alpha, beta)
    self.assertLess(beta, rel)
    print_timer(timer_start, caller='test_suffix_as_tuple_alpha')
  # ----------------------------------------------------------------------------
  def test_suffix_as_tuple_invalid_raises(self):
    """
    test that an unrecognised suffix raises AssertionError
    """
    timer_start = perf_counter()
    with self.assertRaises(AssertionError):
      suffix_as_tuple('rc1')
    print_timer(timer_start, caller='test_suffix_as_tuple_invalid_raises')
  # ----------------------------------------------------------------------------
  def test_parse_version_simple(self):
    """
    test parse_version with a plain version string
    """
    timer_start = perf_counter()
    main, suffix = parse_version('5.3.1')
    self.assertEqual(main, (5, 3, 1))
    print_timer(timer_start, caller='test_parse_version_simple')
  # ----------------------------------------------------------------------------
  def test_parse_version_with_alpha_suffix(self):
    """
    test parse_version with an alpha suffix
    """
    timer_start = perf_counter()
    main, suffix = parse_version('4.9.2-alpha09')
    self.assertEqual(main, (4, 9, 2))
    print_timer(timer_start, caller='test_parse_version_with_alpha_suffix')
  # ----------------------------------------------------------------------------
  def test_compare_version_equal(self):
    """
    test compare_version returns True for identical versions
    """
    timer_start = perf_counter()
    self.assertTrue(compare_version('5.3.1', required_version='5.3.1'))
    print_timer(timer_start, caller='test_compare_version_equal')
  # ----------------------------------------------------------------------------
  def test_compare_version_newer(self):
    """
    test compare_version returns True when version is newer than required
    """
    timer_start = perf_counter()
    self.assertTrue(compare_version('5.4.0', required_version='5.3.1'))
    self.assertTrue(compare_version('6.0.0', required_version='5.3.1'))
    print_timer(timer_start, caller='test_compare_version_newer')
  # ----------------------------------------------------------------------------
  def test_compare_version_older(self):
    """
    test compare_version returns False when version is older than required
    """
    timer_start = perf_counter()
    self.assertFalse(compare_version('5.3.0', required_version='5.3.1'))
    self.assertFalse(compare_version('4.9.2', required_version='5.3.1'))
    print_timer(timer_start, caller='test_compare_version_older')
  # ----------------------------------------------------------------------------
  def test_compare_version_alpha_older_than_release(self):
    """
    test that an alpha version is older than the same release version
    """
    timer_start = perf_counter()
    self.assertFalse(compare_version('5.3.1-alpha01', required_version='5.3.1'))
    print_timer(timer_start, caller='test_compare_version_alpha_older_than_release')
  # ----------------------------------------------------------------------------

#===============================================================================
class print_stat_test_case(unittest.TestCase):
  """
  Tests for print_stat()
  """
  # ----------------------------------------------------------------------------
  def test_print_stat_returns_string(self):
    """
    test that print_stat returns a string
    """
    timer_start = perf_counter()
    x = np.array([1.0, 2.0, 3.0])
    result = print_stat(x, name='test_arr')
    self.assertIsInstance(result, str)
    print_timer(timer_start, caller='test_print_stat_returns_string')
  # ----------------------------------------------------------------------------
  def test_print_stat_min_in_output(self):
    """
    test that min value appears in output when stat includes 'n'
    """
    timer_start = perf_counter()
    x = np.array([1.0, 5.0, 3.0])
    msg = print_stat(x, name='x', stat='n')
    self.assertIn('min', msg)
    print_timer(timer_start, caller='test_print_stat_min_in_output')
  # ----------------------------------------------------------------------------
  def test_print_stat_max_in_output(self):
    """
    test that max value appears in output when stat includes 'x'
    """
    timer_start = perf_counter()
    x = np.array([1.0, 5.0, 3.0])
    msg = print_stat(x, name='x', stat='x')
    self.assertIn('max', msg)
    print_timer(timer_start, caller='test_print_stat_max_in_output')
  # ----------------------------------------------------------------------------
  def test_print_stat_avg_in_output(self):
    """
    test that avg value appears in output when stat includes 'a'
    """
    timer_start = perf_counter()
    x = np.array([2.0, 4.0])
    msg = print_stat(x, name='x', stat='a')
    self.assertIn('avg', msg)
    print_timer(timer_start, caller='test_print_stat_avg_in_output')
  # ----------------------------------------------------------------------------
  def test_print_stat_shape_in_output(self):
    """
    test that shape appears in output when stat includes 'h'
    """
    timer_start = perf_counter()
    x = np.ones((3, 4))
    msg = print_stat(x, name='x', stat='h')
    self.assertIn('shp', msg)
    print_timer(timer_start, caller='test_print_stat_shape_in_output')
  # ----------------------------------------------------------------------------
  def test_print_stat_scientific_format(self):
    """
    test that fmt='e' switches to scientific notation
    """
    timer_start = perf_counter()
    x = np.array([1.0e-10])
    msg = print_stat(x, name='x', stat='n', fmt='e')
    self.assertIn('e', msg.lower())
    print_timer(timer_start, caller='test_print_stat_scientific_format')
  # ----------------------------------------------------------------------------
  def test_print_stat_unit_in_output(self):
    """
    test that unit string appears in output when provided
    """
    timer_start = perf_counter()
    x = np.array([1.0])
    msg = print_stat(x, name='x', stat='n', unit='K')
    self.assertIn('K', msg)
    print_timer(timer_start, caller='test_print_stat_unit_in_output')
  # ----------------------------------------------------------------------------

#===============================================================================
class chk_finite_test_case(unittest.TestCase):
  """
  Tests for chk_finite()
  """
  # ----------------------------------------------------------------------------
  def test_chk_finite_valid_data_passes(self):
    """
    test that chk_finite does not raise for finite data
    """
    timer_start = perf_counter()
    x = xr.DataArray(np.array([1.0, 2.0, 3.0]))
    chk_finite(x)  # should not raise
    print_timer(timer_start, caller='test_chk_finite_valid_data_passes')
  # ----------------------------------------------------------------------------
  def test_chk_finite_raises_on_inf(self):
    """
    test that chk_finite raises ValueError when data contains inf
    """
    timer_start = perf_counter()
    x = xr.DataArray(np.array([1.0, np.inf, 3.0]))
    with self.assertRaises(ValueError):
      chk_finite(x)
    print_timer(timer_start, caller='test_chk_finite_raises_on_inf')
  # ----------------------------------------------------------------------------
  def test_chk_finite_raises_on_nan(self):
    """
    test that chk_finite raises ValueError when data contains nan
    """
    timer_start = perf_counter()
    x = xr.DataArray(np.array([1.0, np.nan, 3.0]))
    with self.assertRaises(ValueError):
      chk_finite(x)
    print_timer(timer_start, caller='test_chk_finite_raises_on_nan')
  # ----------------------------------------------------------------------------
  def test_chk_finite_error_message_contains_name(self):
    """
    test that ValueError message contains the provided name
    """
    timer_start = perf_counter()
    x = xr.DataArray(np.array([np.inf]))
    with self.assertRaises(ValueError) as ctx:
      chk_finite(x, name='my_var')
    self.assertIn('my_var', str(ctx.exception))
    print_timer(timer_start, caller='test_chk_finite_error_message_contains_name')
  # ----------------------------------------------------------------------------

#===============================================================================
if __name__ == '__main__':
  unittest.main()
