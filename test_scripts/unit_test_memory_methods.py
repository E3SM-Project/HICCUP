#!/usr/bin/env python
#===================================================================================================
# Unit testing for hiccup_data_class_memory_methods module
#===================================================================================================
import unittest
from unittest.mock import patch, MagicMock
from time import perf_counter
from hiccup.hiccup_data_class_timer_methods import print_timer
from hiccup.hiccup_data_class_memory_methods import get_mem_usage, print_mem_usage

verbose_default = False # local verbosity default

#===============================================================================
class memory_methods_test_case(unittest.TestCase):
  """
  Tests for hiccup_data_class_memory_methods.py
  """
  # ----------------------------------------------------------------------------
  def _make_mem_info(self, rss_bytes, vms_bytes):
    mem_info = MagicMock()
    mem_info.rss = rss_bytes
    mem_info.vms = vms_bytes
    return mem_info
  # ----------------------------------------------------------------------------
  @patch('hiccup.hiccup_data_class_memory_methods.psutil.Process')
  def test_get_mem_usage_returns_tuple(self, mock_process):
    """
    test that get_mem_usage returns a 2-tuple
    """
    timer_start = perf_counter()
    mock_process.return_value.memory_info.return_value = self._make_mem_info(2e9, 4e9)
    result = get_mem_usage()
    self.assertIsInstance(result, tuple)
    self.assertEqual(len(result), 2)
    print_timer(timer_start, caller='test_get_mem_usage_returns_tuple')
  # ----------------------------------------------------------------------------
  @patch('hiccup.hiccup_data_class_memory_methods.psutil.Process')
  def test_get_mem_usage_converts_bytes_to_gb(self, mock_process):
    """
    test that get_mem_usage correctly converts bytes to GB
    """
    timer_start = perf_counter()
    mock_process.return_value.memory_info.return_value = self._make_mem_info(2_000_000_000, 8_000_000_000)
    rss_gb, vms_gb = get_mem_usage()
    self.assertAlmostEqual(rss_gb, 2.0, places=3)
    self.assertAlmostEqual(vms_gb, 8.0, places=3)
    print_timer(timer_start, caller='test_get_mem_usage_converts_bytes_to_gb')
  # ----------------------------------------------------------------------------
  @patch('hiccup.hiccup_data_class_memory_methods.psutil.Process')
  def test_get_mem_usage_zero_memory(self, mock_process):
    """
    test get_mem_usage with zero memory values
    """
    timer_start = perf_counter()
    mock_process.return_value.memory_info.return_value = self._make_mem_info(0, 0)
    rss_gb, vms_gb = get_mem_usage()
    self.assertEqual(rss_gb, 0.0)
    self.assertEqual(vms_gb, 0.0)
    print_timer(timer_start, caller='test_get_mem_usage_zero_memory')
  # ----------------------------------------------------------------------------
  @patch('hiccup.hiccup_data_class_memory_methods.get_mem_usage', return_value=(1.5, 3.0))
  def test_print_mem_usage_returns_three_tuple(self, _):
    """
    test that print_mem_usage returns a 3-tuple of (msg, rss_gb, vms_gb)
    """
    timer_start = perf_counter()
    result = print_mem_usage(msg='test', use_color=False)
    self.assertIsInstance(result, tuple)
    self.assertEqual(len(result), 3)
    print_timer(timer_start, caller='test_print_mem_usage_returns_three_tuple')
  # ----------------------------------------------------------------------------
  @patch('hiccup.hiccup_data_class_memory_methods.get_mem_usage', return_value=(1.5, 3.0))
  def test_print_mem_usage_rss_vms_values(self, _):
    """
    test that returned rss and vms values match get_mem_usage output
    """
    timer_start = perf_counter()
    msg, rss, vms = print_mem_usage(msg='test', use_color=False)
    self.assertAlmostEqual(rss, 1.5)
    self.assertAlmostEqual(vms, 3.0)
    print_timer(timer_start, caller='test_print_mem_usage_rss_vms_values')
  # ----------------------------------------------------------------------------
  @patch('hiccup.hiccup_data_class_memory_methods.get_mem_usage', return_value=(1.5, 3.0))
  def test_print_mem_usage_custom_msg(self, _):
    """
    test that custom msg string appears in the returned message
    """
    timer_start = perf_counter()
    msg, _, _ = print_mem_usage(msg='my_custom_label', use_color=False)
    self.assertIn('my_custom_label', msg)
    print_timer(timer_start, caller='test_print_mem_usage_custom_msg')
  # ----------------------------------------------------------------------------
  @patch('hiccup.hiccup_data_class_memory_methods.get_mem_usage', return_value=(1.5, 3.0))
  def test_print_mem_usage_custom_indent(self, _):
    """
    test that custom indent is applied to the returned message
    """
    timer_start = perf_counter()
    msg, _, _ = print_mem_usage(indent='>>>', msg='x', use_color=False)
    self.assertTrue(msg.startswith('>>>'))
    print_timer(timer_start, caller='test_print_mem_usage_custom_indent')
  # ----------------------------------------------------------------------------
  @patch('hiccup.hiccup_data_class_memory_methods.get_mem_usage', return_value=(1.5, 3.0))
  def test_print_mem_usage_default_indent(self, _):
    """
    test that default indent is two spaces
    """
    timer_start = perf_counter()
    msg, _, _ = print_mem_usage(msg='x', use_color=False)
    self.assertTrue(msg.startswith('  '))
    print_timer(timer_start, caller='test_print_mem_usage_default_indent')
  # ----------------------------------------------------------------------------
  @patch('hiccup.hiccup_data_class_memory_methods.get_mem_usage', return_value=(2.0, 6.0))
  def test_print_mem_usage_values_in_string(self, _):
    """
    test that memory values are formatted into the returned string
    """
    timer_start = perf_counter()
    msg, _, _ = print_mem_usage(msg='x', use_color=False)
    self.assertIn('2.00', msg)
    self.assertIn('6.00', msg)
    print_timer(timer_start, caller='test_print_mem_usage_values_in_string')
  # ----------------------------------------------------------------------------
  @patch('hiccup.hiccup_data_class_memory_methods.get_mem_usage', return_value=(1.0, 2.0))
  def test_print_mem_usage_color_escape_codes(self, _):
    """
    test that use_color=True adds ANSI escape codes and use_color=False does not
    """
    timer_start = perf_counter()
    msg_color, _, _ = print_mem_usage(msg='x', use_color=True)
    msg_plain, _, _ = print_mem_usage(msg='x', use_color=False)
    self.assertIn('\033[', msg_color)
    self.assertNotIn('\033[', msg_plain)
    print_timer(timer_start, caller='test_print_mem_usage_color_escape_codes')
  # ----------------------------------------------------------------------------
  @patch('hiccup.hiccup_data_class_memory_methods.get_mem_usage', return_value=(1.0, 2.0))
  def test_print_mem_usage_default_msg_uses_caller_name(self, _):
    """
    test that when msg=None the caller function name is used in the output
    """
    timer_start = perf_counter()
    msg, _, _ = print_mem_usage(use_color=False)
    self.assertIn('test_print_mem_usage_default_msg_uses_caller_name', msg)
    print_timer(timer_start, caller='test_print_mem_usage_default_msg_uses_caller_name')
  # ----------------------------------------------------------------------------

#===============================================================================
if __name__ == '__main__':
  unittest.main()
