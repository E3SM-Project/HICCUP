#!/usr/bin/env python
#===================================================================================================
# Unit testing for hiccup_data_class_timer_methods module
#===================================================================================================
import unittest
from unittest.mock import patch
from time import perf_counter
from hiccup.hiccup_data_class_timer_methods import print_timer, print_timer_summary

verbose_default = False # local verbosity default

#===============================================================================
class timer_methods_test_case(unittest.TestCase):
  """
  Tests for hiccup_data_class_timer_methods.py
  """
  # ----------------------------------------------------------------------------
  def test_print_timer_returns_string(self):
    """
    test that print_timer returns a string
    """
    timer_start = perf_counter()
    result = print_timer(perf_counter()-1, caller='test', use_color=False)
    self.assertIsInstance(result, str)
    print_timer(timer_start, caller='test_print_timer_returns_string')
  # ----------------------------------------------------------------------------
  def test_print_timer_caller_in_output(self):
    """
    test that the caller name appears in the returned message
    """
    timer_start = perf_counter()
    msg = print_timer(perf_counter()-1, caller='my_caller', use_color=False)
    self.assertIn('my_caller', msg)
    print_timer(timer_start, caller='test_print_timer_caller_in_output')
  # ----------------------------------------------------------------------------
  def test_print_timer_elapsed_time_in_output(self):
    """
    test that elapsed time in seconds appears in the returned message
    """
    timer_start = perf_counter()
    msg = print_timer(perf_counter()-5, caller='test', use_color=False)
    self.assertIn('sec', msg)
    print_timer(timer_start, caller='test_print_timer_elapsed_time_in_output')
  # ----------------------------------------------------------------------------
  def test_print_timer_elapsed_time_positive(self):
    """
    test that elapsed time is positive
    """
    timer_start = perf_counter()
    start = perf_counter() - 2.0
    msg = print_timer(start, caller='test', use_color=False)
    # extract numeric value from message - should contain a positive number
    import re
    match = re.search(r'(\d+\.\d+)\s+sec', msg)
    self.assertIsNotNone(match)
    self.assertGreater(float(match.group(1)), 0.0)
    print_timer(timer_start, caller='test_print_timer_elapsed_time_positive')
  # ----------------------------------------------------------------------------
  def test_print_timer_shows_minutes_when_over_60s(self):
    """
    test that elapsed time over 60 seconds also shows minutes
    """
    timer_start = perf_counter()
    msg = print_timer(perf_counter()-90, caller='test', use_color=False)
    self.assertIn('min', msg)
    print_timer(timer_start, caller='test_print_timer_shows_minutes_when_over_60s')
  # ----------------------------------------------------------------------------
  def test_print_timer_no_minutes_when_under_60s(self):
    """
    test that elapsed time under 60 seconds does not show minutes
    """
    timer_start = perf_counter()
    msg = print_timer(perf_counter()-5, caller='test', use_color=False)
    self.assertNotIn('min', msg)
    print_timer(timer_start, caller='test_print_timer_no_minutes_when_under_60s')
  # ----------------------------------------------------------------------------
  def test_print_timer_default_indent(self):
    """
    test that default indent is two spaces (newline + two spaces in printed output)
    """
    timer_start = perf_counter()
    # The msg itself (before color/newline) should start with the caller padded string
    msg = print_timer(perf_counter()-1, caller='test', use_color=False)
    self.assertTrue(msg.startswith('test'))
    print_timer(timer_start, caller='test_print_timer_default_indent')
  # ----------------------------------------------------------------------------
  def test_print_timer_color_adds_escape_codes(self):
    """
    test that use_color=True adds ANSI escape codes and use_color=False does not
    """
    timer_start = perf_counter()
    start = perf_counter() - 1
    msg_color = print_timer(start, caller='test', use_color=True)
    msg_plain = print_timer(start, caller='test', use_color=False)
    self.assertIn('\033[', msg_color)
    self.assertNotIn('\033[', msg_plain)
    print_timer(timer_start, caller='test_print_timer_color_adds_escape_codes')
  # ----------------------------------------------------------------------------
  @patch('builtins.print')
  def test_print_timer_print_called(self, mock_print):
    """
    test that print is called when print_msg=True (default)
    """
    timer_start = perf_counter()
    print_timer(perf_counter()-1, caller='test', use_color=False, print_msg=True)
    mock_print.assert_called_once()
    print_timer(timer_start, caller='test_print_timer_print_called')
  # ----------------------------------------------------------------------------
  @patch('builtins.print')
  def test_print_timer_print_suppressed(self, mock_print):
    """
    test that print is not called when print_msg=False
    """
    timer_start = perf_counter()
    print_timer(perf_counter()-1, caller='test', use_color=False, print_msg=False)
    mock_print.assert_not_called()
    print_timer(timer_start, caller='test_print_timer_print_suppressed')
  # ----------------------------------------------------------------------------
  def test_print_timer_default_caller_uses_frame_name(self):
    """
    test that when caller=None the calling function name is used in the output
    """
    timer_start = perf_counter()
    msg = print_timer(perf_counter()-1, use_color=False)
    self.assertIn('test_print_timer_default_caller_uses_frame_name', msg)
    print_timer(timer_start, caller='test_print_timer_default_caller_uses_frame_name')
  # ----------------------------------------------------------------------------
  @patch('builtins.print')
  def test_print_timer_summary_prints_messages(self, mock_print):
    """
    test that print_timer_summary prints each message in timer_msg_all
    """
    timer_start = perf_counter()
    msgs = ['msg_one', 'msg_two', 'msg_three']
    print_timer_summary(timer_msg_all=msgs)
    printed = ' '.join(str(c) for c in mock_print.call_args_list)
    for m in msgs:
      self.assertIn(m, printed)
    print_timer(timer_start, caller='test_print_timer_summary_prints_messages')
  # ----------------------------------------------------------------------------
  @patch('builtins.print')
  def test_print_timer_summary_no_output_when_none(self, mock_print):
    """
    test that print_timer_summary does nothing when both args are None
    """
    timer_start = perf_counter()
    print_timer_summary(timer_start_total=None, timer_msg_all=None)
    mock_print.assert_not_called()
    print_timer(timer_start, caller='test_print_timer_summary_no_output_when_none')
  # ----------------------------------------------------------------------------

#===============================================================================
if __name__ == '__main__':
  unittest.main()
