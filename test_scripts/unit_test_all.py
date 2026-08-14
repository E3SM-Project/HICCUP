#!/usr/bin/env python
import sys
import unittest
from time import perf_counter
from hiccup.hiccup_data_class_timer_methods import print_timer
import unit_test_data_class
import unit_test_state_adjustment
import unit_test_memory_methods
import unit_test_timer_methods
import unit_test_utilities
import unit_test_vertical_remap

timer_start = perf_counter()

loader = unittest.TestLoader()
runner = unittest.TextTestRunner(verbosity=1)

suite_list = []
suite_list.append( loader.loadTestsFromModule(unit_test_data_class) )
suite_list.append( loader.loadTestsFromModule(unit_test_state_adjustment) )
suite_list.append( loader.loadTestsFromModule(unit_test_memory_methods) )
suite_list.append( loader.loadTestsFromModule(unit_test_timer_methods) )
suite_list.append( loader.loadTestsFromModule(unit_test_utilities) )
suite_list.append( loader.loadTestsFromModule(unit_test_vertical_remap) )

all_successful = True
for suite in suite_list:
  result = runner.run(suite)
  if not result.wasSuccessful(): all_successful = False

print_timer(timer_start,caller='total time for all tests')

# Return a non-zero exit code on failure so CI can detect it
if not all_successful: sys.exit(1)
