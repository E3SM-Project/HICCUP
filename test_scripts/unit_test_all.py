#!/usr/bin/env python
import unittest
from time import perf_counter
from hiccup.hiccup_data_class_timer_methods import print_timer
import unit_test_data_class
import unit_test_state_adjustment
import unit_test_memory_methods

timer_start = perf_counter()

loader = unittest.TestLoader()
runner = unittest.TextTestRunner(verbosity=1)

suite_list = []
suite_list.append( loader.loadTestsFromModule(unit_test_data_class) )
suite_list.append( loader.loadTestsFromModule(unit_test_state_adjustment) )
suite_list.append( loader.loadTestsFromModule(unit_test_memory_methods) )

for suite in suite_list: runner.run(suite)

print_timer(timer_start,caller='total time for all tests')
