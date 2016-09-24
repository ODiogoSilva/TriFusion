import unittest

import trifusion.tests.test_load_process_data as load_process_data
import trifusion.tests.test_partitions as partitions
import trifusion.tests.test_process_filters as process_filters
import trifusion.tests.test_process_write as process_write
import trifusion.tests.test_secondary_ops as secondary_ops

loader = unittest.TestLoader()
suite = unittest.TestSuite()

suite.addTests(loader.loadTestsFromModule(load_process_data))
suite.addTests(loader.loadTestsFromModule(partitions))
suite.addTests(loader.loadTestsFromModule(process_filters))
suite.addTests(loader.loadTestsFromModule(process_write))
suite.addTests(loader.loadTestsFromModule(secondary_ops))

runner = unittest.TextTestRunner(verbosity=3)
result = runner.run(suite)
