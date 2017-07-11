import unittest

import unittest_Structure
import unittest_DXFile
import unittest_EnergyGrid
import unittest_MathTools


def set_up_test_suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.TestLoader().loadTestsFromTestCase(unittest_Structure.TestCases))
    suite.addTest(unittest.TestLoader().loadTestsFromTestCase(unittest_DXFile.TestCases))
    suite.addTest(unittest.TestLoader().loadTestsFromTestCase(unittest_EnergyGrid.TestCases))
    suite.addTest(unittest.TestLoader().loadTestsFromTestCase(unittest_MathTools.TestCases))
    return suite


def main():
    runner = unittest.TextTestRunner()
    test_suite = set_up_test_suite()
    runner.run(test_suite)


if __name__ == '__main__':
    main()

