import unittest

import unittest_Structure
import unittest_DXFile


def set_up_test_suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.TestLoader().loadTestsFromTestCase(unittest_Structure.TestCases))
    suite.addTest(unittest.TestLoader().loadTestsFromTestCase(unittest_DXFile.TestCases))
    return suite


def main():
    runner = unittest.TextTestRunner()
    test_suite = set_up_test_suite()
    runner.run(test_suite)


if __name__ == '__main__':
    main()
