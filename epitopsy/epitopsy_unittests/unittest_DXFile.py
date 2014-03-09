import unittest

from epitopsy.DXFile import *

class TestCases(unittest.TestCase):

    def test_DXBox_instantiation(self):
        # just to have something to start from
        dx = DXBox()