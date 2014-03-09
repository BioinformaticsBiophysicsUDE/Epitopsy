import unittest

from epitopsy.Structure import *

class TestCases(unittest.TestCase):

    def test_PDBFile_instantiation(self):
        # just to have something to start from
        pdb = PDBFile()