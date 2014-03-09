import unittest

from epitopsy.Structure import PDBFile


class TestCases(unittest.TestCase):

    def runTest(self):
        pass
        
    def test_PDBFile_void_instation(self):
        # just to have something to start from
        pdb = PDBFile()
