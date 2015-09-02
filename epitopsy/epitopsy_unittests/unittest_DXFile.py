import unittest
import numpy
from epitopsy.DXFile import DXBox, VDWBox, read_dxfile, random_sphere
import os


class TestCases(unittest.TestCase):
    dx = read_dxfile(os.path.dirname(os.path.abspath(__file__))+'sample_data/gibbs_free_energy.dx',
                     'esp')

    def runTest(self):
        pass
        
    def test_DXBox_void_instation(self):
        # just to have something to start from
        dx = DXBox()

    def test_recalculate_radius():
        coords = random_sphere(dx.box.shape)
        offset = dx.box.shape[0]/2
        r = math.sqrt(sum([(coord-offset)**2 for coord in coords]))
        self.assertLess(abs(r-offset),3)

if __name__ == '__main__':
    unittest.main()
