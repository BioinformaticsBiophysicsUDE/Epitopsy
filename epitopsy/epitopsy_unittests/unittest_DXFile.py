import os
import math
import shutil
import tarfile
import tempfile
import unittest
import numpy as np
from epitopsy.DXFile import read_dxfile, random_sphere


class TestCases(unittest.TestCase):
    def setUp(self):
        # create an empty directory
        self.dirname = tempfile.mkdtemp(prefix='epitopsy-unittest-DXFile-',
                                        suffix='')
        # untar reference files
        reference_grids = 'sample_data/EnergyGrid/reference_grids.tar.bz2'
        with tarfile.open(reference_grids, 'r:bz2') as archive:
            archive.extract('protein_esp.dx', self.dirname)
            archive.extract('protein_vdw.dx', self.dirname)
    
    def tearDown(self):
        # delete empty directory
        shutil.rmtree(self.dirname)
    
    def test_read_dxfile(self):
        # integer OpenDX
        vdw = read_dxfile(os.path.join(self.dirname, 'protein_vdw.dx'), 'vdw')
        vdw_sum = int(np.sum(np.abs(vdw.box)))
        self.assertTrue(vdw_sum == 35574,
                        'Read incorrect values from an integer OpenDX grid')
        
        # float OpenDX
        esp = read_dxfile(os.path.join(self.dirname, 'protein_esp.dx'), 'esp')
        esp_sum = np.sum(np.abs(esp.box))
        self.assertTrue(np.isclose(esp_sum, 16581.748724, rtol=1e-5, atol=0),
                        'Read incorrect values from a float OpenDX grid')
        
        # OpenDX file header
        ref = '# Data from 1.4.1\n# \n# POTENTIAL (kT/e)\n# \n'
        self.assertTrue(esp.comment_header.rstrip() == ref.rstrip(),
                        'Read incorrect OpenDX file header')
    
    def test_prepare_for_geometric_matching(self):
        vdw = read_dxfile(os.path.join(self.dirname, 'protein_vdw.dx'), 'vdw')
        vdw.prepare_for_geometric_matching(-15)
        vdw_sum = int(np.sum(np.abs(vdw.box)))
        self.assertTrue(vdw_sum == 6109,
                        'Incorrect prepare_for_geometric_matching')
    
    def test_recalculate_radius(self):
        dx = read_dxfile(os.path.join(self.dirname, 'protein_vdw.dx'), 'vdw')
        coords = random_sphere(dx.box.shape)
        offset = dx.box.shape[0]/2
        r = math.sqrt(sum([(coord-offset)**2 for coord in coords]))
        self.assertLess(abs(r-offset),3)


if __name__ == '__main__':
    unittest.main()

