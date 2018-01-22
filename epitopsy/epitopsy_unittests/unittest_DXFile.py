import os
import math
import shutil
import tarfile
import tempfile
import unittest
import numpy as np
from epitopsy.DXFile import read_dxfile, random_sphere, DXBox


class TestCases(unittest.TestCase):
    def setUp(self):
        # create an empty directory
        self.dirname = tempfile.mkdtemp(prefix='epitopsy-unittest-DXFile-',
                                        suffix='')
        # untar reference files
        reference_grids = 'sample_data/EnergyGrid/reference_grids.tar.bz2'
        with tarfile.open(reference_grids, 'r:bz2') as archive:
            archive.extract('protein1_esp.dx', self.dirname)
            archive.extract('protein1_vdw.dx', self.dirname)
    
    def tearDown(self):
        # delete empty directory
        shutil.rmtree(self.dirname)
    
    def test_read_dxfile(self):
        # integer OpenDX
        vdw = read_dxfile(os.path.join(self.dirname, 'protein1_vdw.dx'), 'vdw')
        vdw_sum = int(np.sum(np.abs(vdw.box)))
        self.assertEqual(vdw_sum, 35548, 'Grid protein1_vdw is different')
        
        # float OpenDX
        esp = read_dxfile(os.path.join(self.dirname, 'protein1_esp.dx'), 'esp')
        esp_sum = np.sum(np.abs(esp.box))
        self.assertTrue(np.isclose(esp_sum, 18250.837397, rtol=1e-5, atol=0),
                        'Grid for protein1_esp is different')
        
        # OpenDX file header
        ref_header = '# Data from 1.4.1\n# \n# POTENTIAL (kT/e)\n#'
        esp_header = esp.comment_header.rstrip()
        self.assertEqual(esp_header, ref_header,
                        'Read incorrect OpenDX file header')
    
    def test_prepare_for_geometric_matching(self):
        vdw = read_dxfile(os.path.join(self.dirname, 'protein1_vdw.dx'), 'vdw')
        vdw.prepare_for_geometric_matching(-15)
        vdw_sum = int(np.sum(np.abs(vdw.box)))
        self.assertEqual(vdw_sum, 6575,
                        'Incorrect prepare_for_geometric_matching')
    
    def test_recalculate_radius(self):
        dx = read_dxfile(os.path.join(self.dirname, 'protein1_vdw.dx'), 'vdw')
        coords = random_sphere(dx.box.shape)
        offset = dx.box.shape[0] / 2
        r = math.sqrt(sum([(coord - offset)**2 for coord in coords]))
        self.assertLess(abs(r - offset), 3)
    
    def test_DXBox_transform_real_to_box_space(self):
        err_fmt = 'Coordinates {} should be at {} in box space, got {}'
        dx = DXBox(box=np.zeros(3*[11]), meshsize=3*[.8], offset=3*[0])
        
        coord = [0.0, 0.0, 0.0]
        ref = (0, 0, 0)
        got = tuple(dx.transform_real_to_box_space(coord))
        self.assertEqual(ref, got, err_fmt.format('center', ref, got))
        
        coord = [+0.2, +0.39999, +0.40001]
        ref = (0, 0, 1)
        got = tuple(dx.transform_real_to_box_space(coord))
        self.assertEqual(ref, got, err_fmt.format(coord, ref, got))
        
        coord = [-0.2, -0.39999, -0.40001]
        ref = (0, 0, -1)
        got = tuple(dx.transform_real_to_box_space(coord))
        self.assertEqual(ref, got, err_fmt.format(coord, ref, got))
    
    def test_DXBox_transform_box_to_real_space(self):
        err_fmt = 'Grid point at {} should be at {} in Cartesian space, got {}'
        dx = DXBox(box=np.zeros(3*[11]), meshsize=3*[.8], offset=3*[0])
        
        coord = (0, 0, 0)
        ref = np.array([0., 0., 0.])
        got = dx.transform_box_to_real_space(coord)
        self.assertTrue(np.allclose(ref, got, rtol=0, atol=1e-5),
                        err_fmt.format(coord, list(ref), list(got)))
        
        coord = (0, -2, 1)
        ref = np.array([0., -1.6, 0.8])
        got = dx.transform_box_to_real_space(coord)
        self.assertTrue(np.allclose(ref, got, rtol=0, atol=1e-5),
                        err_fmt.format(coord, list(ref), list(got)))


if __name__ == '__main__':
    unittest.main()

