# Do 18. Mai 17:09:00 CEST 2017

import os
import unittest
import numpy as np
from epitopsy.tools import MathTools

class TestCases(unittest.TestCase):
    def test_MathTools_get_euler_angles_for_equal_distributed_rotation(self):
        ref_N2 = np.array([[ 105.00000000,5.93173027e-05,  -75.00000000],
                           [ 106.96916088,   37.70676793,  106.96916088]])
        ref_N6 = np.array([[ 118.22134512,6.50270688e-05,  -61.77865488],
                           [  73.03083912,   37.70676793,   73.03083912],
                           [ -81.27658797,   83.5822834 ,   98.72341203],
                           [  84.05287422,   51.88144481,  -95.94712578],
                           [ -74.8977018 ,    8.98885041,  -74.8977018 ],
                           [-119.60887066,   19.95294743,   60.39112934]])
        ref_N128_sum = 27430.78353048778
        ref_N256_sum = 55213.78280619308
        
        result = MathTools.get_euler_angles_for_equal_distributed_rotation(2)
        self.assertTrue(np.allclose(result, ref_N2, rtol=1e-6, atol=0),
                        'Cannot reproduce a spherical lattice with 2 dots')
        result = MathTools.get_euler_angles_for_equal_distributed_rotation(6)
        self.assertTrue(np.allclose(result, ref_N6, rtol=1e-6, atol=0),
                        'Cannot reproduce a spherical lattice with 6 dots')
        result = np.sum(np.abs(
               MathTools.get_euler_angles_for_equal_distributed_rotation(128)))
        self.assertTrue(np.isclose(result, ref_N128_sum, rtol=1e-6, atol=0),
                        'Cannot reproduce a spherical lattice with 128 dots')
        result = np.sum(np.abs(
               MathTools.get_euler_angles_for_equal_distributed_rotation(256)))
        self.assertTrue(np.isclose(result, ref_N256_sum, rtol=1e-6, atol=0),
                        'Cannot reproduce a spherical lattice with 1024 dots')
    
    def test_MathTools_spherical_lattice(self):
        ref_N2 = np.array([[ 0.86602540, -0.5       ,  0.        ],
                           [-0.63858018,  0.5       ,  0.58499175]])
        ref_N6 = np.array([[ 0.5527708 , -0.83333333,  0.        ],
                           [-0.63858018, -0.5       ,  0.58499175],
                           [ 0.08620293, -0.16666667, -0.98223789],
                           [ 0.59992881,  0.16666667,  0.78250089],
                           [-0.85278689,  0.5       , -0.15084599],
                           [ 0.46640329,  0.83333333, -0.29668759]])
        ref_N128_sum = 192.03294495986
        ref_N1024_sum = 1536.01908381982
        
        result = MathTools.spherical_lattice_default(2)
        self.assertTrue(np.allclose(result, ref_N2, rtol=1e-6, atol=0),
                        'Cannot reproduce a spherical lattice with 2 dots')
        result = MathTools.spherical_lattice_default(6)
        self.assertTrue(np.allclose(result, ref_N6, rtol=1e-6, atol=0),
                        'Cannot reproduce a spherical lattice with 6 dots')
        result = np.sum(np.abs(MathTools.spherical_lattice_default(128)))
        self.assertTrue(np.isclose(result, ref_N128_sum, rtol=1e-6, atol=0),
                        'Cannot reproduce a spherical lattice with 128 dots')
        result = np.sum(np.abs(MathTools.spherical_lattice_default(1024)))
        self.assertTrue(np.isclose(result, ref_N1024_sum, rtol=1e-6, atol=0),
                        'Cannot reproduce a spherical lattice with 1024 dots')


if __name__ == '__main__':
    unittest.main()

