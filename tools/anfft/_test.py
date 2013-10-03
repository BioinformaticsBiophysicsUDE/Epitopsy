#   ANFFT is an FFT package for Python, using FFTW.
#   Copyright 2010 Andrew Collette.
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.


"""
    Testsuite for ANFFT
"""
import unittest
import numpy as np

from anfft import *

double_data = None
single_data = None 

class TestComplex(unittest.TestCase):

    """
        Complex dft
    """

    def test_roundtrip(self):
        """ Round-trip value and type checking """

        def cc(data, func, ifunc, *args, **kwds):
            dtype = {4: np.complex64, 8: np.complex128}[data.dtype.itemsize]

            out = func(data, *args, **kwds)
            self.assertEqual(out.dtype, dtype)

            out2 = ifunc(out, *args, **kwds)
            self.assert_(np.max(np.abs((out2-data)/data)) < 1e-4)
            self.assertEqual(out2.dtype, dtype)

        cc(single_data, fft, ifft)
        cc(single_data, fftn, ifftn)
        cc(single_data, fftn, ifftn, 1)
        cc(single_data, fftn, ifftn, 2)

        cc(double_data, fft, ifft)
        cc(double_data, fftn, ifftn)
        cc(double_data, fftn, ifftn, 1)
        cc(double_data, fftn, ifftn, 2)

class TestReal(unittest.TestCase):

    """
        Real dft
    """

    def test_roundtrip(self):
        """ Real-valued roundtrip and types """

        def cc(data, func, ifunc, *args, **kwds):
            cdtype = np.complex128
            rdtype = np.float64

            out = func(data, *args, **kwds)
            self.assertEqual(out.dtype, cdtype)
            self.assertEqual(out.shape[0:-1], data.shape[0:-1])
            self.assertEqual(out.shape[-1], data.shape[-1]//2 + 1)

            out2 = ifunc(out, *args, **kwds)
            self.assertEqual(out2.dtype, rdtype)
            self.assertEqual(out2.shape, data.shape)
            self.assert_(np.max(np.abs((out2-data)/data)) < 1e-4)

        cc(single_data, rfft, irfft)
        cc(single_data, rfftn, irfftn)
        cc(single_data, rfftn, irfftn, 1)
        cc(single_data, rfftn, irfftn, 2)

        cc(double_data, rfft, irfft)
        cc(double_data, rfftn, irfftn)
        cc(double_data, rfftn, irfftn, 1)
        cc(double_data, rfftn, irfftn, 2)


def test(**kwds):
    """ Run test suite.  Keywords are passed to TextTestRunner constructor. """

    global double_data, single_data
    double_data = np.random.random((32,512,512)) + 2.0
    single_data = double_data.astype('f4')

    suite = unittest.TestSuite()
    suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(TestComplex))
    suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(TestReal))
    runner = unittest.TextTestRunner(**kwds)
    runner.run(suite)

    double_data = None
    single_data = None



