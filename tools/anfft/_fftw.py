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
    Internal module for talking to FFTW (via ctypes).
"""

from __future__ import with_statement

import ctypes
from warnings import warn
import os.path as op
import numpy as np
import os, sys

# Define FFTW constants
FFTW_FORWARD = -1
FFTW_BACKWARD = 1
FFTW_MEASURE = 0
FFTW_ESTIMATE = 64
FFTW_WISDOM_ONLY = 1<<21    # Does not seem to work

# Types
from ctypes import c_int, c_uint
from numpy.ctypeslib import ndpointer
c_int_p = ctypes.POINTER(c_int)
c_fftw_plan = ctypes.c_void_p

# Globals.  These may depend on the platform.
PROTO = ctypes.CFUNCTYPE
DLL = ctypes.CDLL
MIN_THREAD = 1024*1024  # Don't use threads for arrays with fewer elements

# CPU count.  16 threads max in case we get something very large by mistake.
NCPUS = -1
try:
    import multiprocessing
    NCPUS = multiprocessing.cpu_count()
except ImportError:
    if sys.platform.startswith('win'):
        try:
            NCPUS = int(os.environ['NUMBER_OF_PROCESSORS'])
        except (KeyError, ValueError):
            pass
    else:
        try:
            NCPUS = int(os.sysconf('SC_NPROCESSORS_ONLN'))
        except (AttributeError, ValueError):
            pass
if NCPUS <= 0:
    NCPUS = 2
elif NCPUS > 16:
    NCPUS = 16

# Custom FFTW path
FFTW_LIBRARY_PATH = os.environ.get('FFTW_PATH')
if FFTW_LIBRARY_PATH is not None:
    FFTW_LIBRARY_PATH = op.abspath(op.expanduser(FFTW_LIBRARY_PATH))
    if not op.exists(FFTW_LIBRARY_PATH):
        raise ValueError("FFTW_LIBRARY_PATH \"%s\" does not exist" % FFTW_LIBRARY_PATH)
    if op.exists(FFTW_LIBRARY_PATH) and not op.isdir(FFTW_LIBRARY_PATH):
        raise ValueError("FFTW_LIBRARY_PATH \"%s\" exists but is not a directory" % FFTW_LIBRARY_PATH)

def find_fftw_library(shortname):
    """ Try to find a the full path to a library.

    Checks FFTW_LIBRARY_PATH exclusively if set, otherwise uses the ctypes
    default mechanism.  Returns None if the library can't be found.
    """
    import re
    from ctypes import util
    if FFTW_LIBRARY_PATH is None:
        return util.find_library(shortname)

    if sys.platform.startswith('win'):
        raise NotImplementedError
    elif sys.platform.startswith('darwin'):
        raise NotImplementedError
    else:
        regexp = re.compile(r'^lib%s.so' % shortname)

    candidates = [x for x in os.listdir(FFTW_LIBRARY_PATH) if regexp.match(x)]
    if len(candidates) == 0:
        return None
    candidates.sort(key=lambda x: len(x))   # Prefer libfoo.so to libfoo.so.X.Y.Z
    return op.abspath(op.join(FFTW_LIBRARY_PATH, candidates[0]))

class BaseAPI(object):

    """
        Base class, for all flavors of the FFTW library.

        Subclasses override the class-level attributes with their
        API prototypes and library names.
    """

    api = {}        # These map (api_name, c_name) -> ctypes_prototype
    thread_api = {}
    shortname = ''  # e.g. 'fftw3f'
    savename = ''   # e.g. 'single.wis'
    have_threads = False    # Autodetected in constructor
    dtype = ''      # Complex dtype for this precision

    def __init__(self):

        fullname = find_fftw_library(self.shortname)
        if fullname is None:
            raise ImportError("Can't find required library %s" % self.shortname)
        library = DLL(fullname, mode=ctypes.RTLD_GLOBAL)
        self._bind_api(self.api, library)

        # The threading interface is optional but recommended
        thread_shortname = self.shortname+'_threads'
        fullname = find_fftw_library(thread_shortname)
        if fullname is not None:
            library = DLL(fullname, mode=ctypes.RTLD_GLOBAL)
            self._bind_api(self.thread_api, library)
            self.api_init_threads()
            self.savename = 'threads-'+self.savename
            self.have_threads = True
        else:
            warn("Can't load threaded library %s... performance will suffer" % thread_shortname)

    def _bind_api(self, api, library):
        for (apiname, funcname), PROTO in api.iteritems():
            setattr(self, apiname, PROTO((funcname, library)))

    def plan_real_dft(self, src, dest, direction, k, planmode):
        """ Plan a real fft """
        forward = (direction == FFTW_FORWARD)
        assert len(src.shape) == len(dest.shape)
        assert src.shape[0:-1] == dest.shape[0:-1]

        rank = len(src.shape)

        if forward:
            real_shape = src.shape
            complex_shape = dest.shape
        else:
            real_shape = dest.shape
            complex_shape = src.shape

        assert complex_shape[-1] == real_shape[-1]//2 + 1

        trans_shape = real_shape[rank-k:]
        trans_num = np.product(real_shape[0:rank-k], dtype='i')

        ct_rank = c_int(len(trans_shape))
        ct_dims = (c_int*len(trans_shape))(*trans_shape)
        ct_howmany = c_int(trans_num)

        ct_real_dist = c_int(np.product(real_shape[rank-k:], dtype='i'))
        ct_complex_dist = c_int(np.product(complex_shape[rank-k:], dtype='i'))

        if self.have_threads:
            self.api_plan_with_nthreads(1 if np.product(real_shape) < MIN_THREAD else NCPUS)

        if forward:
            return self.api_plan_many_dft_r2c(
                    ct_rank, ct_dims, ct_howmany,
                    src, None, 1, ct_real_dist,
                    dest, None, 1, ct_complex_dist,
                    planmode)
        else:
            return self.api_plan_many_dft_c2r(
                    ct_rank, ct_dims, ct_howmany,
                    src, None, 1, ct_complex_dist,
                    dest, None, 1, ct_real_dist,
                    planmode)

    def plan_dft(self, src, dest, direction, k, planmode):
        """ Plan a multidimensional FFT across the last k axes of the array.
        """
        rank = len(src.shape)
        assert 0 < k <= rank, k
        assert src.shape == dest.shape
        shape = src.shape

        trans_shape = shape[rank-k:]
        trans_num = np.product(shape[0:rank-k], dtype='i')

        ct_rank = c_int(len(trans_shape))
        ct_dims = (c_int*len(trans_shape))(*trans_shape)
        ct_howmany = c_int(trans_num)
        ct_dist = c_int(np.product(trans_shape, dtype='i'))

        if self.have_threads:
            self.api_plan_with_nthreads(1 if np.product(src.shape) < MIN_THREAD else NCPUS)

        return self.api_plan_many_dft(
                ct_rank, ct_dims, ct_howmany,
                src, None, 1, ct_dist,
                dest, None, 1, ct_dist,
                direction, planmode)
        
    def get_buffer(self, shape, alignment=0, real=False):
        """ Obtain a buffer of complex type allocated from "fast" memory.

        Inputs: shape tuple, optionally "alignment" which is address mod 16.

        Returns an array of complex type, unless "real" is True.
        """
        # Mechanism inspired by pyfftw
        dtype = np.dtype('f%d' % (self.dtype.itemsize/2)) if real else self.dtype
        nbytes = np.product(shape, dtype='i')*dtype.itemsize
        base = np.empty(nbytes+16, dtype='u1')
        base_offset = base.ctypes.data%16
        if alignment > base_offset:
            offset = alignment-base_offset
        else:
            offset = (alignment+16)-base_offset
        return base[offset:offset+nbytes].view(dtype=dtype).reshape(shape)

    def execute(self, plan):
        """ Execute a plan """
        assert plan is not None
        self.api_execute(plan)

    def execute_dft(self, plan, src, dest):
        """ Apply a plan to a new set of complex arrays """
        assert plan is not None
        self.api_execute_dft(plan, src, dest)

    def execute_real_dft(self, plan, src, dest, direction):
        """ Apply a plan to a new set of real arrays """
        assert plan is not None
        assert src.shape[0:-1] == dest.shape[0:-1]
        if direction == FFTW_FORWARD:
            assert src.shape[-1]//2 + 1 == dest.shape[-1]
            self.api_execute_dft_r2c(plan, src, dest)
        else:
            assert dest.shape[-1]//2 + 1 == src.shape[-1]
            self.api_execute_dft_c2r(plan, src, dest)

    def save_wisdom(self, basedir):
        """ Save wisdom to a file in basedir """
        wisdom = self.api_export_wisdom_to_string()
        with open(op.join(basedir, self.savename), 'wb') as wfile:
            wfile.write(wisdom)

    def load_wisdom(self, basedir):
        """ Try to load wisdom from a file in basedir.
        """
        try:
            wfile = open(op.join(basedir, self.savename), 'rb')
        except IOError:
            warn("Can't load wisdom from \"%s\"" % basedir)
        else:
            try:
                wisdom = wfile.read()
                result = self.api_import_wisdom_from_string(wisdom)
            finally:
                wfile.close()                
            if not result:
                raise SystemError("Error importing wisdom")

wf = ('contiguous', 'writeable')    # Write flags
rf = ('contiguous',)                # Read flags

class SingleAPI(BaseAPI):

    """
        Represents the single-precision library
    """

    api = {
        ('api_plan_many_dft', 'fftwf_plan_many_dft'):
            PROTO(c_fftw_plan, c_int, c_int_p, c_int,
            ndpointer('c8', flags=rf), c_int_p, c_int, c_int,
            ndpointer('c8', flags=wf), c_int_p, c_int, c_int,
            c_int, c_uint),
        ('api_execute', 'fftwf_execute'):
            PROTO(None, c_fftw_plan),
        ('api_export_wisdom_to_string', 'fftwf_export_wisdom_to_string'):
            PROTO(ctypes.c_char_p),
        ('api_import_wisdom_from_string', 'fftwf_import_wisdom_from_string'):
            PROTO(c_int, ctypes.c_char_p),
        ('api_set_timelimit', 'fftwf_set_timelimit'):
            PROTO(None, ctypes.c_double),
        ('api_execute_dft', 'fftwf_execute_dft'):
            PROTO(None, c_fftw_plan,
            ndpointer('c8', flags=rf), ndpointer('c8', flags=wf)),}

    thread_api = {
        ('api_init_threads', 'fftwf_init_threads'):
            PROTO(c_int),
        ('api_plan_with_nthreads', 'fftwf_plan_with_nthreads'):
            PROTO(None, c_int),}

    shortname = 'fftw3f'
    savename = 'single.wis'
    dtype = np.dtype('c8')

class DoubleAPI(BaseAPI):

    """
        Represents the double-precision library
    """

    api = {
        ('api_plan_many_dft', 'fftw_plan_many_dft'):
            PROTO(c_fftw_plan, c_int, c_int_p, c_int,
            ndpointer('c16', flags=rf), c_int_p, c_int, c_int,
            ndpointer('c16', flags=wf), c_int_p, c_int, c_int,
            c_int, c_uint),
        ('api_plan_many_dft_r2c', 'fftw_plan_many_dft_r2c'):
            PROTO(c_fftw_plan, c_int, c_int_p, c_int,
            ndpointer('f8', flags=rf), c_int_p, c_int, c_int,
            ndpointer('c16', flags=wf), c_int_p, c_int, c_int,
            c_uint),
        ('api_plan_many_dft_c2r', 'fftw_plan_many_dft_c2r'):
            PROTO(c_fftw_plan, c_int, c_int_p, c_int,
            ndpointer('c16', flags=rf), c_int_p, c_int, c_int,
            ndpointer('f8', flags=wf), c_int_p, c_int, c_int,
            c_uint),
        ('api_execute', 'fftw_execute'):
            PROTO(None, c_fftw_plan),
        ('api_export_wisdom_to_string', 'fftw_export_wisdom_to_string'):
            PROTO(ctypes.c_char_p),
        ('api_import_wisdom_from_string', 'fftw_import_wisdom_from_string'):
            PROTO(c_int, ctypes.c_char_p),
        ('api_set_timelimit', 'fftw_set_timelimit'):
            PROTO(None, ctypes.c_double),
        ('api_execute_dft', 'fftw_execute_dft'):
            PROTO(None, c_fftw_plan,
            ndpointer('c16', flags=rf), ndpointer('c16', flags=wf)),
        ('api_execute_dft_r2c', 'fftw_execute_dft_r2c'):
            PROTO(None, c_fftw_plan,
            ndpointer('f8', flags=rf), ndpointer('c16', flags=wf)),
        ('api_execute_dft_c2r', 'fftw_execute_dft_c2r'):
            PROTO(None, c_fftw_plan,
            ndpointer('c16', flags=rf), ndpointer('f8', flags=wf)),}

    thread_api = {
        ('api_init_threads', 'fftw_init_threads'):
            PROTO(c_int),
        ('api_plan_with_nthreads', 'fftw_plan_with_nthreads'):
            PROTO(None, c_int),}

    shortname = 'fftw3'
    savename = 'double.wis'
    dtype = np.dtype('c16')


