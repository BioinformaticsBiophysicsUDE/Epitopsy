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

import os
import os.path as op
import atexit
import functools
import numpy as np
import _fftw

single_api = _fftw.SingleAPI()
double_api = _fftw.DoubleAPI()

WISDOM_PATH = op.expanduser('~/.anfft')
#if not op.exists(WISDOM_PATH):
#    os.makedirs(WISDOM_PATH)
#
#single_api.load_wisdom(WISDOM_PATH)
#double_api.load_wisdom(WISDOM_PATH)
#
#atexit.register(functools.partial(single_api.save_wisdom, WISDOM_PATH))
#atexit.register(functools.partial(double_api.save_wisdom, WISDOM_PATH))

plancache = {}
rplancache = {}

def _fftn(data, direction, k, **kwds):
    """ Peform an FFT """

    if data.dtype in (np.float32, np.complex64):
        api = single_api
    elif data.dtype in (np.float64, np.complex128):
        api = double_api
    else:
        raise TypeError("Only single and double precision are supported")

    measure = bool(kwds.get('measure'))

    # Out-of-place transforms (preferred) only work when the data is
    # contiguous and of complex type.  These requirements could also be
    # extended to demand a particular data alignment.
    inplace = (not data.flags.c_contiguous) or \
              (not (data.dtype in (np.complex64, np.complex128)))

    # This is 0,1,2,3 for unaligned and 4, 8 and 16 bytes respectively.
    relative_alignment = 3 if inplace else \
                         sum(int(data.ctypes.data%x==0) for x in (4,8,16))

    key = ( data.shape,
            data.dtype,
            direction,
            k,
            inplace,
            relative_alignment )

    plan, measured = plancache.get(key, (None,None))

    # The output array is aligned to 16 bytes by default
    dest = api.get_buffer(data.shape)

    # We can skip the planning steps entirely if:
    # (1) A plan with this key already exists
    # (2) Measurement not required, or the plan was created by measuring
    if plan is not None and ((not measure) or measured):
        if inplace:
            src = dest
            src[:] = data
        else:
            src = data
        api.execute_dft(plan, src, dest)
        return dest
    
    # No plan, but we can estimate
    if plan is None and (not measure):
        if inplace:
            src = dest
            src[:] = data
        else:
            src = data
        plan = api.plan_dft(src, dest, direction, k, _fftw.FFTW_ESTIMATE)
        if plan is None:
            raise SystemError("Planning failed")
        plancache[key] = (plan, False)
        api.execute(plan)
        return dest

    # No pre-existing plan, and we need to measure

    if inplace:
        src = dest
        src[:] = data
    else:
        # Make sure to measure against an array with the same alignment
        src = api.get_buffer(data.shape, data.ctypes.data%16)
        src[:] = data

    plan = api.plan_dft(src, dest, direction, k, _fftw.FFTW_MEASURE)
    if plan is None:
        raise SystemError("Planning failed")
    plancache[key] = (plan, True)

    if inplace:
        src = dest
        src[:] = data
    else:
        src = data
    api.execute_dft(plan, src, dest)
    return dest

def _rfftn(data, direction, k, **kwds):

    forward = (direction==_fftw.FFTW_FORWARD)

    api = double_api
    measure = bool(kwds.get('measure'))

    destshape = list(data.shape)
    if forward:
        destshape[-1] = (destshape[-1]//2) + 1
    else:
        destshape[-1] = (destshape[-1]-1)*2
    destshape = tuple(destshape)

    # The real-data FFTs destroy their input, so we always make a copy.
    src = api.get_buffer(data.shape, real=forward)
    src[:] = data
    dest = api.get_buffer(destshape, real=(not forward))

    key = (data.shape, direction, k)

    plan, measured = rplancache.get(key,(None,None))

    if plan is not None and ((not measure) or measured):
        api.execute_real_dft(plan, src, dest, direction)
        return dest

    if plan is None and (not measure):
        plan = api.plan_real_dft(src, dest, direction, k, _fftw.FFTW_ESTIMATE)
        rplancache[key] = (plan, False)
        api.execute(plan)
        return dest

    plan = api.plan_real_dft(src, dest, direction, k, _fftw.FFTW_MEASURE)
    rplancache[key] = (plan, True)
    src[:] = data
    api.execute(plan)
    return dest

def fftn(data, k=None, **kwds):
    """ Compute the forward FFT over the last *k* axes (all by default).

    No normalization is performed for the forward FFT.
    """
    if k is None:
        k = len(data.shape)
    return _fftn(data, _fftw.FFTW_FORWARD, k, **kwds)

def ifftn(data, k=None, **kwds):
    """ Compute the inverse FFT over the last *k* axes (all by default).

    Normalization of 1/n is applied where *n* is the product of the
    last *k* dimensions.
    """
    if k is None:
        k = len(data.shape)
    out = _fftn(data, _fftw.FFTW_BACKWARD, k, **kwds)
    out /= np.product(data.shape[len(data.shape)-k:])
    return out

def fft(data, **kwds):
    """ Compute the forward FFT over the last axis only.

    No normalization is performed for the forward FFT.
    """
    return fftn(data, k=1, **kwds)

def ifft(data, **kwds):
    """ Compute the inverse FFT over the last axis only.

    Normalization of 1/n is applied where *n* is the size of the
    final dimension.
    """
    return ifftn(data, k=1, **kwds)

def rfftn(data, k=None, **kwds):
    """ Compute the forward FFT of real data over the last *k* axes.

    The returned array is complex with shape (n1 x n2 ... nd//2 + 1).
    Inputs should be real-valued; if complex, the real part will be taken.
    No normalization is applied to the forward transform.

    Input data must be real.  The calculation is always performed in double
    precision, and a double-precision result is returned.
    """
    if k is None:
        k = len(data.shape)
    return _rfftn(data, _fftw.FFTW_FORWARD, k, **kwds)

def irfftn(data, k=None, **kwds):
    """ Compute the inverse FFT from one-half of a Hermite-symmetric array.
    The calculation is performed over the last *k* axes (all by default).

    The input is complex and corresponds to the portion of an FFT with
    non-negative frequencies (the first half plus DC bin).  The returned
    array is real with shape (n1 x n2 ... 2*(nd-1)).  Normalization of
    1/n is applied, where *n* is the product of the last *k* dimensions
    of the returned array.

    The calculation is always performed in double precision, and a double-
    precision result is returned.
    """
    if k is None:
        k = len(data.shape)
    out = _rfftn(data, _fftw.FFTW_BACKWARD, k, **kwds)
    out /= np.product(out.shape[len(out.shape)-k:])
    return out

def rfft(data, **kwds):
    """ Compute the forward FFT of real data over the last axis only.

    The returned array is complex with shape (n1 x n2 ... nd//2 + 1).
    Inputs should be real-valued; if complex, the real part will be taken.
    No normalization is applied to the forward transform.

    The calculation is always performed in double precision, and a double-
    precision result is returned.
    """
    return rfftn(data, k=1, **kwds)

def irfft(data, **kwds):
    """ Compute the inverse FFT from one-half of a Hermite-symmetric array.
    The calculation is performed over the last axis only.

    The input is complex and corresponds to the portion of an FFT with
    non-negative frequencies (the first half plus DC bin).  The returned
    array is real with shape (n1 x n2 ... 2*(nd-1)).  Normalization of
    1/n is applied, where n is size of the last dimension of the returned
    array.

    The calculation is always performed in double precision, and a double-
    precision result is returned.
    """
    return irfftn(data, k=1, **kwds)


