"""
@author: Christoph Wilms
"""
import numpy as np

# "cimport" is used to import special compile-time information
# about the numpy module (this is stored in a file numpy.pxd which is
# currently part of the Cython distribution).
cimport numpy as np
cimport cython


# We now need to fix a datatype for our arrays. I've used the variable
# DTYPE for this, which is assigned to the usual NumPy runtime
# type info object.
DTYPE_int = np.int
DTYPE_float = np.float
DTYPE_cmplx = np.complex

# "ctypedef" assigns a corresponding compile-time type to DTYPE_t. For
# every type in the numpy module there's a corresponding compile-time
# type with a _t-suffix.
ctypedef np.int_t DTYPE_int_t
ctypedef np.float_t DTYPE_float_t
ctypedef np.complex_t DTYPE_cmplx_t

def fft_shift(np.ndarray[DTYPE_float_t, ndim=3] box_np, int x_center, int y_center,
              int z_center):
    cdef int xdim = box_np.shape[0]
    cdef int ydim = box_np.shape[1]
    cdef int zdim = box_np.shape[2]
    cdef int i,j,k
    cdef int x,y,z
    cdef np.ndarray[DTYPE_float_t, ndim=3] new_box_np = np.zeros([xdim, ydim, zdim], dtype=DTYPE_float)
    cdef double[:,:,:] box = box_np
    cdef double[:,:,:] new_box = new_box_np
    for i in range(xdim):
        for j in range(ydim):
            for k in range(zdim):
                if i < x_center:
                    x = i
                else:
                    x = i - xdim
                if j < y_center:
                    y = j
                else:
                    y = j - ydim
                if k < z_center:
                    z = k
                else:
                    z = k - zdim
                x = x_center - x - 1
                y = y_center - y - 1
                z = z_center - z - 1
                
                new_box[x, y, z] = box[i, j, k]
    
    return new_box_np
