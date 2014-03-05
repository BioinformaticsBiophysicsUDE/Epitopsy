"""
@author: Christoph Wilms
"""
import sys, time
import numpy as np

# "cimport" is used to import special compile-time information
# about the numpy module (this is stored in a file numpy.pxd which is
# currently part of the Cython distribution).
cimport numpy as np
cimport cython

from cython.parallel import prange, parallel
from cpython cimport bool

from libc.math cimport abs as abs_c, exp, round, sqrt, floor

from numpy cimport NPY_INT8 as NPY_int8
from numpy cimport NPY_INT16 as NPY_int16
from numpy cimport NPY_INT32 as NPY_int32
from numpy cimport NPY_INT64 as NPY_int64
from numpy cimport NPY_FLOAT16 as NPY_float16
from numpy cimport NPY_FLOAT32 as NPY_float32
from numpy cimport NPY_FLOAT64 as NPY_float64

int8 = np.dtype(np.int8)
int16 = np.dtype(np.int16)
int32 = np.dtype(np.int32)
int64 = np.dtype(np.int64)
float16 = np.dtype(np.float16)
float32 = np.dtype(np.float32)
float64 = np.dtype(np.float64)

@cython.boundscheck(False)
cdef tuple _count_rotations(np.ndarray[np.float64_t, ndim=3] esp_correlation_np,
        np.ndarray[np.float64_t, ndim=3] shape_correlation_np,
        np.ndarray[np.float64_t, ndim=3] counter_matrix_np,
        np.ndarray[np.float64_t, ndim=3] interaction_energy_np,
        float shape_correlation_cutoff):
    cdef double[:,:,:] shape_correlation = shape_correlation_np
    cdef double[:,:,:] esp_correlation = esp_correlation_np
    cdef double[:,:,:] counter_matrix = counter_matrix_np
    cdef double[:,:,:] interaction_energy = interaction_energy_np,
    cdef int x,y,z
    cdef float counter_add = 1.
    cdef int xdim = shape_correlation.shape[0]
    cdef int ydim = shape_correlation.shape[1]
    cdef int zdim = shape_correlation.shape[2]

    for x in range(xdim):
        for y in range(ydim):
            for z in range(zdim):
                # use absolute value of shape_correlation, so that it is
                # in the
                # range(-shape_correlation_cutoff, shape_correlation_cutoff)
                if abs_c(shape_correlation[x,y,z]) < shape_correlation_cutoff:
                    # increase counter_matrix at the current position
                    counter_matrix[x, y, z] = counter_matrix[x, y, z] + counter_add
                    # get energy if shape correlation is okay
                    # calc exp(-energy), it should be faster to do it here
#                    interaction_energy[x, y, z] = exp(-esp_correlation[x, y, z])
                    interaction_energy[x, y, z] = esp_correlation[x, y, z]
    return interaction_energy_np, counter_matrix_np

def count_rotations(np.ndarray[np.float64_t, ndim=3] shape_correlation,
        np.ndarray[np.float64_t, ndim=3] esp_correlation,
        np.ndarray[np.float64_t, ndim=3] counter_matrix,
        float shape_correlation_cutoff):

    interaction_energy = np.zeros([esp_correlation.shape[0],
        esp_correlation.shape[1], esp_correlation.shape[2]])
    interaction_energy, counter_matrix = _count_rotations(esp_correlation,
        shape_correlation, counter_matrix, interaction_energy,
        shape_correlation_cutoff)
    return interaction_energy, counter_matrix

#@cython.boundscheck(False)
#cdef tuple calculateEnergy(np.ndarray[np.float64_t, ndim=3] esp_correlation_np,
#        np.ndarray[np.float64_t, ndim=3] shape_correlation_np,
#        np.ndarray[np.float64_t, ndim=3] counter_matrix_np,
#        np.ndarray[np.float64_t, ndim=3] interaction_energy_np,
#        int x_center, int y_center, int z_center, float shape_correlation_cutoff):
#    cdef double[:,:,:] shape_correlation = shape_correlation_np
#    cdef double[:,:,:] esp_correlation = esp_correlation_np
#    cdef double[:,:,:] counter_matrix = counter_matrix_np
#    cdef double[:,:,:] interaction_energy = interaction_energy_np
#    cdef int i,j,k
#    cdef int x,y,z
#    cdef float counter_add = 1.
#    cdef int xdim = shape_correlation.shape[0]
#    cdef int ydim = shape_correlation.shape[1]
#    cdef int zdim = shape_correlation.shape[2]
#
#    for i in range(xdim):
#        for j in range(ydim):
#            for k in range(zdim):
#                """
#                correct for fft freq shift
#                """
#                if i < x_center:
#                    x = i
#                else:
#                    x = i - xdim
#                if j < y_center:
#                    y = j
#                else:
#                    y = j - ydim
#                if k < z_center:
#                    z = k
#                else:
#                    z = k - zdim
#
#                ## all are ints, and xdim is usually uneven, so there is
#                # probably an error here, this would produce an offset, but
#                # just of one or two grid points, so not important for the
#                # moment
#                x = x_center - x - 1
#                y = y_center - y - 1
#                z = z_center - z - 1
#
#                if abs_c(shape_correlation[i,j,k]) < shape_correlation_cutoff:
#                    # increase counter_matrix at the current position
#                    counter_matrix[x, y, z] = counter_matrix[x, y, z] + counter_add
#                    # get energy if shape correlation is okay
#                    interaction_energy[x, y, z] = esp_correlation[i, j, k]
#    return interaction_energy_np, counter_matrix_np
#
#
#
#def speed_up_fft(np.ndarray[np.float64_t, ndim=3] shape_correlation,
#        np.ndarray[np.float64_t, ndim=3] esp_correlation,
#        np.ndarray[np.float64_t, ndim=3] counter_matrix,
#        float shape_correlation_cutoff,
#        center = None):
#    cdef int xdim = shape_correlation.shape[0]
#    cdef int ydim = shape_correlation.shape[1]
#    cdef int zdim = shape_correlation.shape[2]
#
#    if center is None:
#        x_center = <int> round(shape_correlation.shape[0] / 2.)
#        y_center = <int> round(shape_correlation.shape[1] / 2.)
#        z_center = <int> round(shape_correlation.shape[2] / 2.)
#    else:
#        x_center = <int> center[0]
#        y_center = <int> center[1]
#        z_center = <int> center[2]
#
#    interaction_energy = np.zeros([xdim, ydim, zdim])
#    interaction_energy, counter_matrix = calculateEnergy(esp_correlation,
#        shape_correlation, counter_matrix, interaction_energy, x_center,
#        y_center, z_center, shape_correlation_cutoff)
#    return interaction_energy, counter_matrix



# This version should work fine, if one changes something, please comment these
# options
@cython.cdivision(True)
@cython.boundscheck(False)
cdef double[:] transform_box_to_real_space(double[:] grid_coord,
        double[:] box_mesh_size,
        double[:] box_offset):
    cdef double[:] new_grid_coord = np.zeros(3)
    cdef int dim = 3
    cdef int i
    for i in range(dim):
        new_grid_coord[i] = grid_coord[i] * box_mesh_size[i] + box_offset[i]

    return new_grid_coord

# This version should work fine, if one changes something, please comment these
# options
@cython.cdivision(True)
@cython.boundscheck(False)
cdef long[:] transform_real_to_box_space(double[:] grid_coord,
        double[:] box_mesh_size,
        double[:] box_offset):
    cdef long[:] new_grid_coord = np.zeros(3,dtype=np.int)
    cdef int dim = 3
    cdef int i
    for i in range(dim):
        new_grid_coord[i] = <int> round((grid_coord[i]-box_offset[i]) / box_mesh_size[i])

    return new_grid_coord

# This version should work fine, if one changes something, please comment these
# options
@cython.cdivision(True)
@cython.boundscheck(False)
cdef bool is_overlapp(double[:,:,:] vdwbox,
        long[:,:] ligand_coords,
        long[:] vdw_radii,
        double solvent_score,
        double protein_score):
    cdef int xdim = vdwbox.shape[0]
    cdef int ydim = vdwbox.shape[1]
    cdef int zdim = vdwbox.shape[2]
    cdef int n_coords = ligand_coords.shape[0]
    cdef int radius = 0
    cdef int i,j,k,l,nx,ny,nz,x,y,z
    for i in range(n_coords):
        x = ligand_coords[i,0]
        y = ligand_coords[i,1]
        z = ligand_coords[i,2]
        radius = vdw_radii[i]
        for j in range(-radius,radius+1):
            for k in range(-radius,radius+1):
                for l in range(-radius,radius+1):
                    nx = x + j
                    ny = y + k
                    nz = z + l
                    if floor(sqrt((nx-x)**2+(ny-y)**2+(nz-z)**2)) <= radius:
                        if 0 <= nx < xdim and 0 <= ny < ydim and 0 <= nz < zdim:
                            if vdwbox[nx,ny,nz] == protein_score:
                                return True

    return False

# This version should work fine, if one changes something, please comment these
# options
@cython.cdivision(True)
@cython.boundscheck(False)
cdef double score_energy(double[:,:,:] potbox,
        long[:,:] ligand_coords,
        double[:] charges):
    cdef int n_coords = ligand_coords.shape[0]
    cdef double energy = 0.
    #cdef double pot_value = 0. # JN: not used in this function, may be removed
    cdef long x,y,z
    cdef int i
    for i in range(n_coords):
        x = ligand_coords[i,0]
        y = ligand_coords[i,1]
        z = ligand_coords[i,2]
        energy += potbox[x,y,z] * charges[i]

    return energy

# This version should work fine, if one changes something, please comment these
# options
@cython.cdivision(True)
@cython.boundscheck(False)
cdef double[:,:] move_ligand(double[:,:] ligand_real_coords,
        double[:] coords):
    cdef double[:] mean = np.zeros(3)  # geometrical center of the ligand
    cdef double[:] shift = np.zeros(3) # a vector from the center of the ligand to the current grid point 
    cdef int n_coords = ligand_real_coords.shape[0]
    cdef double[:,:] new_ligand_real_coords = np.zeros([n_coords,3])
    cdef int i,j
    for i in range(n_coords):
        for j in range(3):
            mean[j] += ligand_real_coords[i,j]

    cdef double n_atoms = <double> n_coords
    for i in range(3):
        mean[i] = mean[i] / n_atoms
        shift[i] = coords[i] - mean[i]

    for i in range(n_coords):
        for j in range(3):
            new_ligand_real_coords[i,j] = ligand_real_coords[i,j] + shift[j]

    return new_ligand_real_coords

# This version should work fine, if one changes something, please comment these
# options
@cython.cdivision(True)
@cython.boundscheck(False)
cdef long[:,:] transform_ligand_to_box(double[:,:] ligand_real_coords,
        double[:] box_mesh_size,
        double[:] box_offset):
    cdef int n_coords = ligand_real_coords.shape[0]
    cdef long[:,:] new_ligand_coords = np.zeros([n_coords,3],dtype=np.int)
    cdef long[:] new_coords = np.zeros(3,dtype=np.int)
    cdef double[:] old_coords = np.zeros(3)
    cdef int i,j

    for i in range(n_coords):
        old_coords = ligand_real_coords[i,:]
        new_coords = transform_real_to_box_space(old_coords,
                box_mesh_size,
                box_offset)
        for j in range(3):
            new_ligand_coords[i,j] = new_coords[j]

    return new_ligand_coords

# This version should work fine, if one changes something, please comment these
# options
@cython.cdivision(True)
@cython.boundscheck(False)
cdef bool out_of_box(long[:,:] ligand_coords, long[:] box_dim):
    cdef int i,j
    cdef int n_coords = ligand_coords.shape[0]
    cdef int x = 0

    for i in range(n_coords):
        for j in range(3):
            x = ligand_coords[i,j]
            if x < 0 or x >= box_dim[j]:
                return True

    return False

def speed_up_brute_force(ligand_pqr,
        np.ndarray[np.float64_t, ndim=3] vdwbox_np,
        espbox,
        np.ndarray[np.float64_t, ndim=3] counter_matrix_np):
    ## general stuff
    cdef int xdim = espbox.box.shape[0]
    cdef int ydim = espbox.box.shape[1]
    cdef int zdim = espbox.box.shape[2]
    cdef int x,y,z
    cdef double energy = 0.
    cdef double solvent_score = 1.
    cdef double protein_score = 0.
    cdef bool overlapp = False
    cdef bool left_the_box = False

    cdef np.ndarray[np.float64_t, ndim=3] interaction_energy_np = np.zeros([xdim,ydim,zdim])
    cdef double[:,:,:] interaction_energy = interaction_energy_np
    cdef double[:,:,:] counter_matrix = counter_matrix_np
    cdef double[:,:,:] vdwbox = vdwbox_np

    ## espbox data
    cdef double[:] box_mesh_size = np.array(espbox.box_mesh_size)
    cdef double[:] box_offset = np.array(espbox.box_offset)
    cdef long[:] box_dim = np.array(espbox.box_dim)
    cdef double[:,:,:] potbox = espbox.box

    ## coord data
    cdef double[:] coord = np.zeros(3)

    all_atoms = ligand_pqr.get_all_atoms()
    coord_list = []
    charge_list = []
    vdw_list = []
    cdef int n_atoms = len(all_atoms)
    for atom in all_atoms:
        coord_list.append(atom.get_coord())
        charge_list.append(atom.get_info_1())
        vdw_list.append(round(atom.get_info_2() / box_mesh_size[0]))

    cdef double[:,:] ligand_real_coords = np.array(coord_list)
    cdef long[:,:] ligand_coords = np.zeros([n_atoms, 3],dtype=np.int)
    cdef double[:] charges = np.array(charge_list)
    cdef long[:] vdw_radii = np.array(vdw_list,dtype=np.int)

    for x in range(xdim):
        for y in range(ydim):
            for z in range(zdim):
                # not inside the protein
                if vdwbox[x, y, z] == solvent_score:
                    # transform box to real space coordinates
                    coord[0] = <double> x
                    coord[1] = <double> y
                    coord[2] = <double> z

                    coord = transform_box_to_real_space(coord,
                            box_mesh_size,
                            box_offset)

                    # move the ligand
                    ligand_real_coords = move_ligand(ligand_real_coords,
                            coord)

                    # transform ligand coordinates
                    ligand_coords = transform_ligand_to_box(ligand_real_coords,
                            box_mesh_size,
                            box_offset)


                    # check if it left the box
                    left_the_box = out_of_box(ligand_coords, box_dim)

                    if not left_the_box:
                        # check for overlapps
                        overlapp = is_overlapp(vdwbox,
                                ligand_coords,
                                vdw_radii,
                                solvent_score,
                                protein_score)

                        # no overlapp
                        if not overlapp:
                            energy = score_energy(potbox, ligand_coords, charges)
                            # increase counter in counter_matrix
                            interaction_energy[x, y, z] = energy
                            counter_matrix[x, y, z] += 1.
    return interaction_energy_np, counter_matrix_np

