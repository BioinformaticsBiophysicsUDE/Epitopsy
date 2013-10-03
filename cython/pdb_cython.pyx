"""
@author: Christoph Wilms
"""
import numpy as np

# "cimport" is used to import special compile-time information
# about the numpy module (this is stored in a file numpy.pxd which is
# currently part of the Cython distribution).
cimport numpy as np
cimport cython

from libc.math cimport sqrt, round, exp, pow

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


def determine_max_diameter(np.ndarray[DTYPE_float_t, ndim=2] atom_coord not None):
    # This function calculates the max diameter of a protein.
    cdef int xdim = atom_coord.shape[0]
    cdef int i,j
    cdef float maxDiameter = 0
    cdef float diam = 0
    cdef float x,y,z
    
    for i from 0 <= i < xdim:
        for j from i+1 <=j < xdim:
            x = atom_coord[i,0] - atom_coord[j,0]
            y = atom_coord[i,1] - atom_coord[j,1]
            z = atom_coord[i,2] - atom_coord[j,2]
            diam = sqrt(x*x + y*y + z*z)
            if diam > maxDiameter:
                maxDiameter = diam
    
    return maxDiameter

def build_vdw_surface(np.ndarray[DTYPE_float_t, ndim=3] vdw_array):
    # This function iterates over the grid and where it finds an element > 0
    # it will use this radius and fill all neighbors, which are 0, with -1.
    # The vdw surface is then defined by all points != 0.
    shape = vdw_array.shape
    cdef int xdim = <int> shape[0]
    cdef int ydim = <int> shape[1]
    cdef int zdim = <int> shape[2]
    vdw_array = _build_vdw_surface_fast(vdw_array, xdim, ydim, zdim)
    return vdw_array

cdef np.ndarray _build_vdw_surface_fast(np.ndarray[DTYPE_float_t, ndim=3] vdw_array_np,
                                       int xdim, int ydim, int zdim):
    # This function iterates over the grid and where it finds an element > 0
    # it will use this radius and fill all neighbors, which are 0, with -1.
    # The vdw surface is then defined by all points != 0.
    # This is 'pure' c code!
    cdef double[:,:,:] vdw_array = vdw_array_np
    cdef int x,y,z
    cdef int i,j,k
    cdef int radius
    cdef float vdw_radius
    cdef float point_distance
    for x from 0 <= x < xdim:
        for y from 0 <= y < ydim:
            for z from 0 <= z < zdim:
                # find centers of atoms
                vdw_radius = vdw_array[x,y,z]
                if vdw_radius > 0:
                    # set the current node to -1
                    vdw_array[x,y,z] = -1.
                    # iterate the environment and set all to -1 which are 0
                    radius = <int> round(vdw_radius)
                    for i from x-radius <= i <= x+radius:
                        for j from y-radius <= j <= y+radius:
                            for k from z-radius <= k <= z+radius:
                                point_distance = sqrt((x-i)*(x-i)
                                                            +(y-j)*(y-j)
                                                            +(z-k)*(z-k))
                                if point_distance <= vdw_radius:
                                    if (i >= 0 and j >= 0 and z >= 0 and
                                        i < xdim and j < ydim and k < zdim):
                                        if vdw_array[i,j,k] == 0.:
                                            vdw_array[i,j,k] = -1.

    # set everything else than 0 to 1:
    for x from 0 <= x < xdim:
        for y from 0 <= y < ydim:
            for z from 0 <= z < zdim:
                if vdw_array[x,y,z] != 0.0:
                    vdw_array[x,y,z] = 1.

    return vdw_array_np
     
cdef float calculate_rmsd(np.ndarray[DTYPE_float_t, ndim=2] coord_1,
        np.ndarray[DTYPE_float_t, ndim=2] coord_2):
    cdef int i
    cdef float rmsd = 0
    cdef float diff_x, diff_y, diff_z
    cdef int xdim = coord_1.shape[0]
    cdef float n = <float> xdim
    for i from 0 <= i < xdim:
        diff_x = coord_1[i,0] - coord_2[i,0]
        diff_y = coord_1[i,1] - coord_2[i,1]
        diff_z = coord_1[i,2] - coord_2[i,2]
        rmsd = rmsd + (diff_x * diff_x + diff_y * diff_y + diff_z * diff_z)
    
    rmsd = sqrt(rmsd / n)
    return rmsd



def get_rmsd(np.ndarray[DTYPE_float_t, ndim=2] reference_coords,
             np.ndarray[DTYPE_float_t, ndim=2] coords):
    cdef float rmsd
    # center on centroid
    av_1 = np.sum(coords,0)/coords.shape[0]
    av_2 = np.sum(reference_coords,0)/reference_coords.shape[0]
    coords_n = coords- av_1
    reference_coords_n = reference_coords - av_2
    # correlation matrix
    a = np.dot(np.transpose(coords_n), reference_coords_n)
    u, d, vt = np.linalg.svd(a)
    rot = np.transpose(np.dot(np.transpose(vt),np.transpose(u)))
    # check if we have found a reflection
    if np.linalg.det(rot)<0:
        vt[2]=-vt[2]
        rot=np.transpose(np.dot(np.transpose(vt), np.transpose(u)))
    tran=av_2-np.dot(av_1, rot)
    # transform coord_2
    new_coord = np.dot(coords, rot) + tran
    rmsd = calculate_rmsd(new_coord, reference_coords)
    return {'rmsd' : rmsd, 'rotation' : rot, 'translation' : tran}

cdef np.ndarray calculate_hydrophic_potential(
        np.ndarray[DTYPE_float_t, ndim=3] potential,
        np.ndarray[DTYPE_float_t, ndim=2] res_center_list,
        np.ndarray[DTYPE_float_t, ndim=1] res_hydro_list,
        box_meshsize, box_dim, box_offset):
    cdef int xdim = box_dim[0]
    cdef int ydim = box_dim[1]
    cdef int zdim = box_dim[2]
    cdef float local_vec[3] 
    cdef float res_vec[3]
    cdef float diff_vec[3]
    cdef float distance
    cdef int x,y,z
    cdef int n
    cdef int list_length = res_center_list.shape[0]
    cdef float hydro_charge
    cdef float hydro_potential
    cdef int vec_length = 3
    cdef float meshsize[3]
    cdef float offset[3]
    cdef int dim[3]
    for i in range(vec_length):
        meshsize[i] = <float> box_meshsize[i]
        offset[i] = <float> box_offset[i]
    
    for n in range(list_length):
        hydro_charge = res_hydro_list[n]
        for i in range(vec_length):
            res_vec[i] = res_center_list[n,i]

        for x in range(xdim):
            for y in range(ydim):
                for z in range(zdim):
                    local_vec[0] = <float> x
                    local_vec[1] = <float> y
                    local_vec[2] = <float> z
                    # transform to real space
                    for j in range(vec_length):
                        local_vec[j] = local_vec[j] * meshsize[j] + offset[j]
                    
                    for k in range(vec_length):
                        diff_vec[k] = res_vec[k] - local_vec[k]
                    distance = sqrt(diff_vec[0] * diff_vec[0]
                                    + diff_vec[1] * diff_vec[1]
                                    + diff_vec[2] * diff_vec[2])
                    
                    potential[x,y,z] = potential[x,y,z] + hydro_charge * exp(-distance)

    return potential

def get_hydrophic_potential(np.ndarray[DTYPE_float_t, ndim=2] res_center_list,
                            np.ndarray[DTYPE_float_t, ndim=1] res_hydro_list,
                            box_meshsize, box_dim, box_offset):
    
    cdef np.ndarray[DTYPE_float_t, ndim=3] potential = np.zeros(box_dim, dtype = DTYPE_float)

    potential = calculate_hydrophic_potential(potential,
            res_center_list, res_hydro_list, box_meshsize, box_dim,
            box_offset)

    return potential


def find_clashing_atoms(np.ndarray[DTYPE_float_t, ndim=2] coord_1,
        np.ndarray[DTYPE_float_t, ndim=1] radii_1,
        np.ndarray[DTYPE_float_t, ndim=2] coord_2,
        np.ndarray[DTYPE_float_t, ndim=1] radii_2,
        float energy_cutoff):
    return _find_clashing_atoms(coord_1, radii_1, coord_2, radii_2,
            energy_cutoff)

cdef int _find_clashing_atoms(np.ndarray[DTYPE_float_t, ndim=2] coord_1,
        np.ndarray[DTYPE_float_t, ndim=1] radii_1,
        np.ndarray[DTYPE_float_t, ndim=2] coord_2,
        np.ndarray[DTYPE_float_t, ndim=1] radii_2,
        float energy_cutoff):
    cdef int n_atoms_1 = coord_1.shape[0]
    cdef int n_atoms_2 = coord_2.shape[0]

    cdef float r_1, r_2, r_vdw, division
    cdef float x[3],u[3],dist, energy
    cdef float e_const = 1.
    cdef int i,j,k

    cdef int atom_clashes = 0

    for i in range(n_atoms_1):
        r_1 = radii_1[i]
        for k in range(3):
            x[k] = coord_1[i,k]
        for j in range(n_atoms_2):
            r_2 = radii_2[j]
            for k in range(3):
                u[k] = coord_2[j,k]
            
            # E = e_const * [ (r_vdw_ij / dist_ij)^12 - 2 * (r_vdw_ij / dist_ij)^6] - 0.008175222784000003 * e
            # with r_vdw_ij = (r_i / 2.) + (r_j / 2.)
            # and e_const = 1, but can also be included in each ()^x as 
            # e* = sqrt(e_i * e_j)
            # cutoff at r_ij < r_c = 2.5 * r_vdw leads to an
            # additional term for the potential, couse otherwise
            # there would be a jump in the potential
            dist = 0
            for k in range(3):
                dist = dist + pow((x[k] - u[k]), 2.)
            dist = sqrt(dist)
            
            r_vdw = r_1 / 2. + r_2 / 2.

            if dist == 0.:
                atom_clashes += 1
                break
            elif dist > 2.5 * r_vdw:
                # we are in the land of the cutoff
                continue
            else:
                division = r_vdw / dist
                energy = e_const * (pow(division, 12.) - 2. * pow(division, 6.)) - 0.008175222784000003 * e_const
                if energy > energy_cutoff:
                    atom_clashes += 1
                    break

    return atom_clashes

