# cython: embedsignature=True
__author__     = "Christoph Wilms"
__copyright__  = "Copyright 2013, Epitopsy"
__credits__    = ["Christoph Wilms"]

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


@cython.cdivision(True)
@cython.boundscheck(False)
cdef double[:] transform_box_to_real_space(double[:] grid_coord,
        double[:] box_mesh_size,
        double[:] box_offset):
    '''
    Calculate the X,Y,Z coordinates in Angstroms of a particular grid point
    whose indices in the Numpy matrix are given by **grid_coord**.
    
    :math:`\\displaystyle \\text{new\\_grid\\_coord}[i] = \\text{grid\_coord}[i]
           \\cdot \\text{box\\_mesh\\_size}[i] + \\text{box\\_offset}[i]`
    
    See :func:`transform_real_to_box_space` for the reverse operation.
    
    :param grid_coord: indices of a grid point from the Numpy matrix
        :attr:`DXFile.DXBox.box`
    :type  grid_coord: list
    :param box_mesh_size: grid resolution in Angstroms, as stored in
        :attr:`DXFile.DXBox.box_mesh_size`
    :type  box_mesh_size: list
    :param box_offset: vector between the grid origin and the DXBox geometrical
        center in Angstroms, as stored in :attr:`DXFile.DXBox.box_offset`
    :type  box_offset: list
    
    :returns: (X,Y,Z) coordinates of a grid point in Angstroms.
    '''
    cdef double[:] new_grid_coord = np.zeros(3)
    cdef int dim = 3
    cdef int i
    for i in range(dim):
        new_grid_coord[i] = grid_coord[i] * box_mesh_size[i] + box_offset[i]
    
    return new_grid_coord


@cython.cdivision(True)
@cython.boundscheck(False)
cdef long[:] transform_real_to_box_space(double[:] grid_coord,
        double[:] box_mesh_size,
        double[:] box_offset):
    '''
    Do the opposite of :func:`transform_box_to_real_space`: from X,Y,Z atomic
    coordinates in Angstroms, find the closest grid point.
    
    :math:`\\displaystyle \\text{new\\_grid\\_coord}[i] = \\frac{\\text{
     grid\\_coord}[i] - \\text{box\\_offset}[i]}{\\text{box\\_mesh\\_size}[i]}`
    
    :param grid_coord: X,Y,Z atomic coordinates in Angstroms
    :type  grid_coord: list
    :param box_mesh_size: grid resolution in Angstroms, as stored in
        :attr:`DXFile.DXBox.box_mesh_size`
    :type  box_mesh_size: list
    :param box_offset: vector between the grid origin and the DXBox geometrical
        center in Angstroms, as stored in :attr:`DXFile.DXBox.box_offset`
    :type  box_offset: list
    '''
    cdef long[:] new_grid_coord = np.zeros(3,dtype=np.int)
    cdef int dim = 3
    cdef int i
    for i in range(dim):
        new_grid_coord[i] = <int> round((grid_coord[i]-box_offset[i]) / box_mesh_size[i])
    
    return new_grid_coord


@cython.cdivision(True)
@cython.boundscheck(False)
cdef bool is_overlapp(double[:,:,:] vdwbox,
        long[:,:] ligand_coords,
        long[:] vdw_radii,
        double solvent_score,
        double protein_score):
    '''
    Check if the molecular probe overlaps with the protein.
    
    Method: for every atom in the ligand, find all grid points within the van
    der Waals radius of that atom and return ``True`` if any of them is equal
    to **protein_score**
    
    :param vdwbox: van der Waals box from an APBS calculation
    :type  vdwbox: np.ndarray
    :param ligand_coords: atomic positions of a discretized molecular probe,
        as obtained from :func:`transform_ligand_to_box`
    :type  ligand_coords: np.array
    :param vdw_radii: atomic radii snapped to the grid (*i.e.* the quantity of
        grid point that can be aligned in the X,Y, and Z directions
        respectively, typically between 1 and 6)
    :type  vdw_radii: np.ndarray
    :param solvent_score: score for the solvent, usually 1.0
    :type  solvent_score: float
    :param protein_score: score for the protein, usually 0.0
    :type  protein_score: float
    
    :returns: ``True`` if there is a least one grid point overlapping,
        ``False`` otherwise
    '''
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


@cython.cdivision(True)
@cython.boundscheck(False)
cdef double score_energy(double[:,:,:] espbox,
        long[:,:] ligand_coords,
        double[:] charges):
    '''
    Multiply the atomic charge by the magnitude of the protein electrostatic
    potential *V* for each atom in the molecular probe. Sum over all atoms.
    
    :math:`E = \\sum_n^N q_n \\left(i,j,k\\right) \\cdot V \\left(i,j,k\\right)`
    
    with i,j,k the coordinates on the grid.
    
    :param espbox: electrostatic box from an APBS calculation
    :type  espbox: np.ndarray
    :param ligand_coords: atomic positions of the molecular probe, rotated and
        snapped to the grid
    :type  ligand_coords: np.ndarray
    :param charges: all atomic charges from the PQR file
    :type  charges: list
    
    :returns: The sum of the electrostatic energy between the ligand and the
        protein, adimensional.
    '''
    cdef int n_coords = ligand_coords.shape[0]
    cdef double energy = 0.
    cdef long x,y,z
    cdef int i
    for i in range(n_coords):
        x = ligand_coords[i,0]
        y = ligand_coords[i,1]
        z = ligand_coords[i,2]
        energy += espbox[x,y,z] * charges[i]
    
    return energy


@cython.cdivision(True)
@cython.boundscheck(False)
cdef double[:,:] move_ligand(double[:,:] ligand_real_coords,
        double[:] coords):
    '''
    Compute the geometrical center of the molecular probe and move it to the
    coordinates of the grid point **coord**.
    
    :param ligand_real_coords: coordinates of every atoms in the molecular probe
    :type  ligand_real_coords: list
    :param coord: coordinates of the grid point, as returned by
        :func:`transform_box_to_real_space`
    :type  coord: list
    
    :returns: The coordinates of the molecular probe whose geometrical center
        coincides with the grid point given as parameter.
    '''
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


@cython.cdivision(True)
@cython.boundscheck(False)
cdef long[:,:] transform_ligand_to_box(double[:,:] ligand_real_coords,
        double[:] box_mesh_size,
        double[:] box_offset):
    '''
    Snap each atom in the molecular probe to the closest grid points. Make use
    of :func:`transform_real_to_box_space`.
    
    :param ligand_real_coords: atomic coordinates of the molecular probe
    :type  ligand_real_coords: np.ndarray
    :param box_mesh_size: grid resolution in Angstroms
    :type  box_mesh_size: list
    :param box_offset: vector between the grid origin and the DXBox geometrical
        center in Angstroms
    :type  box_offset: list
    
    :returns: A list of grid points (discretized molecular probe).
    '''
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


@cython.cdivision(True)
@cython.boundscheck(False)
cdef bool out_of_box(long[:,:] ligand_coords, long[:] box_dim):
    '''
    Check if the molecular probe has left the box.
    
    :param ligand_coords: atomic positions of a discretized molecular probe,
        as obtained from :func:`transform_ligand_to_box`
    :type  ligand_coords: np.ndarray
    :param box_dim: number of grid points in the X,Y,Z directions
    :type  box_dim: list
    
    :returns: ``True`` if the molecular probe left the box, ``False`` otherwise
    '''
    cdef int i,j
    cdef int n_coords = ligand_coords.shape[0]
    cdef int x = 0
    
    for i in range(n_coords):
        for j in range(3):
            x = ligand_coords[i,j]
            if x < 0 or x >= box_dim[j]:
                return True
    
    return False


def explicit_sampling(ligand_pqr,
                      np.ndarray[np.float64_t, ndim=3] vdwbox_np,
                      espbox,
                      np.ndarray[np.float64_t, ndim=3] counter_matrix_np,
                      score_solvent,
                      score_protein):
    '''
    Calculate the microstates energy without FFT. Translate the ligand over
    the protein surface and integrate the energy if there is no overlap with
    the protein core. 
    
    :param ligand_pqr: the molecular probe
    :type  ligand_pqr: :class:`Structure.PQRFile`
    :param vdwbox_np: van der Waals box from an APBS calculation
    :type  vdwbox_np: np.ndarray
    :param espbox: electrostatic box from an APBS calculation
    :type  espbox: :class:`DXFile.DXBox`
    :param counter_matrix_np: list of grid points which do not result in an
        overlap between the ligand and the protein
    :type  counter_matrix_np: np.ndarray
    :param score_solvent: value encoding the solvent in **vdwbox_np**
    :type  score_solvent: float
    :param score_protein: value encoding the protein core in **vdwbox_np**
    :type  score_protein: float
    
    :returns: A tuple containing the microstates energy and an incremented
        **counter_matrix**, both as np.ndarray.
    '''
    ## general stuff
    cdef int xdim = espbox.box.shape[0]
    cdef int ydim = espbox.box.shape[1]
    cdef int zdim = espbox.box.shape[2]
    cdef int x,y,z
    cdef double energy = 0.
    cdef double solvent_score = score_solvent
    cdef double protein_score = score_protein
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

