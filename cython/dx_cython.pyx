"""
@author: Christoph Wilms
"""
import numpy as np

# "cimport" is used to import special compile-time information
# about the numpy module (this is stored in a file numpy.pxd which is
# currently part of the Cython distribution).
cimport numpy as np

from libc.math cimport sqrt, round, ceil

import gzip
import re
import os

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


cdef float _vec_distance(int x1, int y1, int z1, int x2, int y2, int z2):
    return sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1))

def read_dxfile_depreciated(filename, box_type):
    if not(os.path.isabs(filename)):
        filename = os.path.join(os.getcwd(), filename)

    cdef int countX = 0
    cdef int countY = 0
    cdef int countZ = 0
    cdef long counter = 0
    cdef np.ndarray[np.int64_t,ndim=1] dimension = np.zeros(3,dtype=np.int64)
    cdef np.ndarray[np.float64_t,ndim=1] offset = np.zeros(3,dtype=np.float64)
    cdef np.ndarray[np.float64_t,ndim=1] meshsize = np.zeros(3,dtype=np.float64)
    dimensionPattern = re.compile('(\\d+) (\\d+) (\\d+)$')
    lastHeaderLinePattern = re.compile('^object 3')
    objectPattern = re.compile('^object 1')
    valuePattern = re.compile('^[-]?\\d+\\.\\d+[eE][+-]\\d+')
    # this has been tested and it works
    originPattern = re.compile('^origin ([-]?\\d+.\\d+[eE]?[+-]?\\d+?) ([-]?\\d+.\\d+[eE]?[-+]?\\d+?) ([-]?\\d+.\\d+[eE]?[+-]?\\d+?)')

    if filename.endswith('.gz'):
        with gzip.open(filename) as infile:
            content = infile.readlines()
    else:
        with open(filename) as infile:
            content = infile.readlines()

    cdef int delta_counter = 0
    cdef np.ndarray[np.float64_t,ndim=1] box
    for line in content:
        if not line.startswith('#'):
            if line.startswith('object 1'):
                matchedObject = dimensionPattern.search(line)
                dimension[0] = int(matchedObject.group(1))
                dimension[1] = int(matchedObject.group(2))
                dimension[2] = int(matchedObject.group(3))
                box = np.zeros(dimension[0] * dimension[1] * dimension[2],dtype=np.float64)
            elif line.startswith('origin'):
                matchedOrigin = line.split(' ')
                offset[0] = float(matchedOrigin[1])
                offset[1] = float(matchedOrigin[2])
                offset[2] = float(matchedOrigin[3])
            elif line.startswith('delta'):
                    delta_split = line.split()
                    meshsize[delta_counter] = float(delta_split[delta_counter+1])
                    delta_counter += 1

            if valuePattern.match(line):
                numbers_list = line.strip().split(' ')
                for item in numbers_list:
                    box[counter] = float(item)
                    counter += 1

    cdef np.ndarray[np.float64_t,ndim=3] box_3d = box.reshape(dimension)

    return [box_3d, meshsize, offset]

def read_dxfile(filename, box_type):
    if not(os.path.isabs(filename)):
        filename = os.path.join(os.getcwd(), filename)

    cdef int countX = 0
    cdef int countY = 0
    cdef int countZ = 0
    cdef long counter = 0
    cdef np.ndarray[np.int64_t,ndim=1] dimension = np.zeros(3,dtype=np.int64)
    cdef np.ndarray[np.float64_t,ndim=1] offset = np.zeros(3,dtype=np.float64)
    cdef np.ndarray[np.float64_t,ndim=1] meshsize = np.zeros(3,dtype=np.float64)
    dimensionPattern = re.compile('(\\d+) (\\d+) (\\d+)$')
    lastHeaderLinePattern = re.compile('^object 3')
    objectPattern = re.compile('^object 1')
    valuePattern = re.compile('^[-]?\\d+\\.\\d+[eE][+-]\\d+')
    # this has been tested and it works
    originPattern = re.compile('^origin ([-]?\\d+.\\d+[eE]?[+-]?\\d+?) ([-]?\\d+.\\d+[eE]?[-+]?\\d+?) ([-]?\\d+.\\d+[eE]?[+-]?\\d+?)')

    if filename.endswith('.gz'):
        with gzip.open(filename) as infile:
            content = infile.readlines()
    else:
        with open(filename) as infile:
            content = infile.readlines()

    cdef int delta_counter = 0
    cdef np.ndarray[np.float64_t,ndim=1] box
    cdef long i_start
    for i_start, line in enumerate(content):
        if not line.startswith('#'):
            if line.startswith('object 1'):
                matchedObject = dimensionPattern.search(line)
                dimension[0] = int(matchedObject.group(1))
                dimension[1] = int(matchedObject.group(2))
                dimension[2] = int(matchedObject.group(3))
                box = np.zeros(dimension[0] * dimension[1] * dimension[2],dtype=np.float64)
            elif line.startswith('origin'):
                matchedOrigin = line.split(' ')
                offset[0] = float(matchedOrigin[1])
                offset[1] = float(matchedOrigin[2])
                offset[2] = float(matchedOrigin[3])
            elif line.startswith('delta'):
                    delta_split = line.split()
                    meshsize[delta_counter] = float(delta_split[delta_counter+1])
                    delta_counter += 1

            if valuePattern.match(line):
                break

    cdef long i,j
    cdef long n_elements = dimension[0] * dimension[1] * dimension[2]
    cdef long i_end = 1 + i_start + <long> ceil( (<double> n_elements) / 3.)
    for i in range(i_start, i_end):
        line_content = content[i].split(' ')
        for j in range(3):
            if counter < n_elements:
                box[counter] = float(line_content[j])
                counter += 1

    cdef np.ndarray[np.float64_t,ndim=3] box_3d = box.reshape(dimension)

    return [box_3d, meshsize, offset]

def write_box(outfile, np.ndarray[np.float64_t,ndim=3]box_np,int values_per_line):
    cdef int value_count = 0
    cdef int x,y,z
    cdef int xdim = box_np.shape[0]
    cdef int ydim = box_np.shape[1]
    cdef int zdim = box_np.shape[2]
    cdef double[:,:,:] box = box_np
    for x in range(xdim):
        for y in range(ydim):
            for z in range(zdim):
                outfile.write("{0:e} ".format(box[x, y, z]))
                value_count += 1
                if value_count == values_per_line:
                    outfile.write('\n')
                    value_count = 0
    return None

def flood(np.ndarray[np.float64_t, ndim=3] box_np not None,
        np.ndarray[np.int64_t, ndim=2] neighbors):
    """
    The VDW grid is a grid of 0's and 1's, like this:
            1 1 1 1
            1 0 0 1
            1 0 1 1
            1 1 1 1
    The algorithm starts at [0,0,0] and appends the neighbors, which
    have a value equal to a solvent grid node (i.e. 1), to the
    newFront-list. After all neighbor nodes from the first list have
    been checked, the list gets deleted and contains now the new
    neighbors.
    The process is repeated until now new neighbors are found, that is
    there are no new neighbors.

    Then the remaining solvent values are treated as protein nodes
    (i.e. they are set to a peptide score, which is 0). All nodes
    that have been assigned a temp value are then set back to 1.
    """
    cdef double[:,:,:] box = box_np
    cdef int xdim = box.shape[0]
    cdef int ydim = box.shape[1]
    cdef int zdim = box.shape[2]
    cdef int i,j
    cdef int n_iter
    cdef int nx,ny,nz
    cdef int neighbor_elements = neighbors.shape[0]
    cdef int x,y,z
    cdef int counter = 0
    cdef int posx, posy, posz
    cdef np.ndarray[np.int64_t,ndim=2] waterFront = np.zeros([xdim*ydim*zdim,3],dtype = np.int64)
    cdef np.ndarray[np.int64_t,ndim=2] newFront = np.zeros([xdim*ydim*zdim,3],dtype = np.int64)
    cdef int water_counter = 0
    cdef int front_counter = 0

    cdef np.float64_t temp_value = -1
    cdef np.float64_t solventScore = 1
    cdef np.float64_t peptideScore = 0

    # reset counter and array
    waterFront[water_counter,0] = 0
    waterFront[water_counter,1] = 0
    waterFront[water_counter,2] = 0
    water_counter = 1
    while(water_counter > 0):
        if front_counter > 0:
            front_counter = 0
        for j from 0 <= j < water_counter:
            posx = waterFront[j,0]
            posy = waterFront[j,1]
            posz = waterFront[j,2]
            box[posx, posy,posz] = temp_value
            for n_iter from 0 <= n_iter < neighbor_elements:
                nx = posx + neighbors[n_iter,0]
                ny = posy + neighbors[n_iter,1]
                nz = posz + neighbors[n_iter,2]
                if (nx >= 0 and nx < xdim and ny >= 0 and ny < ydim and
                    nz >= 0 and nz < zdim):
                    if box[nx,ny,nz] == solventScore:
                        box[nx, ny, nz] = temp_value
                        newFront[front_counter,0] = nx
                        newFront[front_counter,1] = ny
                        newFront[front_counter,2] = nz
                        front_counter = front_counter + 1
        if water_counter > 0:
            water_counter = front_counter
            waterFront[0:water_counter,:] = newFront[0:water_counter,:]

    for x in range(xdim):
        for y in range(ydim):
            for z in range(zdim):
                if box[x,y,z] == temp_value:
                    # this should be the solvent
                    box[x,y,z] = solventScore
                else:
                    # this is all protein
                    box[x,y,z] = peptideScore

    return box_np


def floodxy(np.ndarray[np.float64_t, ndim=3] box_np not None,
        np.ndarray[np.int64_t, ndim=2] neighbors):
    """
    The VDW grid is a grid of 0's and 1's, like this:
            1 1 1 1
            1 0 0 1
            1 0 1 1
            1 1 1 1
    The algorithm starts at [0,0,0] and appends the neighbors, which
    have a value equal to a solvent grid node (i.e. 1), to the
    newFront-list. After all neighbor nodes from the first list have
    been checked, the list gets deleted and contains now the new
    neighbors.
    The process is repeated until now new neighbors are found, that is
    there are no new neighbors.

    Then the remaining solvent values are treated as protein nodes
    (i.e. they are set to a peptide score, which is 0). All nodes
    that have been assigned a temp value are then set back to 1.
    """
    cdef double[:,:,:] box = box_np
    cdef int xdim = box.shape[0]
    cdef int ydim = box.shape[1]
    cdef int zdim = box.shape[2]
    cdef int i,j
    cdef int n_iter
    cdef int nx,ny,nz
    cdef int neighbor_elements = neighbors.shape[0]
    cdef int x,y,z
    cdef int counter = 0
    cdef int posx, posy, posz
    cdef np.ndarray[np.int64_t,ndim=2] waterFront = np.zeros([xdim*ydim,2],dtype = np.int64)
    cdef np.ndarray[np.int64_t,ndim=2] newFront = np.zeros([xdim*ydim,2],dtype = np.int64)
    cdef int water_counter = 0
    cdef int front_counter = 0

    cdef np.float64_t temp_value = -1
    cdef np.float64_t solventScore = 1
    cdef np.float64_t peptideScore = 0

    # reset counter and array
    for posz in range(zdim):
        waterFront[water_counter,0] = 0
        waterFront[water_counter,1] = 0
        water_counter = 1
        while(water_counter > 0):
            if front_counter > 0:
                front_counter = 0
            for j from 0 <= j < water_counter:
                posx = waterFront[j,0]
                posy = waterFront[j,1]
                box[posx, posy,posz] = temp_value
                for n_iter from 0 <= n_iter < neighbor_elements:
                    nx = posx + neighbors[n_iter,0]
                    ny = posy + neighbors[n_iter,1]
                    nz = posz
                    if (nx >= 0 and nx < xdim and ny >= 0 and ny < ydim):
                        if box[nx,ny,nz] == solventScore:
                            box[nx, ny, nz] = temp_value
                            newFront[front_counter,0] = nx
                            newFront[front_counter,1] = ny
                            front_counter = front_counter + 1
            if water_counter > 0:
                water_counter = front_counter
                waterFront[0:water_counter,:] = newFront[0:water_counter,:]

    for x in range(xdim):
        for y in range(ydim):
            for z in range(zdim):
                if box[x,y,z] == temp_value:
                    # this should be the solvent
                    box[x,y,z] = solventScore
                else:
                    # this is all protein
                    box[x,y,z] = peptideScore

    return box_np


def get_sas(np.ndarray[np.float64_t, ndim=3] box,
        protein_score, solvent_score, probe_radius):
    return _get_sas(box, <float> protein_score, <float> solvent_score,
                    <float> probe_radius)

cdef np.ndarray _get_sas(np.ndarray[np.float64_t, ndim=3] box,
        float protein_score, float solvent_score, float probe_radius):
    cdef int xdim = box.shape[0]
    cdef int ydim = box.shape[1]
    cdef int zdim = box.shape[2]
    cdef list neighbors_list = []
    cdef int i,j,k, n
    cdef int x,y,z
    cdef int nx,ny,nz
    cdef int probe_int = int(round(probe_radius))
    cdef float temp_score = -10
    cdef float surface_score = -11
    for i from -probe_int <= i <= probe_int:
        for j from -probe_int <= j <= probe_int:
            for k from -probe_int <= k <= probe_int:
                if i == 0 and j == 0 and k ==0:
                    continue
                else:
                    if _vec_distance(i,j,k,0,0,0) <= probe_radius:
                        neighbors_list.append([i,j,k])

    cdef np.ndarray[np.int64_t, ndim=2] neighbors = np.array(neighbors_list, dtype=np.int64)
    cdef int neighbor_elements = neighbors.shape[0]

    for x in range(xdim):
        for y in range(ydim):
            for z in range(zdim):
                if box[x,y,z] == solvent_score:
                    # assume that this is a valid solvent grid point
                    box[x,y,z] = temp_score
                    # check neighbors, if they are also solvent
                    for n in range(neighbor_elements):
                        nx = x + neighbors[n,0]
                        ny = y + neighbors[n,1]
                        nz = z + neighbors[n,2]
                        if (nx >= 0 and nx < xdim
                            and ny >= 0 and ny < ydim
                            and nz >= 0 and nz < zdim):
                            if box[nx,ny,nz] == protein_score:
                                # point of x,y,z is to close to the protein
                                box[x,y,z] = surface_score
                                break

    for x in range(xdim):
        for y in range(ydim):
            for z in range(zdim):
                if box[x,y,z] == temp_score:
                    # this should be the solvent
                    box[x,y,z] = solvent_score
                else:
                    # this is all protein
                    box[x,y,z] = protein_score

    return box


def determine_actual_surface_nodes(np.ndarray[np.float64_t, ndim=3] box,
        score, np.ndarray[np.int64_t, ndim = 2] neighbors):
    surface_candidates = np.nonzero(box != score)
    surface_nodes = []
    cdef float solventScore = score
    cdef int neighbor_elements = len(neighbors)
    cdef int i
    cdef int j
    cdef int nx, ny, nz
    cdef int posx, posy, posz
    cdef int xdim = box.shape[0]
    cdef int ydim = box.shape[1]
    cdef int zdim = box.shape[2]
    for i from 0 <= i < len(surface_candidates[:][0]) :
        for j from 0 <= j < neighbor_elements:
            nx = neighbors[j,0]
            ny = neighbors[j,1]
            nz = neighbors[j,2]
            posx = surface_candidates[0][i]+nx
            posy = surface_candidates[1][i]+ny
            posz = surface_candidates[2][i]+nz
            if (posx >= 0 and posx < xdim
                and posy >= 0 and posy < ydim
                and posz >= 0 and posz < zdim):
                if box[posx, posy, posz] == solventScore:
                    surface_nodes.append([surface_candidates[0][i],
                                        surface_candidates[1][i],
                                        surface_candidates[2][i]])
                    break
    return surface_nodes


def determineActualInnerPeptideNodes(np.ndarray[np.float64_t, ndim=3] box,
        solventScore_ext, peptideScore_ext,
        np.ndarray[np.int64_t, ndim = 2] neighbors):
    innerNodes = []
    cdef int peptideScore = <int> peptideScore_ext
    cdef int solventScore = <int> solventScore_ext
    cdef int neighbor_elements = len(neighbors)
    cdef int xdim = box.shape[0]
    cdef int ydim = box.shape[1]
    cdef int zdim = box.shape[2]
    cdef int x,y,z
    cdef int i
    cdef int nx, ny, nz
    for x from 0 <= x < xdim:
        for y from 0 <= y < ydim:
            for z from 0 <= z < zdim:
                if box[x, y, z] == peptideScore:
                    innerNode = True
                    # there should be no boundary error, otherwise there
                    # is something wrong!
                    for i from 0 <= i < neighbor_elements:
                        nx = neighbors[i, 0]
                        ny = neighbors[i, 1]
                        nz = neighbors[i, 2]
                        if box[nx + x, ny + y, nz + z] == solventScore:
                            innerNode = False
                            break
                    if innerNode:
                        innerNodes.append([x, y, z])
    return innerNodes

def extendSurface(int iterations, np.ndarray[np.float64_t, ndim=3] box,
        solventScore_ext, peptideScore_ext,
        np.ndarray[np.int64_t, ndim = 2] neighbors):
    cdef int i, s, n
    cdef int neighbor_elements = len(neighbors)
    cdef int nx, ny, nz
    cdef int xdim = box.shape[0]
    cdef int ydim = box.shape[1]
    cdef int zdim = box.shape[2]
    cdef float solventScore = <float> solventScore_ext
    cdef float peptideScore = <float> peptideScore_ext
    cdef np.ndarray[np.int64_t, ndim = 2] surface_nodes
    cdef np.ndarray[np.int64_t, ndim = 2] next_neighbors = np.zeros([6, 3], dtype = np.int64)
    # [-1, 0, 0]
    next_neighbors[0, 0] = -1
    # [0, -1, 0]
    next_neighbors[1, 1] = -1
    # [0, 0, -1]
    next_neighbors[2, 2] = -1
    # [1, 0, 0]
    next_neighbors[3, 0] = 1
    # [0, 1, 0]
    next_neighbors[4, 1] = 1
    # [0, 0, 1]
    next_neighbors[5, 2] = 1
    if iterations>0:
        surface_nodes = np.array(determine_actual_surface_nodes(box, solventScore_ext, neighbors),dtype=np.int64)
        for i from 0<= i < iterations:
            newSurfaceNodes = []
            for s from 0<= s < len(surface_nodes):
                # check only next neighbors
                for n from 0 <= n < 6:
                    nx = surface_nodes[s,0] + next_neighbors[n,0]
                    ny = surface_nodes[s,1] + next_neighbors[n,1]
                    nz = surface_nodes[s,2] + next_neighbors[n,2]
                    if (nx >= 0 and nx < xdim and ny >= 0 and ny < ydim and
                        nz >= 0 and nz < zdim):
                        if box[nx, ny, nz] == solventScore:
                                box[nx, ny, nz] = peptideScore
                                newSurfaceNodes.append([nx, ny, nz])
            if newSurfaceNodes == []:
                newSurfaceNodes.append([])
            surface_nodes = np.array(newSurfaceNodes,dtype=np.int64)
    return box

cdef np.ndarray switch_value(np.ndarray[np.float64_t, ndim =3] box,
                                          int xdim, int ydim, int zdim,
                                          np.float64_t value_to_switch,
                                          np.float64_t new_value):
    for x from 0 <= x < xdim:
        for y from 0 <= y < ydim:
            for z from 0 <= z < zdim:
                if box[x,y,z] == value_to_switch:
                    box[x,y,z] = new_value

    return box

def extend_nice_surface(int iterations, np.ndarray[np.float64_t, ndim=3] box,
        solventScore_ext, peptideScore_ext, dummy_score_ext,
        np.ndarray[np.int64_t, ndim = 2] neighbors):
    cdef int i,n,x,y,z
    cdef int nx, ny, nz
    cdef int xdim = box.shape[0]
    cdef int ydim = box.shape[1]
    cdef int zdim = box.shape[2]
    cdef float solventScore = <float> solventScore_ext
    cdef float peptideScore = <float> peptideScore_ext
    cdef float dummy_score = <float> dummy_score_ext
    cdef np.ndarray[np.int64_t, ndim = 2] next_neighbors = np.zeros([6, 3], dtype = np.int64)
    # [-1, 0, 0]
    next_neighbors[0, 0] = -1
    # [0, -1, 0]
    next_neighbors[1, 1] = -1
    # [0, 0, -1]
    next_neighbors[2, 2] = -1
    # [1, 0, 0]
    next_neighbors[3, 0] = 1
    # [0, 1, 0]
    next_neighbors[4, 1] = 1
    # [0, 0, 1]
    next_neighbors[5, 2] = 1

    for i from 0 <= i < iterations:
        for x from 0 <= x < xdim:
            for y from 0 <= y < ydim:
                for z from 0 <= z < zdim:
                    if box[x,y,z] == solventScore:
                        for n from 0 <= n < 6:
                            nx = x + next_neighbors[n,0]
                            ny = y + next_neighbors[n,1]
                            nz = z + next_neighbors[n,2]
                            if (nx >= 0 and nx < xdim and ny >= 0 and ny < ydim and
                                nz >= 0 and nz < zdim):
                                if box[nx,ny,nz] == peptideScore:
                                    box[x,y,z] = dummy_score
        box = switch_value(box, xdim, ydim, zdim, dummy_score, peptideScore)

    return box


def find_protein_surface(protein_score,
                 solvent_score,
                 surface_score,
                 np.ndarray[np.float64_t, ndim=3] box):
    return _find_protein_surface(<float> protein_score, <float> solvent_score,
                         <float> surface_score, box)

cdef double[:,:,:] check_protein_neighbor(double[:,:,:] box, int xdim, int ydim,
        int zdim, int x, int y, int z, float solvent_score,
        float surface_score):
    cdef int i,j,k
    cdef int posx,posy,posz
    for i from -1 <= i <= 1:
        for j from -1 <= j <= 1:
            for k from -1 <= k <= 1:
                if i == 0 and j == 0 and k == 0:
                    continue
                else:
                    posx = x + i
                    posy = y + j
                    posz = z + k

                    if (posx >= 0 and posx < xdim
                        and posy >= 0 and posy < ydim
                        and posz >= 0 and posz < zdim):
                        if box[posx,posy,posz] == solvent_score:
                            box[x,y,z] = surface_score
                            return box

    return box

cdef np.ndarray _find_protein_surface(float protein_score, float solvent_score,
                  float surface_score,
                  np.ndarray[np.float64_t, ndim=3] box_np):
    cdef int x,y,z,i,j,k
    cdef double[:,:,:] box = box_np
    cdef int posx, poxy, posz
    cdef int xdim = box.shape[0]
    cdef int ydim = box.shape[1]
    cdef int zdim = box.shape[2]
    # use stop to interupt the loop over the neighbors
    cdef int stop = 0
    cdef float score

    for x in range(xdim):
        for y in range(ydim):
            for z in range(zdim):
                # only check points from the proteins interior
                if box[x,y,z] == protein_score:
                    box = check_protein_neighbor(box, xdim, ydim, zdim, x, y, z,
                            solvent_score, surface_score)
#                    stop = 0
#                    # check neighbor grid points
#                    for i from -1 <= i <= 1:
#                        if stop == 1:
#                            break
#                        for j from -1 <= j <= 1:
#                            if stop == 1:
#                                break
#                            for k from -1 <= k <= 1:
#                                if i == 0 and j == 0 and k == 0:
#                                    continue
#                                else:
#                                    posx = x + i
#                                    posy = y + j
#                                    posz = z + k
#
#                                    if (posx >= 0 and posx < xdim
#                                        and posy >= 0 and posy < ydim
#                                        and posz >= 0 and posz < zdim):
#                                        if box[posx,posy,posz] == solvent_score:
#                                            box[x,y,z] = surface_score
#                                            stop = 1
#                                            break
    return box_np


def find_solvent_surface(protein_score,
                 solvent_score,
                 surface_score,
                 np.ndarray[np.float64_t, ndim=3] box):
    return _find_solvent_surface(<float> protein_score, <float> solvent_score,
                         <float> surface_score, box)

cdef double[:,:,:] check_surface_neighbor(double[:,:,:] box, int xdim, int ydim,
        int zdim, int x, int y, int z, float protein_score,
        float surface_score):
    cdef int i,j,k
    cdef int posx,posy,posz
    for i from -1 <= i <= 1:
        for j from -1 <= j <= 1:
            for k from -1 <= k <= 1:
                if i == 0 and j == 0 and k == 0:
                    continue
                else:
                    posx = x + i
                    posy = y + j
                    posz = z + k

                    if (posx >= 0 and posx < xdim
                        and posy >= 0 and posy < ydim
                        and posz >= 0 and posz < zdim):
                        if box[posx,posy,posz] == protein_score:
                            box[x,y,z] = surface_score
                            return box

    return box

cdef np.ndarray _find_solvent_surface(float protein_score, float solvent_score,
                  float surface_score,
                  np.ndarray[np.float64_t, ndim=3] box_np):
    cdef int x,y,z,i,j,k
    cdef double[:,:,:] box = box_np
    cdef int posx, poxy, posz
    cdef int xdim = box.shape[0]
    cdef int ydim = box.shape[1]
    cdef int zdim = box.shape[2]
    # use stop to interupt the loop over the neighbors
    cdef int stop = 0
    cdef float score

    for x in range(xdim):
        for y in range(ydim):
            for z in range(zdim):
                # only check points from the proteins interior
                if box[x,y,z] == solvent_score:
                    box = check_surface_neighbor(box, xdim, ydim, zdim, x, y, z,
                            protein_score, surface_score)

    return box_np