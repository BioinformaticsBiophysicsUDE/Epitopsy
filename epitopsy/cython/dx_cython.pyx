# @author: Christoph Wilms

from epitopsy.tools.MathTools import spherical_lattice_default
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
    
    # load file
    if filename.endswith('.gz'):
        with gzip.open(filename) as infile:
            content = infile.readlines()
    else:
        with open(filename) as infile:
            content = infile.readlines()
    
    # read metadata
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
    
    # read grid values
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


def flood_fill(np.ndarray[np.float64_t, ndim=3] box_np not None,
              seed_faces = True):
    '''
    A 6-way flood-fill algorithm using a queue for 3D volumetric maps.
    
    Example:
    
    .. code-block:: none
        
        initial      seeding      flooding     filling
        1 1 1 1 1    2 1 1 1 2    2 2 2 2 2    1 1 1 1 1
        1 0 0 1 1    1 0 0 1 1    2 0 0 2 2    1 0 0 1 1
        1 0 1 0 1    1 0 1 0 1    2 0 1 0 2    1 0 0 0 1
        1 0 0 0 1 => 1 0 0 0 1 => 2 0 0 0 2 => 1 0 0 0 1
        1 0 0 1 1    1 0 0 1 1    2 0 0 2 2    1 0 0 1 1
        1 0 0 0 1    1 0 0 0 1    2 0 0 0 2    1 0 0 0 1
        1 1 1 1 1    2 1 1 1 2    2 2 2 2 2    1 1 1 1 1
    
    where 1's represent the solvent and 0's the protein interior.
    
    The algorithm starts by seeding the box 8 corners (or the 6 faces) if they
    are not occupied by the protein, resulting in N seeds. The algorithm walks
    through these N seeds and checks their 6 neighbors; if they contain
    solvent, they are added to the seeds queue, resulting in M new seeds. The
    algorithm repeats the process for these M seeds until no new seed is found.
    Solvent voxels that could not be flooded are considered inaccessible and
    filled with the protein score.
    '''
    cdef np.float64_t proteinScore = 0
    cdef np.float64_t solventScore = 1
    cdef np.float64_t floodedScore = 2
    cdef double[:,:,:] box = box_np
    cdef int xdim = box.shape[0]
    cdef int ydim = box.shape[1]
    cdef int zdim = box.shape[2]
    cdef int i, j
    cdef int x, y, z    # seed coordinates
    cdef int nx, ny, nz # neighbor coordinates
    
    # shift vectors for 6-way neighbors in 3D space
    cdef np.ndarray[np.int64_t,ndim=2] neighbors = np.array([[+1, 0, 0],
                                                             [ 0,+1, 0],
                                                             [ 0, 0,+1],
                                                             [-1, 0, 0],
                                                             [ 0,-1, 0],
                                                             [ 0, 0,-1]])
    cdef int neighbor_size = neighbors.shape[0]
    
    # seed queue
    cdef np.ndarray[np.int64_t,ndim=2] seed_queue = np.zeros(
                                      [xdim * ydim * zdim, 3], dtype=np.int64)
    cdef int seed_queue_index_start = 0
    cdef int seed_queue_index_end = 0
    cdef int seed_queue_index_end_backtrack = 0
    
    if seed_faces:
        # find seeds in the 6 faces
        for face in [1, 2, 3]:
            if face == 1:
                range_x = [0, xdim-1]
            else:
                range_x = range(xdim)
            if face == 2:
                range_y = [0, ydim-1]
            else:
                range_y = range(ydim)
            if face == 3:
                range_z = [0, zdim-1]
            else:
                range_z = range(zdim)
            for x in range_x:
                for y in range_y:
                    for z in range_z:
                        if box[x,y,z] == 1:
                            box[x,y,z] = floodedScore
                            seed_queue[seed_queue_index_end, 0] = x
                            seed_queue[seed_queue_index_end, 1] = y
                            seed_queue[seed_queue_index_end, 2] = z
                            seed_queue_index_end += 1
        if seed_queue_index_end == 0:
            raise ValueError('Solvent box faces are occupied by the protein')
    else:
        # find seeds in the 8 corners
        for x,y,z in [(0,0,0), (xdim-1,0,0), (0,ydim-1,0), (0,0,zdim-1),
                      (xdim-1,ydim-1,0), (xdim-1,0,zdim-1), (0,ydim-1,zdim-1),
                      (xdim-1,ydim-1,zdim-1)]:
            if box[x,y,z] == 1:
                box[x,y,z] = floodedScore
                seed_queue[seed_queue_index_end,:] = (x,y,z)
                seed_queue_index_end += 1
        if seed_queue_index_end == 0:
            raise ValueError('Solvent box corners are occupied by the protein')
    
    # flood box (assign floodedScore to solvent voxels accessible from seeds)
    while True:
        seed_queue_index_end_backtrack = seed_queue_index_end
        # scan seeds found in the last iteration
        for i in range(seed_queue_index_start, seed_queue_index_end):
            x = seed_queue[i, 0]
            y = seed_queue[i, 1]
            z = seed_queue[i, 2]
            # scan 6-way neighbors
            for j in range(neighbors.shape[0]):
                # position of neighbor
                nx = x + neighbors[j, 0]
                ny = y + neighbors[j, 1]
                nz = z + neighbors[j, 2]
                # skip if neighbor outside box boundaries
                # if neighbor is a solvent, add to the seed queue
                if (nx >= 0 and nx < xdim and ny >= 0 and ny < ydim and
                    nz >= 0 and nz < zdim):
                    if box[nx,ny,nz] == solventScore:
                        box[nx, ny, nz] = floodedScore
                        seed_queue[seed_queue_index_end, 0] = nx
                        seed_queue[seed_queue_index_end, 1] = ny
                        seed_queue[seed_queue_index_end, 2] = nz
                        seed_queue_index_end += 1
        
        # exit loop if no neighbor remains
        if seed_queue_index_end_backtrack == seed_queue_index_end:
            break
        else:
            seed_queue_index_start = seed_queue_index_end_backtrack
    
    # fill box (assign proteinScore to solvent voxels which were not flooded)
    for x in range(xdim):
        for y in range(ydim):
            for z in range(zdim):
                if box[x,y,z] == floodedScore:
                    box[x,y,z] = solventScore
                else:
                    box[x,y,z] = proteinScore
    
    return box_np


def flood(np.ndarray[np.float64_t, ndim=3] box_np not None,
        np.ndarray[np.int64_t, ndim=2] neighbors):
    """
    The VDW grid is a grid of 0's and 1's, like this:
    
    .. code-block:: none
    
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
    
    .. code-block:: none
    
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
                            if (nx >= 0 and nx < xdim and ny >= 0 and
                                ny < ydim and nz >= 0 and nz < zdim):
                                if box[nx,ny,nz] == peptideScore:
                                    box[x,y,z] = dummy_score
        box = switch_value(box, xdim, ydim, zdim, dummy_score, peptideScore)

    return box


def find_protein_surface(protein_score, solvent_score, surface_score,
                         np.ndarray[np.float64_t, ndim=3] box):
    """
    Find the layer below the solvent accessible surface and assign
    :attr:`score_of_surface` to its grid points. This layer is
    inaccessible to the solvent and the ligands.
    The grid points in the solvent and protein interior
    are assigned :attr:`solventScore` resp. :attr:`peptideScore`.
    
    Cython wrapper for C function :func:`_find_protein_surface`.
    """
    return _find_protein_surface(<float> protein_score, <float> solvent_score,
                                 <float> surface_score, box)


cdef double[:,:,:] check_protein_neighbor(double[:,:,:] box, int xdim, int ydim,
       int zdim, int x, int y, int z, float solvent_score, float surface_score):
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
    return box_np


def find_solvent_surface(protein_score, solvent_score, surface_score,
                         np.ndarray[np.float64_t, ndim=3] box):
    """
    Find the solvent accessible surface and assign **surface_score**
    to its grid points. The grid points in the solvent and protein interior
    are assigned **solvent_score** resp. **protein_score**.
    
    Cython wrapper for C function :func:`_find_solvent_surface`.
    """
    return _find_solvent_surface(<float> protein_score, <float> solvent_score,
                                 <float> surface_score, box)


cdef double[:,:,:] check_surface_neighbor(double[:,:,:] box, int xdim, int ydim,
       int zdim, int x, int y, int z, float protein_score, float surface_score):
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
                  float surface_score, np.ndarray[np.float64_t, ndim=3] box_np):
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


def simulated_annealing(box, runs=50, Tini=1.0, gamma=0.25, block_size=100,
                        convergence_T=1e-4, convergence_rejections=None,
                        step_size=1, adaptative_step_size=True,
                        acceptance_ratio_upper=70, acceptance_ratio_lower=30,
                        seed=None, return_full_chains=False):
    '''
    Local minima search by simulated annealing.
    
    From several starting positions uniformly distributed on a sphere centered
    around the protein, run a Markov chain with a temperature schedule. The
    probability of transition from the current state to the next step is 1 when
    the next state has lower energy, or is proportionnal to the exponential of
    the energy difference divided by the temperature otherwise. The temperature
    decreases steadily until the Markov chain is unable to progress anymore.
    
    More formally, the transition probability :math:`P(x_1 \\to w_2, T)`
    between two states :math:`x_1` and :math:`x_2` with energy difference
    :math:`\\Delta E = E_2 - E_1` at temperature T is given by the Metropolis
    criterion:
    
    .. math::
    
         P(x_1 \\to x_2, T) = \\begin{cases}
                                    e^{-\\Delta E/T} & \\Delta E > 0 \\\\
                                    1                & \\Delta E \\leq 0
                              \\end{cases}
    
    The annealing schedule is geometric. After each block of :math:`n`
    iterations at temperature :math:`T_i`, temperature :math:`T_{i+1}` is
    calculated as :math:`T_{i+1} = T_i \cdot (1 - \\gamma)` with
    :math:`\\gamma` the geometric temperature factor. The chain stops after
    reaching a limit temperature :math:`T_{\\text{min}}`; the maximal number
    of iterations is:
    
    .. math::
    
         n \\cdot \\left \\lceil
                             \\dfrac{\\log_{10}\\left(T_{\\text{min}}\\right)
                                   - \\log_{10}\\left(T_{\\text{ini}}\\right)}
                                    {\\log_{10}\\left(1-\\gamma\\right)}
                  \\right \\rceil
    
    The chain may stop sooner if no jump was made in :math:`3 \\cdot n` steps.
    
    To decrease the time spent sampling far away from the protein surface and
    avoid getting trapped into noise too early, the algorithm provides an
    optional adaptative step size mechanism. The step size increases
    when too many moves are accepted (*i.e.* the chain is sampling far away
    from a local minimum) and decreases when too many moves are rejected
    (*i.e.* the chain is sampling near or inside a local minimum). This
    approach was adapted from the variable step size method used to speed-up
    simulated annealing in continuous sample spaces.
    
    When the grid box is not cubic, the sphere will be stretched into an
    spheroid (oblate or prolate). This deformation does not preserve the
    uniform distribution, but the bias only becomes visible when one box
    dimension is twice as large as another dimension. At this point however,
    the protein topology is probably very complex, and a different
    distribution should probably be used.
    
    The simulated annealing algorithm was published in:
    
    * S. Kirkpatrick, C. D. Gelatt and M. P. Vecchi, Optimization by simulated
      annealing, *Science (New York, N.Y.)* **1983**, *220(4598)*: 671-680.
      DOI: `10.1126/science.220.4598.671
      <http://dx.doi.org/10.1126/science.220.4598.671>`_
    * V. Cerny, Thermodynamical approach to the traveling salesman problem: An
      efficient simulation algorithm, *Journal of Optimization Theory and Applications*
      **1985**, *45(1)*: 41--51. DOI: `10.1007/BF00940812
      <http://dx.doi.org/10.1007/BF00940812>`_
    
    The variable step size (VSS) method was published in:
    
    * Jon M. Sutter and John H. Kalivas, Convergence of generalized simulated
      annealing with variable step size with application towards parameter
      estimations of linear and nonlinear models, *Analytical Chemistry*
      **1991**, *63(20)*: 2383--2386. DOI: `10.1021/ac00020a034
      <http://dx.doi.org/10.1021/ac00020a034>`_
    
    Discussion on VSS:
    
    * Edwin E. Tucker, John H. Kalivas and Jon M. Sutter, Exchange of Comments
      on Convergence of generalized simulated annealing with variable step
      size with application towards parameter estimations of linear and
      nonlinear models, *Analytical Chemistry* **1992**, *64(10)*: 1199--1200.
      DOI: `10.1021/ac00034a021 <http://dx.doi.org/10.1021/ac00034a021>`_
    
    Application of VSS:
    
    * Kimberly Hitchcock, John H. Kalivas and Jon M. Sutter,
      Computer-generated multicomponent calibration designs for optimal
      analysis sample predictions, *Journal of Chemometrics* **1992**, *6(2)*:
      85--96. DOI: `10.1002/cem.1180060206
      <http://dx.doi.org/10.1002/cem.1180060206>`_
    
    :param box: energy grid
    :type  box: np.array
    :param runs: number of independent runs
    :type  runs: int
    :param Tini: initial temperature
    :type  Tini: float
    :param gamma: geometric temperature factor
    :type  gamma: float
    :param block_size: number of iterations between each temperature decrease
    :type  block_size: int
    :param convergence_T: stop when this temperature is reached
    :type  convergence_T: float
    :param convergence_rejections: stop when the chain did not move for that
       many steps, leave ``None`` to stop after 3 \* **block_size** or use 0
       to stop only once :math:`T_{\\text{min}}` is reached
    :type  convergence_rejections: int
    :param step_size: step size
    :type  step_size: int
    :param adaptative_step_size: if ``True``, the step size is changed
        dynamically, default is ``True``
    :type  adaptative_step_size: bool
    :param acceptance_ratio_upper: if **adaptative_step_size** is ``True``,
       minimum number of accepted moves in a temperature block before the step
       size is decremented
    :type  acceptance_ratio_upper: int
    :param acceptance_ratio_lower: if **adaptative_step_size** is ``True``,
       maximum number of accepted moves in a temperature block before the step
       size is incremented
    :type  acceptance_ratio_lower: int
    :param seed: numpy random generator seed (positive integer), value will be
       incremented after each run, so that run #1 with seed 500 is identical to
       run #11 with seed 490
    :type  seed: int
    :param return_full_chains: return the full Markov paths, default is to
       return only the converged solutions
    :type  return_full_chains: np.array
    
    :returns: List of local minima grid coordinates (n,l,m);
       if **return_full_chains** is ``True`` however, the list will take the
       form [(run, n, l, m, n_new, l_new, m_new, iterations in curent block,
       step size, accepted moves in current block), ...]
    :rtype: :class:`numpy.ndarray[:,3]` or :class:`numpy.ndarray[:,10]`
    '''
    if seed is None:
        seed = np.random.randint(500)
    if convergence_rejections is None:
        convergence_rejections = 3 * block_size
    return c_simulated_annealing(box, runs, Tini, gamma, block_size,
          convergence_T, convergence_rejections, step_size,
          adaptative_step_size, acceptance_ratio_upper,
          acceptance_ratio_lower, seed, return_full_chains)


cdef np.ndarray[np.int16_t] c_simulated_annealing(
    np.ndarray[np.float64_t, ndim=3] box, int runs, float Tini, float gamma,
    int block_size, float convergence_T, int convergence_rejections,
    int step_size, int adaptative_step_size, int acceptance_ratio_upper,
    int acceptance_ratio_lower, int seed, int return_full_chains):
    '''
    See :func:`simulated_annealing`.
    '''
    
    cdef int max_x = box.shape[0]
    cdef int max_y = box.shape[1]
    cdef int max_z = box.shape[2]
    cdef int x, y, z
    cdef int xnew, ynew, znew
    cdef int xshift, yshift, zshift
    cdef int delta
    cdef int delta_ini
    cdef int step
    cdef int accepted_moves_in_block
    cdef int consecutive_rejections
    cdef int acc
    cdef float T
    cdef float DeltaE
    cdef int acc_stop
    
    # maximal umber of iterations
    cdef int max_iter = int(block_size * np.ceil(np.log10(convergence_T)
                                               / np.log10(1 - gamma)))
    cdef int iterations
    
    # Markov chain stops after n consecutive rejections
    if convergence_rejections <= 0: # infinite
        acc_stop = max_iter
    else:
        acc_stop = convergence_rejections
    
    # if full Markov chains are requested, create a large 2D Numpy array
    # to store run, 
    cdef np.ndarray[np.int16_t, ndim=2] output
    cdef int output_incr = 0
    if return_full_chains:
        # values: run ID, x, y, z, xnew, ynew, znew,
        # iterations in curent block, delta, accepted moves in current block
        output = np.zeros([runs * max_iter, 10], dtype=np.int16)
    else:
        # values: xend, yend, zend
        output = np.zeros([runs, 3], dtype=np.int16)
    
    # GSA vs. VSGSA
    if adaptative_step_size:
        # initial step size is 1/30th the grid box dimensions, yet at least 3
        delta_ini = max(3, int(np.round((max_x + max_y + max_z) / 3. / 30.)))
    else:
        delta_ini = step_size
    
    # initial positions: uniform distribution on a sphere or spheroid
    # (depending on the box shape), stretched to less than half the
    # box dimensions
    sphere_lattice = np.array(spherical_lattice_default(runs), dtype=float)
    sphere_lattice[:,0] = np.fix(max_x / 2.5 * sphere_lattice[:,0]) + max_x//2
    sphere_lattice[:,1] = np.fix(max_y / 2.5 * sphere_lattice[:,1]) + max_y//2
    sphere_lattice[:,2] = np.fix(max_z / 2.5 * sphere_lattice[:,2]) + max_z//2
    cdef np.ndarray[np.int16_t, ndim=2] sphere = np.array(sphere_lattice,
                                                          dtype=np.int16)
    
    for run in range(runs):
        np.random.seed(seed + run)
        T = Tini
        x, y, z = sphere[run]
        delta = delta_ini
        step = 0
        accepted_moves_in_block = 0
        consecutive_rejections = 0
        iterations = 0
        # main loop, stop when the minimum temperature is reached or when the
        # system doesn't evolve anymore
        while T > convergence_T and consecutive_rejections < acc_stop:
            step += 1
            # decrease temperature if end of a temperature block
            if step == block_size:
                T *= (1.0 - gamma)
                if adaptative_step_size:
                    # too many jumps were accepted => we are probably far away
                    # from any minimum of energy and should increment the step
                    # size to reach the nearest minimum faster;
                    # too many jumps were rejected => we are probably in an
                    # energy minimum and should decrement the step size to
                    # stay near it
                    if accepted_moves_in_block > acceptance_ratio_upper:
                        delta += 1
                    elif accepted_moves_in_block < acceptance_ratio_lower:
                        delta = max(1, delta - 1)
                step = 0
                accepted_moves_in_block = 0
            
            # draw next position at random, excluding vector [0,0,0] (no jump)
            xshift, yshift, zshift = (0, 0, 0)
            while xshift == 0 and yshift == 0 and zshift == 0:
                xshift, yshift, zshift = np.random.randint(-delta, delta+1, 3)
            xnew = x + xshift
            ynew = y + yshift
            znew = z + zshift
            
            # store full Markov state if requested
            if return_full_chains:
                output[output_incr] = np.array([run + 1, x, y, z,
                                    xnew, ynew, znew, iterations//block_size,
                                    delta, accepted_moves_in_block])
                output_incr += 1
            
            # check if move is accepted
            consecutive_rejections += 1
            if xnew >= 0 and xnew < max_x and ynew >= 0 and ynew < max_y \
                                          and znew >= 0 and znew < max_z:
                Delta_E = box[xnew, ynew, znew] - box[x, y, z]
                if Delta_E <= 0 or np.math.exp(-Delta_E / T) > np.random.rand():
                    x = xnew
                    y = ynew
                    z = znew
                    accepted_moves_in_block += 1
                    consecutive_rejections = 0
            iterations += 1
        
        if not return_full_chains:
            output[run] = np.array([x,y,z])
    
    if return_full_chains:
        output = output[0:output_incr,:]
    return output

