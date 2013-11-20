'''
Created on 10.01.2012

@author: chris
'''
import os
import sys
import re
import gzip

import numpy as np


from epitopsy.cython import dx_cython

class DXBox:
    '''
    This class stores the output of an APBS calculation.
    '''
    X = 0
    Y = 1
    Z = 2

    # next neighbors are in direct contact!
    next_neighbors = np.zeros([6, 3], dtype = np.int)
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

    # neighbors are all around the origin/seed
    all_neighbors = np.zeros([26, 3], dtype = np.int)
    counter = 0
    for i in range(-1,2):
        for j in range(-1,2):
            for k in range(-1,2):
                if i == 0 and j == 0 and k == 0:
                    continue
                else:
                    all_neighbors[counter,0] = i
                    all_neighbors[counter,1] = j
                    all_neighbors[counter,2] = k
                    counter = counter + 1

    # all neighbors in xy!
    all_neighbors_xy = np.zeros([8, 3], dtype = np.int)
    counter = 0
    for i in range(-1,2):
        for j in range(-1,2):
            if i == 0 and j == 0:
                continue
            else:
                all_neighbors_xy[counter,0] = i
                all_neighbors_xy[counter,1] = j
                all_neighbors_xy[counter,2] = 0
                counter = counter + 1
    # [-1, -1, 0]
    all_neighbors_xy[0, 0] = -1
    all_neighbors_xy[0, 1] = -1
    # [-1, 0, 0]
    all_neighbors_xy[1, 0] = -1
    # [-1, 1, 0]
    all_neighbors_xy[2, 0] = -1
    all_neighbors_xy[2, 1] = 1
    # [0, -1, 0]
    all_neighbors_xy[3, 1] = -1
    # [0, 1, 0]
    all_neighbors_xy[4, 1] = 1
    # [1, -1, 0]
    all_neighbors_xy[5, 0] = 1
    all_neighbors_xy[5, 1] = -1
    # [1, 0, 0]
    all_neighbors_xy[6, 0] = 1
    # [1, 1, 0]
    all_neighbors_xy[7, 0] = 1
    all_neighbors_xy[7, 1] = 1

    def __init__(self, box = None, meshsize = None, offset = None):
        '''
        Creates a new DXBox Object.
        @param box Three-dimensional array with van-der-Waals surface or
        electrostatic potential of a protein structure
        '''
        if box is None and meshsize is None and offset is None:
            pass
        else:
            self.box = box
            self.box_dim = list(self.box.shape)
            self.box_mesh_size = meshsize
            self.box_offset = offset
            self.box_center = (np.ceil(np.array(self.box_dim) / 2.)
                               * np.array(self.box_mesh_size)
                               + np.array(self.box_offset))

        self.file_path = None

    def getDimensionX(self):
        if self.box_dim is not None:
            return self.box_dim[0]
        else:
            return None

    def getDimensionY(self):
        if self.box_dim is not None:
            return self.box_dim[1]
        else:
            return None

    def getDimensionZ(self):
        if self.box_dim is not None:
            return self.box_dim[2]
        else:
            return None

    def setDimensions(self, dim, y = 0, z = 0):
        if isinstance(dim, list):
            self.box_dim = dim
        else:
            self.box_dim = [dim, y, z]

    def  getDimensions(self):
        return self.box_dim

    def getMaxSideLength(self):
        if self.box_dim is not None:
            return max(self.box_dim)
        else:
            return None

    # I switched the variables compared to Niko,  because that way
    # I just need one function
    def makeCubic(self, fillInValue, sidelength = 'noValue'):
        if sidelength == 'noValue':
            sidelength = self.getMaxSideLength()
        self.enlargeBox([sidelength, sidelength, sidelength], fillInValue)

    def enlargeBox(self, dim, fillInValue):
        '''
        This method enlarges the box.
        '''
        raise AttributeError('Not implemented')

    def has_same_size_as(self, other_dxbox):
        if isinstance(other_dxbox, DXBox):
            dim = other_dxbox.getDimensions()
        else:
            raise AttributeError('Given object is not a DXBox!')

        return self.box_dim == other_dxbox.box_dim

    def translate(self, dx, dy = 0, dz = 0):
        '''
         * Translates all nodes by the given distance along each axis. The box is continuous,  thus nodes which are translated
         * over the right edge of the grid are reintroduced at the left edge.
         * <p>
         * Example:
         * <p>
         * With the grid x-Dimension being 97 and the node at x-positon 96 being translated by 3 along the x-axis,  the nodes
         * new x-position will be 2.
         * <code>x = (x+dx)%dimx</code>
         * TODO: changed,  check if new method is correct!
        '''
        if isinstance(dx, list):
            dz = dx[self.Z]
            dy = dx[self.Y]
            dx = dx[self.X]
        # shift in x-direction
        self.box = np.roll(self.box, dx, 0)
        # shift in y-direction
        self.box = np.roll(self.box, dy, 1)
        # shift in z-direction
        self.box = np.roll(self.box, dz, 2)

    def getNodesByValue(self, value):
        '''
        * Gets a list of all nodes which hold the given numerical value.
        '''
        print('Use np.nonzero(self.box == value)!')

    '''
     * Gets a list of all nodes which holds one of the given numerical values.
     * @param value array with values to searched for
     * @return Vector-List of all nodes with one of the supplied values
    '''

    def getNodesWithoutValue(self, value):
        '''
         * Gets a list of all nodes which DO NOT hold the given numerical value.
         * @param value value to NOT searche for
         * @return Vector-List of all nodes without supplied value
        '''
        print('Use np.nonzero(self.box != value)!')


    def multiply(self, factor_dxbox):
        if self.has_same_size_as(factor_dxbox):
            self.box = self.box * factor_dxbox.box
        else:
            print("Box sizes mismatch!!!")

    def subtract(self, substrahend_dxbox):
        if self.has_same_size_as(substrahend_dxbox):
            self.box = self.box - substrahend_dxbox.box
        else:
            print("Box sizes mismatch!!!")

    def determineEpitopeNodes(self, nodes, atoms, meshSize, radius):
        '''
        * Method takes a list of nodes and returns these which are
        * epitope-associated.
        * Nodes which reside within the supplied distance from at least
        * one of the supplied residues are considered epitope-associated.
        * <p>
        * Alogrithm:
        * compare every surface-associated node to every atom of every
        * epitope-associated residue:
        * If this node is within a given range of an atom,  it is
        * epitope-associated and added to the list.
        * @param nodes Vector-list of epitope candidate nodes
        * @param atoms Vector-list atoms which describe the epitope
        * @param meshSize Array with mesh spacing of this grid.
        * @param radius Radius of atoms within which nodes are considered epitope-associated
        * @return Vector-list of epitope-associated nodes
        *****
        Remark: I don't know what Niko is doing here exactly,  but I'll
        find that out ... eventually ;)
        *****
        '''
        # node list from e.g. self.getNodesWithoutValue(...)
        # -> simple 2d-list
        tolerableDistance = 15
        epitopeNodes = []
        atomVectors = []
        # cast as array,  dont't know the return type yet
        epitopeCenter = np.array()
        # TODO: check if this works
        for i in atoms.size():
            atomVectors.append(i.getCoords())
        centerNodeX = np.ceil(self.getDimensionX() / 2)
        centerNodeY = np.ceil(self.getDimensionY() / 2)
        centerNodeZ = np.ceil(self.getDimensionZ() / 2)
        # TODO: don't know what jama of java does here
        gridCenterNode = np.array([centerNodeX, centerNodeY, centerNodeZ])
        meshDim = np.array(meshSize)
        nodeID = 0
        # TODO: len or size,  don't know yet
        while nodeID < len(nodes):
            nodePositionVector = np.dot((nodes[i] - gridCenterNode), meshDim)
            nodeHighlighted = False
            tooFarAway = False
            atomID = 0
            while atomID < len(atomVectors) and not(nodeHighlighted) and not(tooFarAway):
                diffVector = epitopeCenter - nodePositionVector
                if np.sqrt(np.dot(diffVector, diffVector)) <= tolerableDistance:
                    currentAtomCoord = np.array(atomVectors[atomID])
                    diffVector2 = currentAtomCoord - nodePositionVector
                    if np.sqrt(np.dot(diffVector2, diffVector2)):
                        epitopeNodes.append(nodes[nodeID])
                        nodeHighlighted = True
                    atomID = atomID = 1
                else:
                    tooFarAway = True
            nodeID = nodeID + 1
        return epitopeNodes

    def count_nodes_with_value(self, value):
        '''
        Returns number occurrences of a given value within this grid.
        '''
        box_coord = np.nonzero(self.box == value)

        return len(box_coord[0])

    def count_nodes_without_value(self, value):
        box_coord = np.nonzero(self.box != value)

        return len(box_coord[0])


    def getUniqueValues(self):
        '''
        Returns unique values within this grid.
        *****
        Remark: Nikos code is very complex here (with thresshold and so on),
        so I don't know if I miss something,  though the name of the function
        is very clear
        '''
        unique = []
        for x in range(self.getDimensionX()):
            for y in range(self.getDimensionY()):
                for z in range(self.getDimensionZ()):
                    if not(self.box[x][y][z] in unique):
                        unique.append(self.box[x][y][z])
        return unique

    def getRotateBox(self, theta, phi, psi, fillInValue):
        '''
         * Return a rotated copy of this grid with nodes rotated around the grid center.
         * <p>
         * This can be seen as rotating the protein and thus the specific information (van-der-Waals surface
         * or electrostatic potential) within this grid.

         I use euler-angles, so it is actually phi, theta, psi (x, y, z).
         You can look up euler angles at wikipedia to understand what
         I do here, it is different from Niko's!
        '''
        # preparation start:
        centerNodeVector = np.array([np.ceil(self.getDimensionX() / 2), np.ceil(self.getDimensionY() / 2), np.ceil(self.getDimensionZ() / 2)])
        '''
        http://mathworld.wolfram.com/EulerAngles.html

        phi = 0 - 360 about the z-axis
        theta = 0 - 180 about the new x'-axis
        psi = 0 - 360 about the new z'-axis
        '''
        eulerRotation = np.zeros((3, 3))
        theta = np.deg2rad(theta)
        phi = np.deg2rad(phi)
        psi = np.deg2rad(psi)
        eulerRotation[0][0] = np.cos(psi) * np.cos(phi) - np.cos(theta) * np.sin(phi) * np.sin(psi)
        eulerRotation[0][1] = np.cos(psi) * np.sin(phi) + np.cos(theta) * np.cos(phi) * np.sin(psi)
        eulerRotation[0][2] = np.sin(psi) * np.sin(theta)
        eulerRotation[1][0] = -np.sin(psi) * np.cos(phi) - np.cos(theta) * np.sin(phi) * np.cos(psi)
        eulerRotation[1][1] = -np.sin(psi) * np.sin(phi) + np.cos(theta) * np.cos(phi) * np.cos(psi)
        eulerRotation[1][2] = np.cos(psi) * np.sin(theta)
        eulerRotation[2][0] = np.sin(theta) * np.sin(phi)
        eulerRotation[2][1] = -np.sin(theta) * np.cos(phi)
        eulerRotation[2][2] = np.cos(theta)
        # preparation end!

        rb = np.zeros((self.getDimensionX(), self.getDimensionY(), self.getDimensionZ()))
        for x in range(self.dimensionX):
            for y in range(self.dimensionY):
                for z in range(self.dimensionZ):
                    rb[x][y][z] = fillInValue
        for x in range(self.dimensionX):
            for y in range(self.dimensionY):
                for z in range(self.dimensionZ):
                    if self.box[x][y][z] != fillInValue:
                        v = np.array([x, y, z])
                        gridPos = np.array([0.0, 0.0, 0.0])
                        # translate:
                        v = v - centerNodeVector
                        v = np.dot(eulerRotation, v)
                        v = v + centerNodeVector
                        gridPos[self.X] = round(v[self.X])
                        gridPos[self.Y] = round(v[self.Y])
                        gridPos[self.Z] = round(v[self.Z])
                        try:
                            self.setRotatedNode(rb, gridPos[0], gridPos[1],
                                                gridPos[2], self.box[x][y][z])
                        except:
                            print("ROTATING NODE: {0}/{1}/{2}\n".format(x, y, z))
                            print("rotated node with value {0}\n".format(self.box[x][y][z]))
                            print("to   position {0} {1} {2}\n".format(gridPos[0],
                                                                       gridPos[1],
                                                                       gridPos[2]))
        return rb


    def determineConsensoursDimensions(self, grid1, grid2):
        '''
        Returns the consensus grid size of two grids.
        '''
        cons_size = []
        cons_size.append(max(grid1.box_dim[0], grid2.box_dim[0]))
        cons_size.append(max(grid1.box_dim[1], grid2.box_dim[1]))
        cons_size.append(max(grid1.box_dim[2], grid2.box_dim[2]))
        return cons_size

    def add_box(self, box):
        if not(self.has_same_size_as(box)):
            print("ERROR: Box dimensions differ!")
        else:
            self.box = self.box + box

    def divide_by_scalar(self, d):
        self.box = self.box / d

    def transform_real_to_box_space(self, atom_coord):
        '''
        This method transforms coordinates from the real to the box space.
        '''
        return np.around((np.array(atom_coord) - np.array(self.box_offset)) / np.array(self.box_mesh_size)).astype(int)

    def transform_box_to_real_space(self, grid_coord):
        '''
        This method transforms coordinates from the box to the real space.
        '''
        return np.array(grid_coord) * np.array(self.box_mesh_size) + np.array(self.box_offset)

    def indices_in_box(self, indices):
        '''
        This method checks, if the given indices lie inside the box, i.e
        they lie in: 0 <= index < dimensionXYZ.
        '''
        x = indices[0]
        y = indices[1]
        z = indices[2]
        if (not 0 <= x < self.box_dim[0] or not 0 <= y < self.box_dim[1] or
            not 0 <= z < self.box_dim[2]):
            return False
        else:
            return True


    def write_deprecated(self, filename, values_per_line = 3):
        '''
        Write an array to a dx-file.
        '''
        if filename.endswith('.gz'):
            with gzip.open(filename, 'wb') as outfile:
                self._write_deprecated(outfile, values_per_line)

        else:
            with open(filename, 'w') as outfile:
                self._write_deprecated(outfile, values_per_line)

    def _write_deprecated(self,outfile, values_per_line):
        valueCount = 0
        header = self._generate_header()
        footer = self._generate_footer()
        outfile.write(header)
        for x in range(self.box_dim[0]):
            for y in range(self.box_dim[1]):
                for z in range(self.box_dim[2]):
                    outfile.write("{0:e} ".format(self.box[x, y, z]))
                    valueCount = valueCount + 1
                    if valueCount == values_per_line:
                        outfile.write('\n')
                        valueCount = 0
        outfile.write('\n')
        outfile.write(footer)

    def write(self, filename, values_per_line = 3):
        '''
        Write an array to a dx-file.
        '''
        if filename.endswith('.gz'):
            with gzip.open(filename, 'wb') as outfile:
                self._write(outfile, values_per_line)

        else:
            with open(filename, 'w') as outfile:
                self._write(outfile, values_per_line)

    def _write(self,outfile, values_per_line):
        header = self._generate_header()
        footer = self._generate_footer()
        outfile.write(header)
        dx_cython.write_box(outfile, self.box, values_per_line)
        outfile.write('\n')
        outfile.write(footer)

    def _generate_header(self):
        header = "#DX generated by Epitopsy\n"
        header = header + "object 1 class gridpositions counts {0} {1} {2}\n".format(self.box_dim[0], self.box_dim[1], self.box_dim[2])
        header = header + "origin {0:e} {1:e} {2:e}\n".format(self.box_offset[0], self.box_offset[1], self.box_offset[2])
        header = header + "delta {0:e} 0.000000e+00 0.000000e+00\n".format(self.box_mesh_size[0])
        header = header + "delta 0.000000e+00 {0:e} 0.000000e+00\n".format(self.box_mesh_size[1])
        header = header + "delta 0.000000e+00 0.000000e+00 {0:e}\n".format(self.box_mesh_size[2])
        header = header + "object 2 class gridconnections counts {0} {1} {2}\n".format(self.box_dim[0], self.box_dim[1], self.box_dim[2])
        header = header + "object 3 class array type double rank 0 items {0} data follows\n".format(self.box_dim[0] * self.box_dim[1] * self.box_dim[2])
        return header

    def _generate_footer(self):
        '''
        This method generates the information at the end of a dx file.
        '''
        footer = "attribute \"dep\" string \"positions\"\n"
        footer = footer + "object \"regular positions regular connections\" class field\n"
        footer = footer + "component \"positions\" value 1\n"
        footer = footer + "component \"connections\" value 2\n"
        footer = footer + "component \"data\" value 3\n"
        return footer

def read_dxfile(filename, box_type):
    [box_3d, meshsize, offset] = dx_cython.read_dxfile(filename, box_type)
    if box_type == "vdw":
        dxBox = VDWBox(box_3d, meshsize, offset)
        dxBox.file_path = filename
    elif box_type == "smol":
        dxBox = VDWBox(box_3d, meshsize, offset)
        dxBox.file_path = filename
    elif box_type == "esp":
        dxBox = DXBox(box_3d, meshsize, offset)
        dxBox.file_path = filename
    return dxBox

class DXReader:
    '''
    classdocs
    '''
    def __init__(self):
        '''
        Constructor
        '''
        self.dxFilename = None
        self.VDW = 'vdw'
        self.ESP = 'esp'
        self.SMOL = 'smol'
        self.X = 0
        self.Y = 1
        self.Z = 2

    def parse_old(self, filename, boxType):
        '''
        Method to parse supplied DX file. This one is rather slow!

        Returns a VDWBox for 'vdw' and 'smol' box types, otherwise it returns
        a DXBox.
        '''
        countX = 0
        countY = 0
        countZ = 0
        dimension = [0, 0, 0]
        readingHeader = True
        header = ""
        footer = ""
        offset = [0.0, 0.0, 0.0]
        meshsize = [0.0, 0.0, 0.0]
        # this was difficult!
        dimensionPattern = re.compile('(\\d+) (\\d+) (\\d+)$')
        lastHeaderLinePattern = re.compile('^object 3')
        objectPattern = re.compile('^object 1')
        valuePattern = re.compile('^[-]?\\d+\\.\\d+[eE][+-]\\d+')
        # this has been tested and it works
        originPattern = re.compile('^origin ([-]?\\d+.\\d+[eE]?[+-]?\\d+?) ([-]?\\d+.\\d+[eE]?[-+]?\\d+?) ([-]?\\d+.\\d+[eE]?[+-]?\\d+?)')
        deltaPattern = re.compile('^delta ([-]?\\d+.\\d+[eE]?[+-]?\\d+?) ([-]?\\d+.\\d+[eE]?[-+]?\\d+?) ([-]?\\d+.\\d+[eE]?[+-]?\\d+?)')
        with open(filename, 'r') as infile:
            for line in infile:
                if readingHeader == True:
                    if objectPattern.match(line):
                        matchedObject = dimensionPattern.search(line)
                        dimension[0] = int(matchedObject.group(1))
                        dimension[1] = int(matchedObject.group(2))
                        dimension[2] = int(matchedObject.group(3))
                        if boxType == self.VDW:
                            # box = np.zeros(dimension, dtype = int)
                            box = np.zeros(dimension)
                        elif boxType == self.ESP:
                            box = np.zeros(dimension)
                        elif boxType == self.SMOL:
                            box = np.zeros(dimension)
                    if originPattern.match(line):
                        '''
                        Origin = Offset!
                        '''
                        # don't use regular expressions here, they are a mess!
                        matchedOrigin = line.split(' ')
                        offset[0] = float(matchedOrigin[1])
                        offset[1] = float(matchedOrigin[2])
                        offset[2] = float(matchedOrigin[3])
                    if line[0:5] == 'delta':
                        '''
                        Read mesh size.
                        delta = meshsize!
                        '''
                        for i in range(3):
                            delta_split = line.split()
                            meshsize[i] = float(delta_split[i + 1])
                            line = infile.next()
                    header = header + line + '\n'
                    if lastHeaderLinePattern.match(line):
                        readingHeader = False
                else:
                    if valuePattern.match(line):
                        nodes = line.split(' ')
                        # -1, because the last element is '\n'
                        for i in range(len(nodes) - 1):
                            '''
                            reads the elements in the following way:
                                1. z-elements first [0][0][increment this one]
                                2. if first z is done, read y
                                    [0][increment this one][]
                                3. x is the last one
                                    [increment this one][][]

                            '''
                            if countZ == dimension[self.Z]:
                                countZ = 0
                                if countY == (dimension[self.Y] - 1):
                                    countY = 0
                                    countX = countX + 1
                                else:
                                    countY = countY + 1
                            if boxType == self.VDW:
                                box[countX][countY][countZ] = float(nodes[i])
                            elif boxType == self.ESP:
                                box[countX][countY][countZ] = float(nodes[i])
                            countZ = countZ + 1
                    else:
                        footer = footer + line + '\n'

        if boxType == self.VDW or boxType == self.SMOL:
            dxBox = VDWBox(box, meshsize, offset)
        else:
            dxBox = DXBox(box, meshsize, offset)

        return dxBox

    def parse(self, filename, boxType):
        return read_dxfile(filename, boxType)

    def _parse_depreciated(self, filename, boxType):
        '''
        Method to parse supplied DX file. This one is twice as fast as the one
        above (but this depends on the size of the DXFile of course)!
        '''
        if not(os.path.isabs(filename)):
            filename = os.path.join(os.getcwd(), filename)

        countX = 0
        countY = 0
        countZ = 0
        counter = 0
        dimension = [0, 0, 0]
        offset = [0.0, 0.0, 0.0]
        meshsize = [0.0, 0.0, 0.0]
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

        delta_counter = 0
        for line in content:
            if not line.startswith('#'):
                if line.startswith('object 1'):
                    matchedObject = dimensionPattern.search(line)
                    dimension[0] = int(matchedObject.group(1))
                    dimension[1] = int(matchedObject.group(2))
                    dimension[2] = int(matchedObject.group(3))
                    box = np.zeros(dimension[0] * dimension[1] * dimension[2])
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

        box = box.reshape(dimension)

        if boxType == self.VDW:
            dxBox = VDWBox(box, meshsize, offset)
            dxBox.file_path = filename
        elif boxType == self.SMOL:
            dxBox = VDWBox(box, meshsize, offset)
            dxBox.file_path = filename
        elif boxType == self.ESP:
            dxBox = DXBox(box, meshsize, offset)
            dxBox.file_path = filename

        return dxBox


class VDWBox(DXBox):
    """
    classdocs
    """
    SOLVENT_ID_VALUE = 99.0
    INNER_ID_VALUE = 98
    SURFACE_ID_VALUE = 97

    def __init__(self, box, box_mesh_size, box_offset):
        """
        Constructor
        """
        self.peptideScore = 0.0
        self.surfaceScore = 0.0
        self.solventScore = 1.0
        self.epitopeScore = 0.0
        self.protein_interior = 0.0
        # this score is assigned after 'find_surface'
        self.score_of_surface = 2.0
        DXBox.__init__(self, box, box_mesh_size, box_offset)

    def assignNodeIDs(self):
        self.applySurfaceScore(self.SURFACE_ID_VALUE)
        self.applyPeptideScore(self.INNER_ID_VALUE)
        self.applySolventScore(self.SOLVENT_ID_VALUE)

    def applyPeptideScore(self, score):
        """
         * Applies a given score to peptide.
         * Method will retrieve actual peptide(peptide nodes without solvent access) if surface- and peptide score are identical.
         * If surface-score is different from peptide score(peptide nodes are thus distinguishable from others) given score
         * is applied to those.
         * @param score
        """
        if self.peptideScore != score:
            peptideNodes = []
            if self.surfaceScore == self.peptideScore:
                peptideNodes = self.determineActualInnerPeptideNodes()
            else:
                peptideNodes = self.getPeptideNodes()
            self.setElementsRange(peptideNodes, score)
            self.peptideScore = score

    def applySolventScore(self, score):
        if self.solventScore != score:
            for x in range(self.dimensionX):
                for y in range(self.dimensionY):
                    for z in range(self.dimensionZ):
                        if self.getElement(x, y, z) == self.solventScore:
                            self.setElement(score, x, y, z)

    def applySurfaceScore(self, score):
        """
         * Applies a given score to surface.
         * Method will retrieve actual surface(peptide nodes with solvent access) if surface- and peptide score are identical.
         * If surface-score is different from peptide score(surface nodes are thus distinguishable from others) given score
         * is applied to those.
         * @param score
        """
        if self.surfaceScore != score:
            surfaceNodes = []
            if self.surfaceScore == self.peptideScore:
                surfaceNodes = self.determineActualSurfaceNodes()
            else:
                surfaceNodes = self.getSurfaceNodes()
            self.setElementsRange(surfaceNodes, score)
            self.surfaceScore = score
            if self.epitopeScore == self.surfaceScore:
                self.epitopeScore = score

    def clone(self):
        cp = VDWBox(self.box.copy(), self.box_mesh_size, self.box_offset)
        cp.peptideScore = self.peptideScore
        cp.solventScore = self.solventScore
        return cp

    def createPotentialArea(self, areaWidth):
        """
         * Artificially enlarges the protein surface into surrounding solvent. Then sets score within original and
         * artificial surface to 1. Everything else is set to 0.
         * @param areaWidth
        """
        self.assignNodeIDs()
        self.applySurfaceScore(1.0)
        self.smoothSurface(areaWidth, 1.0)
        self.applySolventScore(0.0)
        self.applyPeptideScore(0.0)

    def enlargeBox(self, dim):
        # TODO: Does this work? probably not ;)
        DXBox.enlargeBox(self, dim, self.solventScore)

    def determineActualInnerPeptideNodes(self):
        """
         * This Method determines and return an array containing all nodes without solvent access.
         * Algorithm works as follows:
         * If a node is non solvent-accessible(surface or inner node),
         * and has no neighbors which are solvent,
         * the node is an inner node.
         * @return
        """
        return dx_cython.determineActualInnerPeptideNodes(self.box,
            self.solventScore, self.peptideScore, self.all_neighbors)

    def determineActualSurfaceNodes(self):
        """
         * This Method determines and returns an array containing all nodes with direct solvent access.
         * @return
        """
        return dx_cython.determine_actual_surface_nodes(self.box, self.solventScore, self.all_neighbors)

    def determineEpitopeSurfaceNodes(self, atoms, meshSize, rangeCutoff):
        """
         * TODO: made use of a more generic method to determine epitope nodes
         * Method resides in super-class DXBox ->determine Epitope nodes.
        """
        self.determineEpitopeNodes(self, self.getSurfaceNodes(), atoms, meshSize, rangeCutoff)

    def flood_old_slow(self):
        temp_value = -1.0
        #for z in range(self.box.shape[2]):
        #    self.floodfill(0, 0, 0, temp_value)
        waterFront = [np.array([0, 0, 0])]
        newFront = []
        tempValue = -1.0
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
        for i in range(self.box.shape[2]):
            waterFront = [np.array([0, 0, i])]
            while(len(waterFront) > 0):
                if len(newFront) > 0:
                    newFront = []
                for node in waterFront:
                    neighbors = self.all_neighbors_xy
                    for neighbor in neighbors:
                        n = node + neighbor
                        try:
                            if (n[0] >= 0 and n[0] < self.dimensionX
                                and n[1] >= 0 and n[1] < self.dimensionY):
                                if self.box[n[0]][n[1]][n[2]] == self.solventScore:
                                    self.box[n[0]][n[1]][n[2]] = tempValue
                                    newFront.append(n)
                        except Exception:
                            print("Error out of bounds!")
                if len(waterFront) > 0:
                    waterFront = []
                waterFront = newFront[:]
        for x in range(self.box.shape[0]):
            for y in range(self.box.shape[1]):
                for z in range(self.box.shape[2]):
                    if self.box[x][y][z] == self.solventScore:
                        self.box[x][y][z] = self.peptideScore
                    elif self.box[x][y][z] == temp_value:
                        self.box[x][y][z] = self.solventScore

    def flood(self):
        """
        This method is the same as flood_old_slow, but it has been written in
        cython and is therefore at least twice as fast!
        """
        self.box = dx_cython.floodxy(self.box, self.all_neighbors_xy)

    def calculate_sas(self):
        '''Calculate the solvent accessible surface.
        '''
        self.box = dx_cython.flood(self.box, self.next_neighbors)
        self.box = dx_cython.get_sas(self.box, self.peptideScore, self.solventScore, 1.4)
        self.box = dx_cython.flood(self.box, self.next_neighbors)

    def find_surface(self):
        """
        This method finds the surface of the protein and assigns
        'score_of_surface' to these points. The box should have been flooded
        before using this function!

        In contrast, find_solvent surface yields the surface that lies in the
        surface. Both surfaces should lie directly next to each other.
        """
        self.box = dx_cython.find_protein_surface(self.peptideScore,
                                          self.solventScore,
                                          self.score_of_surface,
                                          self.box)

    def find_solvent_surface(self):
        """
        This method finds the surface that lies in the solvent and
        encapsulates the protein and assigns 'score_of_surface' to these
        points. The box should have been flooded before using this function!

        In contrast, find_surface yields the surface that lies in the protein.
        Both surfaces should lie directly next to each other.
        """
        self.box = dx_cython.find_solvent_surface(self.peptideScore,
                                          self.solventScore,
                                          self.score_of_surface,
                                          self.box)

    def prepare_for_geometric_matching(self, interior):
        """
        This method prepares the box for geometric matching. It sets the
        solvent to 0, the surface to 1 and the protein interior to the given
        value. This value should be -15 for the fixed protein and 1 for the
        rotated according to Katchalski-Katzir.

        The values for 'solventScore', 'peptideScore', etc. are not updated!
        """
        # set surface to self.score_of_surface
        self.find_surface()
        surface_points = np.nonzero(self.box == self.score_of_surface)
        solvent_points = np.nonzero(self.box == self.solventScore)
        interior_points = np.nonzero(self.box == self.peptideScore)
        # set protein interior to interior
        self.box[interior_points] = interior
        # set solvent to 0
        self.box[solvent_points] = 0
        # set surface to 1
        self.box[surface_points] = 1
        # done!

    def getPeptideNodes(self):
        """
         * This method return a vector containing all nodes with peptide score.
         * @return
        """
        self.getNodesByValue(self.peptideScore)

    def getSolventNodes(self):
        self.getNodesByValue(self.solventScore)

    def getSurfaceNodes(self):
        """
         * This method returns a Vector containing all nodes with surface score.
         * @return
        """
        self.getNodesByValue(self.surfaceScore)

    def extendSurface(self, iterations):
        """
         * Method to artificially enlarge the protein surface into surrounding solvent.
         * @param iterations
        """
        self.box = dx_cython.extendSurface(int(iterations), self.box,
            self.solventScore, self.surfaceScore, self.next_neighbors)

    def extend_nice_surface(self, iterations, dummy_score = 1000):
        if(dummy_score == self.solventScore
           or dummy_score == self.surfaceScore):
            raise ValueError('Somehow the dummy_score in extend_nice_surface was the same as the solvent or the protein score!')

        self.box = dx_cython.extend_nice_surface(int(iterations), self.box,
            self.solventScore, self.surfaceScore, dummy_score,
            self.next_neighbors)

    def getNeighborList(self, seed, radius = 1):
        """
         * Generates a List of nodes surrounding a seed-node (0, 0, 0) within a given radius.
         * @param radius
         * @return
        """
        if not(isinstance(seed, list)):
            # received just the radius,  therefore seed is set to [0, 0, 0]
            radius = seed
            seed = [0., 0., 0.]
        neighbors = []
        # cast radius as int
        radius = int(round(radius))
        # +1, because otherwise you'll lose a number!
        # e.g. range(-1,1) = [-1, 0]
        radiusRange = range(-radius, radius + 1)
        for x in radiusRange:
            for y in radiusRange:
                for z in radiusRange:
                    # do not add the seed itself
                    if [x, y, z] != [0, 0, 0]:
                        neighbor_coord = np.array([x, y, z]) + np.array(seed)
                        neighbors.append(list(neighbor_coord))
        return neighbors[:]

    def highlightEpitope(self, atoms, meshSize, rangeCutoff, score):
        epitopeNodes = self.determineEpitopeSurfaceNodes(atoms, meshSize, rangeCutoff)
        print("{0} nodes markes as epitope associated".format(len(epitopeNodes)))
        self.setElementsRange(epitopeNodes, score)
        self.epitopeScore = score

    def invertScore(self):
        """
         * Method to invert scores of protein/peptide and solvent.
        """
        for x in range(self.dimensionX):
            for y in range(self.dimensionY):
                for z in range(self.dimensionZ):
                    if self.getElement(x, y, z) == self.peptideScore:
                        self.setElement(self.solventScore, x, y, z)
                    elif self.getElement(x, y, z) == self.solventScore:
                        self.setElement(self.peptideScore, x, y, z)
        dummy = self.peptideScore
        self.peptideScore = self.solventScore
        self.solventScore = dummy
    def getRotatedBox(self, theta, phi, psi):
        rb = VDWBox(DXBox().getRotatedBox (self, theta, phi, psi, self.solventScore), self.getDimensions(), self.getMeshSize(), self.getOffset())
        rb.solventScore = self.solventScore
        rb.peptideScore = self.peptideScore
        rb.surfaceScore = self.surfaceScore
        rb.epitopeScore = self.epitopeScore
        rb.setDimensions(self.getDimensions())
        return rb
    def setRotatedNode(self, rb, x, y, z, value):
        rb[x][y][z] = rb[x][y][z] + value
    """
     * Method to alter score for a larger set of nodes.
     * @param value
    """
    def setElementsRange(self, nodes, value):
        for node in nodes:
            self.setElement(value, node)
    def smoothSurface(self, iterations, gradient = 1.0):
        if iterations > 0:
            if self.surfaceScore == self.peptideScore:
                currentNodeSet = self.determineActualSurfaceNodes()
            else:
                currentNodeSet = self.getSurfaceNodes()
            if self.epitopeScore != self.surfaceScore:
                currentNodeSet = list(np.array(currentNodeSet) + np.array(self.getNodesByValue(self.epitopeScore)))
            # smooth surface for n times OR until smoothing value is gets lower than molecules interior value
            for i in range(iterations):
                newNodeSet = []
                for node in currentNodeSet:
                    neighbors = self.getNeighborList(node, 1)
                    for neighbor in neighbors:
                        element = self.getElement(neighbor)
                        if element == self.solventScore:
                            nodeValue = self.getElement(node)
                            self.setElement(round(nodeValue / gradient, 2), neighbor)
                            newNodeSet.append(neighbor)
                currentNodeSet = newNodeSet