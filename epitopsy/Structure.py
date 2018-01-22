__author__     = "Christoph Wilms, Jean-Noel Grad"
__copyright__  = "Copyright 2012, Epitopsy and Biopython"
__date__       = "2012-01-11"
__credits__    = ["Christoph Wilms", "Jean-Noel Grad"]
__license__    = '''
                 Biopython License Agreement

Permission to use, copy, modify, and distribute this software and its
documentation with or without modifications and for any purpose and
without fee is hereby granted, provided that any copyright notices
appear in all copies and that both those copyright notices and this
permission notice appear in supporting documentation, and that the
names of the contributors or copyright holders not be used in
advertising or publicity pertaining to distribution of the software
without specific prior permission.

THE CONTRIBUTORS AND COPYRIGHT HOLDERS OF THIS SOFTWARE DISCLAIM ALL
WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL THE
CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY SPECIAL, INDIRECT
OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
OR PERFORMANCE OF THIS SOFTWARE.
'''

import os
import gzip
import time
import numpy as np
import subprocess
from Bio.PDB.Atom import Atom as BioAtom
from Bio.PDB.Residue import Residue
from Bio.PDB.Chain import Chain
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
#from Bio.PDB import PDBParser, PDBIO, Select
from epitopsy.DXFile import DXBox, VDWBox
from epitopsy.cython import pdb_cython
from epitopsy.tools import MathTools

three2oneletter = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'ASN': 'N',
                       'GLN': 'Q',
                       'LYS': 'K', 'MLY':'K', 'M2L':'K', 'M3L':'K', 'ACK':'K', 'KCX':'K',
                       'ILE': 'I', 'PRO': 'P',
                       'THR': 'T', 'PHE': 'F', 'ALA': 'A', 'GLY': 'G',
                       'HIS': 'H', 'HIP':'H', 'HSP':'H',
                       'LEU': 'L',
                       'ARG': 'R', 'ARM':'R', 'SRM':'R', 'DRM':'R',
                       'TRP': 'W',
                       'VAL': 'V', 'GLU': 'E', 'TYR': 'Y',
                       'MET': 'M', 'MSE':"M"}

class Structure_Template(object):
    '''
    This is a template class from which :class:`PDBFile` and :class:`PQRFile`
    will be derived.

    This class is not supposed to work with multiple models (for example NMR
    structures)!
    '''
    # dictionary to convert three letter codes to one letter codes and vice
    # versa
    three2oneletter = three2oneletter

    one2threeletter = {'C':'CYS', 'D':'ASP', 'S':'SER', 'N':'ASN',
                       'Q':'GLN', 'K':'LYS', 'I':'ILE', 'P':'PRO',
                       'T':'THR', 'F':'PHE', 'A':'ALA', 'G':'GLY',
                       'H':'HIS', 'L':'LEU', 'R':'ARG', 'W':'TRP',
                       'V':'VAL', 'E':'GLU', 'Y':'TYR', 'M':'MET'}

    all_amino_acids = ['C', 'D', 'S', 'N', 'Q', 'K', 'I', 'P', 'T', 'F', 'A',
                       'G', 'H', 'L', 'R', 'W', 'V', 'E', 'Y', 'M']

    kyte_doolittle = {'A': 1.8, 'C': 2.5, 'D': -3.5, 'E': -3.5, 'F': 2.8,
                      'G': -0.4, 'H': -3.2, 'I': 4.5, 'K': -3.9, 'L': 3.8,
                      'M': 1.9, 'N': -3.5, 'P': -1.6, 'Q': -3.5, 'R': -4.5,
                      'S': -0.8, 'T': -0.7, 'V': 4.2, 'W': -0.9, 'Y': -1.3}

    def __init__(self, structure_path=None):
        self.structure_path = structure_path

        self.structure = None
        # 'pdb', 'pqr', ?
        self.what_am_i = None
        self.remarks = []

    def get_all_atom_coords(self):
        '''
        Navigate through all ATOM and HETATM records in the structure and
        return their coordinates. See :meth:`get_amino_atoms_coords` and
        :meth:`get_hetero_atoms_coords` to return all ATOM resp. HETATM
        coordinates.

        :returns: An iterable of the requested atom coordinates.
        :rtype: list(:class:`np.array[3]`)

        Example::

            >>> a = Structure.PDBFile('4N6W.pdb')
            >>> for ATOM in a.get_all_atom_coords():
            ...     print ATOM
            [-42.316  38.638  20.425]
            [-41.82   39.723  19.581]
            [-41.93   41.072  20.291]
            [-40.962  41.836  20.357]
            [-42.577  39.762  18.255]
            [-42.085  40.845  17.3  ]
            [-40.299  40.776  17.008]
            [-40.081  42.082  15.79 ]
            [-43.122  41.359  20.805]
            [-43.34   42.492  21.695]
            [...]

        '''
        atoms_coord = []
        for atom in self.structure.get_atoms():
            atoms_coord.append(atom.get_coord())
        return atoms_coord

    def get_coords_from_atom_list(self, atom_list):
        '''
        Navigate through all atoms in **atom_list** and return their
        coordinates.

        :param atom_list: atoms to read
        :type  atom_list: list(:class:`Atom`)

        :returns: An iterable of the requested atom coordinates.
        :rtype: list(:class:`np.array[3]`)

        Example::

            >>> a = Structure.PDBFile('4N6W.pdb')
            >>> b = a.get_hetero_atoms()[1:3] # retrieving two atom objects
            >>> print b
            [<Atom FE>, <Atom CAC>]
            >>> print a.get_coords_from_atom_list(b) # getting coordinates
            [array([-28.139,  20.588,  11.849]), array([-30.158,  21.541,   5.145])]
            >>> print a.get_hetero_atoms_coords()[1:3] # just checking
            [array([-28.139,  20.588,  11.849]), array([-30.158,  21.541,   5.145])]

        '''
        atom_coord = []
        for atom in atom_list:
            atom_coord.append(atom.get_coord())
        return atom_coord

    def get_hetero_atoms(self):
        '''
        Navigate through all HETATM records in the structure.

        :returns: All HETATM in the structure.
        :rtype: list(:class:`Atom`)

        Example::

            >>> a = Structure.PDBFile('4N6W.pdb')
            >>> for HETATM in a.get_hetero_atoms():
            ...     print HETATM
            <Atom FE>
            <Atom FE>
            <Atom CAC>
            <Atom CA>
            <Atom CB>
            <Atom CBC>
            <Atom CG>
            <Atom CGC>
            <Atom OA1>
            <Atom OA2>
            [...]

        '''
        het_atoms = []
        for chain in self.structure:
            for res in chain:
                if res.id[0].upper() == 'H':
                    for atom in res:
                        het_atoms.append(atom)

        return het_atoms

    def get_hetero_atoms_coords(self):
        '''
        Navigate through all HETATM records in the structure and return
        their coordinates.

        :returns: An iterable of the requested atom coordinates.
        :rtype: list(:class:`np.array[3]`)

        Example::

            >>> a = Structure.PDBFile('4N6W.pdb')
            >>> for HETATM in a.get_hetero_atoms_coords():
            ...     print HETATM
            [-28.886  22.883   9.035]
            [-28.139  20.588  11.849]
            [-30.158  21.541   5.145]
            [-30.083  20.273   5.97 ]
            [-29.335  20.464   7.282]
            [-27.967  21.06    7.026]
            [-29.217  19.094   7.918]
            [-28.359  19.09    9.162]
            [-29.115  22.097   4.693]
            [-31.292  22.037   4.902]
            [...]

        '''
        het_coord = []
        for chain in self.structure:
            for res in chain:
                if res.id[0].upper() == 'H':
                    for atom in res:
                        het_coord.append(atom.get_coord())

        return het_coord

    def get_amino_atoms(self):
        '''
        Navigate through all ATOM records in the structure.

        :returns: All ATOM in the structure.
        :rtype: list(:class:`Atom`)

        Example::

            >>> a = Structure.PDBFile('4N6W.pdb')
            >>> for ATOM in a.get_amino_atoms():
            ...     print ATOM
            <Atom N>
            <Atom CA>
            <Atom C>
            <Atom O>
            <Atom CB>
            <Atom CG>
            <Atom SD>
            <Atom CE>
            <Atom N>
            <Atom CA>
            [...]

        '''
        aa_atoms = []
        for chain in self.structure:
            for res in chain:
                if res.id[0].upper() != 'H':
                    for atom in res:
                        aa_atoms.append(atom)

        return aa_atoms

    def get_amino_atoms_coords(self):
        '''
        Navigate through all ATOM records in the structure and return
        their coordinates.

        :returns: An iterable of the requested atom coordinates.
        :rtype: list(:class:`np.array[3]`)

        Example::

            >>> a = Structure.PDBFile('4N6W.pdb')
            >>> for ATOM in a.get_amino_atoms_coords():
            ...     print ATOM
            [-42.316  38.638  20.425]
            [-41.82   39.723  19.581]
            [-41.93   41.072  20.291]
            [-40.962  41.836  20.357]
            [-42.577  39.762  18.255]
            [-42.085  40.845  17.3  ]
            [-40.299  40.776  17.008]
            [-40.081  42.082  15.79 ]
            [-43.122  41.359  20.805]
            [-43.34   42.492  21.695]
            [...]

        '''
        aa_coord = []
        for chain in self.structure:
            for res in chain:
                if res.id[0].upper() != 'H':
                    for atom in res:
                        aa_coord.append(atom.get_coord())

        return aa_coord

    def get_all_atoms(self):
        '''
        Navigate through all ATOM and HETATM records in the structure. See
        :meth:`get_amino_atoms` and :meth:`get_hetero_atoms` to return
        all ATOM resp. HETATM entries.

        :returns: All atoms in the structure.
        :rtype: list(:class:`Atom`)

        Example::

            >>> a = Structure.PDBFile('4N6W.pdb')
            >>> for ATOM in a.get_amino_atoms_coords():
            ...     print ATOM
            <Atom N>
            <Atom CA>
            <Atom C>
            <Atom O>
            <Atom CB>
            <Atom CG>
            <Atom SD>
            <Atom CE>
            <Atom N>
            <Atom CA>
            [...]

        '''
        all_atoms = []
        for atom in self.structure.get_atoms():
            all_atoms.append(atom)
        return all_atoms

    def get_info_1(self, atoms=None):
        '''
        Get information slot 1:

        * pdb: occupancy
        * pqr: charge

        :param atoms: atom selection, or ``None`` for all atoms in the
           structure
        :type  atoms: list(:class:`Atom`)

        :returns: Slot 1.
        :rtype: list(float)
        '''
        if atoms is None:
            atoms = self.get_all_atoms()

        info_1_list = []
        for atom in atoms:
            info_1_list.append(atom.get_info_1())
        return info_1_list

    def get_info_2(self, atoms=None):
        '''
        Get information slot 2:

        * pdb: B-factor
        * pqr: atomic radius

        :param atoms: atom selection, or ``None`` for all atoms in the
           structure
        :type  atoms: list(:class:`Atom`)

        :returns: Slot 2.
        :rtype: list(float)
        '''
        if atoms is None:
            atoms = self.get_all_atoms()

        info_2_list = []
        for atom in atoms:
            info_2_list.append(atom.get_info_2())
        return info_2_list

    def get_chain_ids(self):
        '''
        Find all chains in the file and return their id ('A', 'B', ...).

        :returns: All chain ids in the structure.
        :rtype: list(str)

        Example::

            >>> a = Structure.PDBFile('4N6W.pdb')
            >>> a.get_chain_ids()
            ['A']

        '''
        chain_list = []
        for chain in self.structure:
            chain_list.append(chain.id)

        return chain_list

    def get_first_res_id(self):
        '''
        Find the number of the first residue in this structure. Useful when
        the residue numbering doesn't start at 1.

        :returns: Residue index of the first amino acid.
        :rtype: int

        Example::

            >>> a = Structure.PDBFile('4N6W.pdb')
            >>> a.get_first_res_id()
            1

        '''
        chain_id = self.structure.child_dict.keys()[0]
        return self.structure[chain_id].child_list[0].id[1]

    def get_atoms_of_type(self, atom_type):
        '''
        Find all atoms of type **atom_type** in the structure.

        :param atom_type: type of atom to find (e.g. 'CA')
        :type  atom_type: str

        :returns: An iterable of all atoms of a certain type in the structure.
        :rtype: list(:class:`Atom`).

        Example::

            >>> a = Structure.PDBFile('4N6W.pdb')
            >>> a.get_atoms_of_type('FE')
            [<Atom FE>, <Atom FE>]

        '''
        all_atoms_of_type = []
        for atom in self.structure.get_atoms():
            if atom.get_name().upper() == atom_type.upper():
                all_atoms_of_type.append(atom)
        return all_atoms_of_type


    def transform(self, T):
        '''
        Transform the structure coordinates with the given matrix.

        :param T: transformation matrix
        :type  T: :class:`np.array[3,3]`
        '''
        for atom in self.structure.get_atoms():
            atomcoord = atom.get_coord()
            newcoord = np.dot(T, atomcoord)
            atom.set_coord(newcoord)

    def translate(self, transVector):
        '''
        Translate structure coordinates by the given vector.

        :param transvector: translation vector
        :type  transvector: :class:`np.array[3]`
        '''
        for atom in self.structure.get_atoms():
            atomcoord = atom.get_coord()
            newcoord = atomcoord + np.array(transVector)
            atom.set_coord(newcoord)

    def translate_x(self, dist):
        '''
        Translate protein structure along the X axis.

        :param dist: translation distance
        :type  dist: float
        '''
        self.translate([dist, 0, 0])

    def translate_y(self, dist):
        '''
        Translate protein structure along the Y axis.

        :param dist: translation distance
        :type  dist: float
        '''
        self.translate([0, dist, 0])

    def translate_z(self, dist):
        '''
        Translate protein structure along the Z axis.

        :param dist: translation distance
        :type  dist: float
        '''
        self.translate([0, 0, dist])

    def translate_origin_and_rotate(self, phi, theta, psi):
        '''
        This methods centers the structure at the origin, rotates it with
        angle_x around the x axis (angle_y around y axis, etc.) and moves
        it back to where it was.

        :param phi: Euler angle for rotation
        :type  phi: float
        :param theta: Euler angle for rotation
        :type  theta: float
        :param psi: Euler angle for rotation
        :type  psi: float
        '''
        transvector = self.determine_geometric_center()
        self.translate(-transvector)
        self.pRotate(phi, theta, psi)
        self.translate(transvector)

    def move_to_new_position(self, new_coord):
        '''
        Move the structure geometric center to the supplied coordinates.

        :param new_coord: new coordinates
        :type  new_coord: :class:`np.array[3]`
        '''
        old_coord = self.determine_geometric_center()
        self.translate(np.array(new_coord) - np.array(old_coord))

    def rotate_and_move_to_new_position(self, phi, theta, psi, new_coord):
        '''
        Center the structure at the origin, rotate it and
        move it to the new position **new_coord**.

        :param phi: Euler angle for rotation
        :type  phi: float
        :param theta: Euler angle for rotation
        :type  theta: float
        :param psi: Euler angle for rotation
        :type  psi: float
        :param new_coord: new coordinates
        :type  new_coord: :class:`np.array[3]`
        '''
        self.translate(-self.determineGeometricCenter())
        self.pRotate(phi, theta, psi)
        self.translate(np.array(new_coord))

    def rotate(self, angle, axis):
        '''
        Rotate protein structure using the Rodrigues' rotation formula.

        :param degree: angle by which to rotate (in degrees)
        :type  degree: float
        :param axis: axis around which to rotate (x = [1,0,0],
           y = [0,1,0], z = [0,0,1])
        :type  axis: :class:`np.array[3]`

        Example::

            >>> a = Structure.PDBFile('4N6W.pdb')
            >>> print a.get_all_atom_coords()[0] # before rotation
            [-42.316,  38.638,  20.425]
            >>> a.rotate(+30.0, [1,0,0]) # x axis
            >>> print a.get_all_atom_coords()[0] # after rotation
            [-42.316  23.24898955  37.00756887]

        '''
        # TODO: check the algorithm
        axis = np.array(axis)
        vectorLength = np.sqrt(np.dot(axis, axis))
        normVector = axis / vectorLength
        # convert angle from [degree] to [rad]
        theta = np.deg2rad(angle)
        # unit matrix
        unitM = np.diag(np.ones(3))
        # cross product matrix
        crossProductmatrix = np.array([[0, -normVector[2], normVector[1]], [normVector[2], 0, -normVector[0]], [-normVector[1], normVector[0], 0]])
        # tensor product
        tensorProduct = np.zeros((3, 3))
        tensorProduct[0][0] = normVector[0] * normVector[0]
        tensorProduct[0][1] = normVector[0] * normVector[1]
        tensorProduct[0][2] = normVector[0] * normVector[2]
        tensorProduct[1][0] = normVector[1] * normVector[0]
        tensorProduct[1][1] = normVector[1] * normVector[1]
        tensorProduct[1][2] = normVector[1] * normVector[2]
        tensorProduct[2][0] = normVector[2] * normVector[0]
        tensorProduct[2][1] = normVector[2] * normVector[1]
        tensorProduct[2][2] = normVector[2] * normVector[2]
        rotationMatrix = unitM * np.cos(theta) + np.sin(theta) * crossProductmatrix + (1 - np.cos(theta)) * tensorProduct
        self.transform(rotationMatrix)

    def rotateX(self, degree):
        '''
        Rotate protein structure around the X axis.

        :param degree: angle by which to rotate (in degrees)
        :type  degree: float
        '''
        self.rotate(degree, [1, 0, 0])

    def rotateY(self, degree):
        '''
        Rotate protein structure around the Y axis.

        :param degree: angle by which to rotate (in degrees)
        :type  degree: float
        '''
        self.rotate(degree, [0, 1, 0])

    def rotateZ(self, degree):
        '''
        Rotate protein structure around the Z axis.

        :param degree: angle by which to rotate (in degrees)
        :type  degree: float
        '''
        self.rotate(degree, [0, 0, 1])

    def pRotate(self, phi, theta, psi):
        '''
        Apply Euler angle rotation to the structure.

        :param phi: Euler angle for rotation
        :type  phi: float
        :param theta: Euler angle for rotation
        :type  theta: float
        :param psi: Euler angle for rotation
        :type  psi: float
        '''
        euler_rotation = MathTools.calculate_rotation_matrix(phi, theta, psi)
        self.transform(euler_rotation)

    def rotate_by_matrix(self, rot_matrix):
        '''
        Transform the structure coordinates with the given rotation matrix.

        :param rot_matrix: transformation matrix
        :type  rot_matrix: :class:`np.array[3,3]`
        '''
        self.transform(rot_matrix)

    def determineCenterOfMass(self):
        '''
        Determine the center of mass for the protein structure.
        '''
        raise RuntimeError('not implemented ... missing atom weights!')

    def determine_geometric_center(self):
        '''
        Determine the structure geometric center.

        :returns: Center coordinates.
        :rtype: :class:`np.array[3]`

        Example::

            >>> a = Structure.PDBFile('4N6W.pdb')
            >>> print a.determine_geometric_center()
            [-27.25387547  25.13667925  13.83051887]

        '''
        all_coords = np.array(self.get_all_atom_coords())
        return np.mean(all_coords, 0)

    def determine_center_of_extremes_of_atoms(self, atoms):
        '''
        Determine the center of the extremes of an atom selection.

        :param atoms: atom selection
        :type  atoms: list(:class:`Atom`)

        :returns: Center coordinates.
        :rtype: :class:`np.array[3]`
        '''
        centerVector = np.zeros(3)
        extremes = self.determine_coordinate_extremes(atoms)
        centerVector[0] = (extremes[0][1] + extremes[0][0]) / 2.0
        centerVector[1] = (extremes[1][1] + extremes[1][0]) / 2.0
        centerVector[2] = (extremes[2][1] + extremes[2][0]) / 2.0
        return centerVector

    def determine_center_of_extremes(self):
        '''
        Determine the center of the extremes.

        :returns: Center coordinates.
        :rtype: :class:`np.array[3]`

        Example::

            >>> a = Structure.PDBFile('4N6W.pdb')
            >>> print a.determine_center_of_extremes()
            [-26.8235  24.865  13.304]

        '''
        centerVector = np.zeros(3)
        extremes = self.determine_coordinate_extremes()
        centerVector[0] = (extremes[0][1] + extremes[0][0]) / 2.0
        centerVector[1] = (extremes[1][1] + extremes[1][0]) / 2.0
        centerVector[2] = (extremes[2][1] + extremes[2][0]) / 2.0
        return centerVector

    def determine_max_diameter(self, atoms=None):
        '''
        :param atoms: atom selection, or all atoms in the structure if
           ``None``
        :type  atoms: list(:class:`Atom`)

        :returns: Maximum diameter (Angstroms).
        :rtype: float

        Example::

            >>> a = Structure.PDBFile('4N6W.pdb')
            >>> print a.determine_max_diameter()
            52.4878883362

        '''
        atom_coord = []
        if atoms is None:
            atoms = self.get_all_atoms()
            for atom in atoms:
                atom_coord.append(atom.get_coord())

        else:
            atom_coord = self.get_all_atom_coords()

        atom_coord = np.array(atom_coord)

        maxDiameter = pdb_cython.determine_max_diameter(atom_coord)
        return maxDiameter

    def determine_radius(self, atoms = None):
        '''
        :param atoms: atom selection, or all atoms in the structure if
           ``None``
        :type  atoms: list(:class:`Atom`)

        :returns: Protein radius (Angstrom).
        :rtype: float
        '''
        if atoms is None:
            atoms = self.get_all_atoms()

        radius = 0.0
        center = self.determine_goometric_center_of_atoms(atoms)

        for atom in atoms:
            current_vec = atom.get_coord() - center
            current_radius = np.linalg.norm(current_vec, current_vec)
            if current_radius > radius:
                radius = current_radius
        return radius

    def center(self):
        '''
        Translate the geometric center to the origin.
        '''
        center_vector = self.determine_geometric_center()
        self.translate(-center_vector)

    def determine_coordinate_extremes(self, atoms = None):
        '''
        :param atoms: atom selection, or all atoms in the structure if
           ``None``
        :type  atoms: list(:class:`Atom`)

        :returns: Center coordinates.
        :rtype: :class:`np.array[3,2]`
        '''
        if atoms is None:
            atoms = self.get_all_atoms()

        atom_coord = []

        for atom in atoms:
            atom_coord.append(atom.get_coord())

        atom_coord = np.array(atom_coord)
        min_vec = np.amin(atom_coord, 0)
        max_vec = np.amax(atom_coord, 0)
        return np.array(zip(min_vec, max_vec))

    def get_PCA_base(self, selection='CA'):
        '''
        Compute the unit vectors of a basis that maximizes the spread of atomic
        coordinates on the x-axis, then the y-axis, and then the z-axis.

        :param selection: subset of atoms to use for the PCA, default is
           all C-alpha atoms ('CA'), for non-proteins please choose between
           non-hydrogens ('noH') and all ('all')
        :type  selection: str

        :returns: Set of basis vectors
        :rtype: :class:`np.array[3,3]`
        '''
        if selection == 'CA':
            coord = np.array([a.coord for c in self.structure for r in c
                                      for a in r if a.name == 'CA'])
        elif selection == 'noH':
            coord = np.array([a.coord for c in self.structure for r in c
                                      for a in r if a.elem != 'H'])
        elif selection == 'all':
            coord = np.array([a.coord for c in self.structure for r in c
                                      for a in r])
        else:
            raise ValueError('Selection not understood: "{}"'.format(selection))
        return MathTools.PCA_base(coord)
    
    def apply_PCA_projection(self, base=None, selection='CA'):
        '''
        Compute atomic positions in a new basis set. Useful to center and
        rotate a protein before an APBS calculation using a rectangular box,
        so as to reduce the box dimension.

        :param base: Set of basis vectors (optional), if ``None``,
           perform a PCA on the molecule first
        :type  base: :class:`np.array[3,3]`
        :param selection: subset of atoms to use for the PCA, default is
           all C-alpha atoms ('CA'), for non-proteins please choose between
           non-hydrogens ('noH') and all ('all')
        :type  selection: str
        '''
        if base is None:
            base = self.get_PCA_base(selection=selection)
        coord = np.array(self.get_all_atom_coords())
        coord -= coord.mean(axis=0)
        coord_new = MathTools.PCA_projection(base, coord)
        i = 0
        for c in self.structure:
            for r in c:
                for a in r:
                    a.coord = coord_new[i,:]
                    i += 1

    def get_radius_of_gyration(self):
        '''
        Calculate the radius of gyration.

        :returns: Radius of gyration.
        :rtype: float

        Example::

            >>> a = Structure.PDBFile('4N6W.pdb')
            >>> a.get_radius_of_gyration()
            28.848387249760471

        '''
        all_coords = np.array(self.get_all_atom_coords())
        geo_center = np.mean(all_coords,0)
        r_g = 0

        diff = all_coords - geo_center
        r = np.sqrt(diff[:,0]**2 + diff[:,1]**2 + diff[:,2]**2)

        r_g = np.amax(r)

        return r_g

    def clone(self, chain_id_list=None, res_id_list=None, res_type_list=None,
              atom_types_list=None):
        '''
        Return a clone of self. Through the list parameters specific items can
        be selected, namely the list of residues or certain types of residues
        (mutually exclusive). One-letter codes for the residues will be
        translated to three-letter codes. Use ``None`` to select all.

        :param chain_id_list: chains to copy to the clone
        :type  chain_id_list: list(str)
        :param res_id_list: residues in each chain to copy to the clone
        :type  res_id_list: list(str)
        :param res_type_list: residue types to copy to the new clone
        :type  res_type_list: list(str)
        :param atom_types_list: atom types to copy to the new clone
        :type  atom_types_list: list(str)

        :returns: A clone of the structure object.
        :rtype: :class:`PDBFile` or :class:`PQRFile` or :class:`LatFile`
        :raises AttributeError: if both **res_id_list** and **res_type_list**
           were used, or if **self.what_am_i** is empty.

        Example::

            >>> a = Structure.PDBFile('4N6W.pdb')
            >>> b = a.clone()
            >>> print b
            <epitopsy.Structure.PDBFile object at 0x2088a90>

        '''
        if res_id_list is not None and res_type_list is not None:
            raise AttributeError("Please use eiter specific residue id's or residue types! This method has returned None!")

        # translate one letter codes to three letter codes
        if res_type_list is not None:
            for i, res in enumerate(res_type_list):
                if len(res) == 1:
                    res_type_list[i] = self.one2threeletter[res]

        if self.what_am_i == 'pdb':
            new_clone = PDBFile()
        elif self.what_am_i == 'pqr':
            new_clone = PQRFile()
        elif self.what_am_i == 'lat':
            new_clone = LatFile()
        else:
            raise AttributeError("I don't know what i am, pdb, pqr? This method returned None!")

        # this method only works with one model!
        new_structure = Structure(0)
        new_model = Model(0)
        new_structure.add(new_model)
        # iterate through structure and copy selected (or all) elements
        for chain in self.structure:
            if chain_id_list is not None and chain.id in chain_id_list:
                # chain.id is in the given list
                new_chain = Chain(chain.id)
                new_model.add(new_chain)
            elif chain_id_list is None:
                # we want to clone all chains
                new_chain = Chain(chain.id)
                new_model.add(new_chain)
            else:
                # we do not want this chain
                continue

            for res in chain:
                if res_id_list is not None and res.id[1] in res_id_list:
                    # we choose them by their id's
                    new_res = Residue(res.id, res.resname, '')
                    new_chain.add(new_res)
                elif res_type_list is not None and res.resname in res_type_list:
                    # we choose them by their type
                    new_res = Residue(res.id, res.resname, '')
                    new_chain.add(new_res)
                elif res_id_list is None and res_type_list is None:
                    # we want to clone all residues
                    new_res = Residue(res.id, res.resname, '')
                    new_chain.add(new_res)
                else:
                    # we do not want this residue
                    continue

                for atom in res:
                    if(atom_types_list is not None
                       and atom.id in atom_types_list):
                        new_atom = Atom(atom.name, atom.get_coord(),
                                        atom.get_info_1(), atom.get_info_2(),
                                        ' ', atom.fullname, atom.serial_number)
                        if self.what_am_i == 'pdb':
                            new_atom.set_element(atom.element)
                            new_atom.set_charge(atom.charge)
                        new_res.add(new_atom)
                    elif atom_types_list is None:
                        new_atom = Atom(atom.name, atom.get_coord(),
                                        atom.get_info_1(), atom.get_info_2(),
                                        ' ', atom.fullname, atom.serial_number)
                        if self.what_am_i == 'pdb':
                            new_atom.set_element(atom.element)
                            new_atom.set_charge(atom.charge)
                        new_res.add(new_atom)
                    else:
                        # we do not want this atom
                        continue

        new_clone.structure = new_structure[0]

        return new_clone

    def get_residue_id_list(self, chain_id=None):
        '''
        Display all residue numbers found in the structure.

        :param chain_id: which chain in the structure should be used
           (optional), if ``None`` use all available chains
        :type  chain_id: str

        :returns: a list with all residue ID's of the structure

        Example::

            >>> a = Structure.PDBFile('4N6W.pdb')
            >>> print a.get_residue_id_list()
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, [...], 444]

        '''
        res_id_list = []
        if chain_id is None:
            # if no chain id has been supplied take all chains
            for chain in self.structure:
                for res in chain:
                    res_id_list.append(res.get_id()[1])
        else:
            for res in self.structure[chain_id]:
                res_id_list.append(res.get_id()[1])

        return res_id_list

    def get_res_id_aa_dict(self, chain_id):
        '''
        Display all residues from chain **chain_id** in a dictionary with
        residue number as key and amino acid one-letter code as value. The
        advantage of a dictionary over a list is that gaps in the sequence
        numbering are preserved. All non-amino acids are ignored.

        :param chain_id: which chain in the PDB file should be used
        :type  chain_id: str

        :returns: Residues number / one-letter code pairs.
        :rtype: dict(int,str)

        Example::

            >>> a = Structure.PDBFile('4N6W.pdb')
            >>> a.get_res_id_aa_dict('A')
            Encountered the following non amino acids: ['FE', 'FLC', 'HOH']
            {1: 'M', 2: 'S', 3: 'L', 4: 'S', 5: 'N', 6: 'S', 7: 'S', 8: 'K', 9: 'V', 10: 'S', [...],  187: 'E'}

        '''
        res_map = {}
        non_amino_acid = False
        non_amino_acid_list = []
        for res in self.structure[chain_id]:
            res_id = res.get_id()[1]
            res_name = res.resname
            # only add the amino acid, if it is one of proteogen amino acids
            # in that way no calciums etc. will be added
            if res_name in self.three2oneletter:
                aa = self.three2oneletter[res_name]
                res_map[res_id] = aa
            else:
                non_amino_acid = True
                non_amino_acid_list.append(res_name)

        if non_amino_acid is True:
            print('Encountered the following non amino acids: {0}'
                  .format(sorted(set(non_amino_acid_list))))

        return res_map

    def contains_chain_break(self, chain_id=None):
        '''
        Tell if a chain break exist in **chain_id**, or in the whole structure
        if omitted. HETATM are skipped.

        :param chain_id: the chain id where to look for a break (optional),
           uses all available chains if ``None``
        :type  chain_id: str

        :returns: ``False`` if no chain break was detected, or a formatted
           string if a break was detected.
        :rtype: bool or str

        '''
        if chain_id is None:
            # if no chain id has been supplied take all chains
            for chain in self.structure:
                res_id_list = []
                for i,res in enumerate(chain):
                    if res.id[0] != "H":
                        res_id = res.get_id()[1]
                        if i > 0 and res_id_list:
                            if (res_id_list[-1] + 1) != res_id:
                                    return 'Break between {0} and {1}'.format(
                                                      res_id_list[-1], res_id)
                        res_id_list.append(res_id)
        else:
            res_id_list = []
            for i,res in enumerate(self.structure[chain_id]):
                if res.id[0] != "H":
                    res_id = res.get_id()[1]
                    if i > 0:
                        if res_id_list[-1] + 1 != res_id:
                                return 'Break between {0} and {1}'.format(
                                                  res_id_list[-1], res_id)
                    res_id_list.append(res_id)

        ## no chain break
        return False


    def get_res_id_array_mapping(self):
        '''
        Remove gaps in the residue sequence and return the new mapping
        in a dictionary, with the old residue id's as key and the new
        id's as value. The new mapping starts at zero and is suitable
        for use as an index for an array.

        :returns: Mapping vector.
        :rtype: dict(int,int)

        Example::

            >>> a = Structure.PDBFile('4N6W.pdb')
            >>> a.get_res_id_array_mapping()
            {1: 0, 2: 1, 3: 2, 4: 3, 5: 4, 6: 5, 7: 6, 8: 7, 9: 8, 10: 9, ..., 444: 326}

        '''
        id_dict = {}
        chain = self.structure.child_list[0].id
        for i,res in enumerate(self.structure[chain]):
            id_dict[res.id[1]] = i

        return id_dict

    def get_residue_names_from_res_id_list(self, res_id_list, chain_id=None):
        '''
        Display the three-letter code of the residues given in *res_id_list*
        from chain *chain_id* (if ``None``, take the first chain).

        :param res_id_list: residue id's from which one wants the names
        :type  res_id_list: list(int)
        :param chain_id: if ``None`` uses the all available chains (optional)
        :type  chain_id: str

        :returns: Three-letter code of residues matching the given criteria.
        :rtype: list(str)

        Example::

            >>> a = Structure.PDBFile('4N6W.pdb')
            >>> a.get_residue_names_from_res_id_list([1,2])
            ['MET', 'SER']

        '''
        if chain_id is None:
            # if no chain id has been supplied take the first one
            chain_id = self.structure.child_dict.keys()[0]

        res_name = []
        for res in self.structure[chain_id].child_list:
            if res.get_id()[1] in res_id_list:
                res_name.append(res.resname)

        return res_name

    def get_residues(self, chain_id=None, res_id_list=None):
        '''
        Display all residues from **res_id_list** contained in chain
        **chain_id** as :class:`Bio.PDB.Residue`. If **res_id_list** is
        ``None``, display all residues from chain **chain_id**.
        If **chain_id** is ``None``, return an error excepted when there is
        only one chain in the PDB file.

        :param chain_id: if ``None``, it uses the all available chains (optional)
        :type  chain_id: str
        :param res_id_list: residue id's from which one wants the objects
        :type  res_id_list: list

        :returns: a list of residues matching the given criteria
        :rtype: list(:class:`Bio.PDB.Residue`)

        Example::

            >>> a = Structure.PDBFile('4N6W.pdb')
            >>> a.get_residues('A')[0:2] # residue id starts at 0 in a list
            [<Residue MET het=  resseq=1 icode= >, <Residue SER het=  resseq=2 icode= >]

        '''
        if chain_id is None:
            if len(self.structure.child_list) > 1:
                raise ValueError("This structure has more than one chain: '{0}".format(self.structure.child_list))
            else:
                chain_id = self.structure.child_list[0].id

        if res_id_list is None:
            res_list = self.structure[chain_id].child_list
            return res_list

        res_list = []
        for res in self.structure[chain_id].child_list:
            if res.get_id()[1] in res_id_list:
                res_list.append(res)

        return res_list

    def get_atoms_by_id(self, atom_id_list):
        '''
        Get atoms by their corresponding number from the pdb-file.

        :param atom_id_list: atom id numbers
        :type  atom_id_list: list(int)

        :returns: Qtoms matching the given criteria.
        :rtype: list(:class:`Atom`)

        Example::

            >>> a = Structure.PDBFile('4N6W.pdb')
            >>> a.get_atoms_by_id([1,2,3])
            [<Atom N>, <Atom CA>, <Atom C>]

        '''
        atom_list = []
        for atom in self.get_all_atoms():
            if atom.serial_number in atom_id_list:
                atom_list.append(atom)
        return atom_list

    def get_atoms_close_to_reference(self, reference_point, max_radius,
                                     min_radius=0, atoms=None):
        '''
        Find all atoms lying in the range [**min_radius**, **max_radius**] to
        **reference_point**.

        :param reference_point: coordinates of the reference point
        :type  reference_point: :class:`np.array[3]`
        :param max_radius: maximal distance from the reference (Angstroms)
        :type  max_radius: float
        :param min_radius: minimal distance from the reference (Angstroms)
        :type  min_radius: float
        :param atoms: list of the atoms to which the results should be
           restricted (optional), if ``None`` it uses all
           atoms from the protein
        :type  atoms: list(:class:`Atom`)

        :returns: Atoms close to the given reference point.
        :rtype: list(:class:`Atom`)

        Example::

            >>> a = Structure.PDBFile('4N6W.pdb')
            >>> # retrieve the iron catalytic center
            >>> fe_list = a.get_atoms_of_type('FE')
            >>> coord = a.get_coords_from_atom_list([fe_list[1]])
            >>> # get all 6 atoms chelating Fe, excepted Fe
            >>> a.get_atoms_close_to_reference(coord, 2.5, 0.1)
            [<Atom NE2>, <Atom NE2>, <Atom OD2>, <Atom OD1>, <Atom OG1>, ...]

        Biopython equivalent::

            >>> from Bio.PDB import PDBParser
            >>> from Bio.PDB import NeighborSearch
            >>> p = PDBParser()
            >>> a = p.get_structure('4N6W.pdb', '4N6W.pdb')
            >>> fe_list = []
            >>> atom_list = []
            >>> for residue in struct[0]['A']:
            ...     for atom in residue:
            ...         if atom.name == 'FE':
            ...             fe_list.append(atom)
            ...         elif atom.name[0] <> 'H':
            ...             atom_list.append(atom)
            >>> ns = NeighborSearch(atom_list)
            >>> rd = 2.5
            >>> coord = fe_list[1].get_coord()
            >>> print sorted(ns.search(coord, rd, 'A'))
            [<Atom NE2>, <Atom NE2>, <Atom OD2>, <Atom OD1>, <Atom OG1>, ...]

        '''
        if atoms is None:
            atoms = self.get_all_atoms()

        reference_point = np.array(reference_point)

        nearby_atoms = []

        for atom in atoms:
            atom_coord = atom.get_coord()
            diff_vec = atom_coord - reference_point
            distance = np.linalg.norm(diff_vec)
            if min_radius <= distance <= max_radius:
                nearby_atoms.append(atom)

        return nearby_atoms

    def find_chain_contacts(self, chain1, chain2, max_distance):
        '''
        Find atoms of **chain1** which are within **max_distance** of
        **chain2**.

        :param chain1: first chain id
        :type  chain1: str
        :param chain2: second chain id
        :type  chain2: str
        :param max_distance: maximal distance to include atoms in the search
        :type  max_distance: float

        :returns: Atoms from **chain1** for which distance to
           **chain2** is smaller than **max_distance**.
        :rtype: list(:class:`Atom`)

        '''
        contact_atoms = []
        # get chains:
        chain_a_found = False
        chain_b_found = False
        for chain in self.structure:
            if chain.id == chain1:
                chain_a = chain
                chain_a_found = True
            elif chain.id == chain2:
                chain_b = chain
                chain_b_found = True
        if chain_b_found == False or chain_a_found == False:
            print("Could not find chain {0} and/or chain {1}!.".format(chain1, chain2))
            return
        atoms_a = chain_a.get_atoms()
        atoms_b = chain_b.get_atoms()
        for atom_a in atoms_a:
            for atom_b in atoms_b:
                diff_vec = atom_a.get_coord() - atom_b.get_coord()
                if np.linalg.norm(diff_vec) <= max_distance:
                    contact_atoms.append(atom_a)
        return contact_atoms

    def get_sequence(self):
        '''
        Display the one-letter code sequence of each chain of the PDB file in
        a dictionary. Non standard amino acids will not be returned.

        :returns: Amino acid sequences by chains.
        :rtype: dict(str,str)

        Example::

            >>> a = Structure.PDBFile('4N6W.pdb')
            >>> a.get_sequence()
            {'A': 'MSL[...]LSE'}

        '''
        seq_dict = {}
        for chain in self.structure:
            chain_id = chain.id
            for res in chain:
                res_name = res.resname
                if res_name not in self.three2oneletter:
                    continue
                else:
                    if chain_id in seq_dict:
                        seq_dict[chain_id] = seq_dict[chain_id] + self.three2oneletter[res_name]
                    else:
                        seq_dict[chain_id] = self.three2oneletter[res_name]

        return seq_dict

    def snap_vdw_to_box(self, box_mesh_size, box_dim, box_offset,
                        warning=True, vdw_radii=2, increase_vdw_by=0):
        '''
        This method snaps a structure to a given dxbox. If the structure is
        a pqr it uses the supplied vdw radii otherwise it uses the variable
        'vdw_radii' for each atom.
        
        If any coordinate lies outside the box an error will be printed to the
        standard output.

        :param box_mesh_size: mesh size
        :type  box_mesh_size: :class:`np.array[3]`
        :param box_dim: box dimesions
        :type  box_dim: :class:`np.array[3]`
        :param box_offset: box origin
        :type  box_offset: :class:`np.array[3]`
        :param warning: print an error if the structure does not fit
           completely into the given box dimensions
        :type  warning: bool
        :param vdw_radii: default vdw radius if not a PQR file (Angstroms)
        :type  vdw_radii: float
        :param increase_vdw_by: increase all vdw radii by a constant (Angstroms)
        :type  increase_vdw_by: float

        :returns: OpenDX box with 0's outside the molecule and 1's inside
        :returntype: :class:`np.array[:,:,:]`
        '''
        vdw_box = np.zeros(box_dim)
        mesh_size = float(box_mesh_size[0])
        if not (box_mesh_size[0] == box_mesh_size[1] == box_mesh_size[2]):
            print ('WARNING: mesh size different depending on the direction, '
                   'using only the x-axis to snap vdw radii to grid space')

        # list which stores the coordinates of all atoms
        coordinates = []

        # list which contains vdw_radii
        vdw_radii_list = []

        # respect the mesh size!
        vdw_radii = vdw_radii / mesh_size

        # we need a dxbox, because it lets us transform the coordinates
        dxbox = DXBox(vdw_box, box_mesh_size, box_offset)

        # error warning
        error_warning = False

        for atom in self.structure.get_atoms():
            box_coord = list(dxbox.transform_real_to_box_space(list(atom.get_coord())))
            # we do not want to leave the box!
            if(0 <= box_coord[0] < box_dim[0] and 0 <= box_coord[1] < box_dim[1]
               and 0 <= box_coord[2] < box_dim[2]):
                # check if coordinates already in the list
                if box_coord not in coordinates:
                    coordinates.append(box_coord)

                    if self.what_am_i == 'pqr':
                        # it is a pqr and there are vdw radii information
                        new_radii = (atom.get_info_2() + increase_vdw_by) / mesh_size
                        if new_radii == 0.0:
                            # for example hydrogen sometimes has 0.0, in
                            # sepecial cases we need to make sure, that
                            # there will be a grid point
                            # the value just needs to be smaller then the
                            # mesh size
                            new_radii = mesh_size / 10.

                        vdw_radii_list.append(new_radii)
                    else:
                        # not a pqr and therfore no vdw radii information
                        vdw_radii_list.append((vdw_radii + increase_vdw_by))
                else:
                    position_index = coordinates.index(box_coord)

                    if self.what_am_i == 'pqr':
                        # it is a pqr and there are vdw radii information
                        if (vdw_radii_list[position_index] <
                            (atom.get_info_2() + increase_vdw_by) / mesh_size):
                            # take the bigger value
                            vdw_radii_list[position_index] = (atom.get_info_2() + increase_vdw_by) / mesh_size
                        else:
                            # the stored one is already bigger
                            pass
                    else:
                        # not a pqr and therfore no vdw radii information
                        vdw_radii_list[position_index] = vdw_radii

            else:
                error_warning = True

        if warning is True and error_warning is True:
            print('{0} ({1}) has left the box!'.format(self.what_am_i,
                                                       self.structure_path))

        coordinates = np.array(coordinates)
        # set vdw radii
        vdw_box[(coordinates[:][:, 0], coordinates[:][:, 1], coordinates[:][:, 2])] = vdw_radii_list

        vdw_box = pdb_cython.build_vdw_surface(vdw_box)
        #vdw_box[np.nonzero(vdw_box != 0)] = 1

        return vdw_box

    def get_rmsd_rotation_translation_from_superposition(self, pdb_to_rotate,
                                                         atom_types=None):
        '''
        Fit the structure **pdb_to_rotate** onto self.
        Return a dictionary containing a *rmsd* value, a *rotation* matrix
        and a *translation* value.
        The parameter **atom_types** restricts the fitting on a particular set
        of atoms, if ``None`` is supplied, all atoms will be used.

        I guess the units are Angstroem.

        :param pdb_to_rotate: object to superimpose onto this object
        :type  pdb_to_rotate: :class:`Structure`
        :param atom_types: restrict the fitting to specific types (e.g.
           ``['CA']``, ``['CA', 'N']``), or fit all atoms if ``None``
        :type  atom_types: list(str)

        :returns: A dictionary with following keys: ``'rmsd'``: root mean
           square deviation, ``'rotation'``: rotation matrix,
           ``'translation'``: translation vector
        :rtype: dict
        :raises ValueError: if there is an atom mismatch between the two pdb's
           or if at least one of the atom lists is empty
        '''
        # check if atom_types is a list
        if not isinstance(atom_types, list) and atom_types is not None:
            atom_types = [atom_types]

        if atom_types is None:
            my_atom_coords = self.get_all_atom_coords()
            ref_pdb_coords = pdb_to_rotate.get_all_atom_coords()
        else:
            my_atoms = []
            ref_pdb_atoms = []
            # add all given types
            for atom_type in atom_types:
                my_atoms.extend(self.get_atoms_of_type(atom_type))
                ref_pdb_atoms.extend(pdb_to_rotate.get_atoms_of_type(atom_type))

            # iterate over atoms to add coords
            my_atom_coords = []
            ref_pdb_coords = []
            for atom in my_atoms:
                my_atom_coords.append(atom.get_coord())

            for atom in ref_pdb_atoms:
                ref_pdb_coords.append(atom.get_coord())

        # check for empty lists
        if len(my_atom_coords) == 0 or len(ref_pdb_coords) == 0:
            raise ValueError("At least one of the atom lists is empty!")

        if len(my_atom_coords) != len(ref_pdb_coords):
            raise ValueError("Atom missmatch between the two pdb's!")

        result_dict = pdb_cython.get_rmsd(np.array(my_atom_coords),
                                          np.array(ref_pdb_coords))

        return result_dict

    def superimpose_given_pdb_onto_self(self, pdb_to_superimpose,
                                        atom_types=None):
        '''
        Fit the structure **pdb_to_superimpose** onto self.
        The atomic coordinates are updated in the process.

        :param pdb_to_superimpose: object to superimpose onto this object
        :type  pdb_to_superimpose: :class:`Structure`
        :param atom_types: restrict the fitting to specific types (e.g.
           ``['CA']``, ``['CA', 'N']``), or fit all atoms if ``None``
        :type  atom_types: list(str)

        :raises ValueError: if there is an atom mismatch between the two pdb's
           or if at least one of the atom lists is empty
        '''
        super_instructions = self.get_rmsd_rotation_translation_from_superposition(pdb_to_superimpose, atom_types)
        trans_vector = super_instructions['translation']
        rot_matrix = super_instructions['rotation']

        for atom in pdb_to_superimpose.structure.get_atoms():
            atom.transform(rot_matrix, trans_vector)

    def superimpose_self_onto_given_pdb(self, pdb_to_superimpose,
                                        atom_types=None):
        '''
        Fit self onto the structure **pdb_to_superimpose** and update self's
        atomic coordinates.

        :param pdb_to_superimpose: object on which to superimpose self
        :type  pdb_to_superimpose: :class:`Structure`
        :param atom_types: restrict the fitting to specific types (e.g.
           ``['CA']``, ``['CA', 'N']``), or fit all atoms if ``None``
        :type  atom_types: list(str)

        '''
        super_instructions = pdb_to_superimpose.get_rmsd_rotation_translation_from_superposition(self, atom_types)
        trans_vector = super_instructions['translation']
        rot_matrix = super_instructions['rotation']

        for atom in self.structure.get_atoms():
            atom.transform(rot_matrix, trans_vector)

    def get_dxbox_dim(self, box_mesh_size, extend=None, cubic_box=False,
                      nlev=4):
        '''
        Return the dimensions of a DXbox. The edges are calculated using the
        protein maximal diameter in each direction, **extend** if given,
        and the grid resolution **box_mesh_size**, using the formula:

        :math:`a[i] = \\frac{\\displaystyle protein\\_diameter[i] +
        2 \\cdot extend[i]}{\\displaystyle box\\_mesh\\_size[i]}`

        with *i* = {x,y,z} the coordinates. If **cubic_box** is ``True``,
        all edges have the same length.

        :param box_mesh_size: grid resolution (Angstroms)
        :type  box_mesh_size: :class:`np.array[3]`
        :param extend: extension of the box dimensions (Angstroms)
        :type  extend: float
        :param cubix_box: use a cubic box if ``True``
        :type  cubic_box: bool
        :param nlev: depth of the multilevel hierarchy used by the multigrid
           solver (optional)
        :type  nlev: int

        :returns: Dimensions of the DXbox.
        :rtype: :class:`np.array[3]`
        '''
        box_mesh_size = np.array(box_mesh_size)
        if cubic_box is True:
            max_dim = self.determine_max_diameter()
            box_dim = np.array([max_dim, max_dim, max_dim])
        else:
            extremes = self.determine_coordinate_extremes()
            box_dim = np.zeros(3)
            box_dim[0] = abs((extremes[0][1] - extremes[0][0]))
            box_dim[1] = abs((extremes[1][1] - extremes[1][0]))
            box_dim[2] = abs((extremes[2][1] - extremes[2][0]))

        if extend is not None:
            # add extend in each direction
            box_dim = box_dim + 2 * float(extend)

        box_dim = np.ceil(box_dim / box_mesh_size)
        from epitopsy.APBS import fix_grid_size
        box_dim = fix_grid_size(box_dim, nlev)

        return box_dim

    def get_dxbox_offset(self, box_mesh_size, box_dim, box_center):
        '''
        Returns the offset for the given dimensions.

        :param box_mesh_size: dimensions of the mesh (Angstroms)
        :type  box_mesh_size: :class:`np.array[3]`
        :param box_dim: dimensions of the box (Angstroms)
        :type  box_dim: :class:`np.array[3]`
        :param box_center: center of the box (Angstroms)
        :type  box_center: :class:`np.array[3]`

        :returns: The box offset.
        :rtype: :class:`np.array[3]`
        '''
        box_offset = np.array(box_center) - (np.array(box_mesh_size) *
                                             np.array(box_dim) // 2)
        return box_offset

    def get_hydrophobic_potential(self, box_mesh_size, box_dim, box_offset):
        '''
        Calculate the hydrophobic potential of a protein.
        This method uses a simplified model for the hydrophobic potential.
        The charges are taken from the Kyte and Doolittle hydrophobicity
        scale. For each residue the center is calculated and the potential
        is modelled as:

            :math:`\\phi = \\sum \\left ( \\text{hydrophobic\\_score} \
                                          \\times e^{-\\text{distance}} \
                                 \\right )`

        .. seealso::

            Kyte, Doolittle, *A simple method for displaying the hydropathic
            character of a protein*, *J. Mol. Biol.* **1982**, *157*, 105-132.

        :param box_mesh_size: grid resolution (Angstroms)
        :type  box_mesh_size: :class:`np.array[3]`
        :param box_dim: dimensions of the box (Angstroms)
        :type  box_dim: :class:`np.array[3]`
        :param box_offset: offset of the box (Angstroms)
        :type  box_offset: :class:`np.array[3]`

        :returns: Hydrophobic potential grid.
        :rtype: :class:`np.array[:,:,:]`
        '''
        start_time = time.time()
        def get_residue_center(res):
            coord_list = []
            for atom in res:
                coord_list.append(atom.get_coord())

            return np.mean(coord_list,0)

        all_residues = self.get_residues()
        res_center_list = []
        res_hydro_list = []
        for res in all_residues:
            aa = self.three2oneletter[res.resname]
            hydro_charge = self.kyte_doolittle[aa]
            center = get_residue_center(res)
            res_center_list.append(center)
            res_hydro_list.append(hydro_charge)

        res_center_list = np.array(res_center_list)
        res_hydro_list = np.array(res_hydro_list)

        potential = pdb_cython.get_hydrophic_potential(res_center_list,
                res_hydro_list,
                box_mesh_size,
                box_dim,
                box_offset)

        end_time = time.time()
        #print('it took {0} s'.format(end_time - start_time))
        return potential

    def get_vdw_hull(self, box_mesh_size, box_dim, box_offset, vdw_radii = 2,
                     increase_vdw_by=0):
        '''
        Get the protein van der Waals hull.

        :param box_mesh_size: grid resolution (Angstroms)
        :type  box_mesh_size: :class:`np.array[3]`
        :param box_dim: dimensions of the box (Angstroms)
        :type  box_dim: :class:`np.array[3]`
        :param box_offset: offset of the box (Angstroms)
        :type  box_offset: :class:`np.array[3]`
        :param vdw_radii: default atom radii, if this is not a pqr
        :type  vdw_radii: float
        :param increase_vdw_by: extend each radius by this value
        :type  increase_vdw_by: float

        :returns: Grid with 1's at the hull and 0's everywhere else.
        :rtype: :class:`np.array[:,:,:]`
        '''
        # vdw_array has 1's inside and 0's outside
        vdw_array = self.snap_vdw_to_box(box_mesh_size, box_dim, box_offset,
                                         vdw_radii, increase_vdw_by)
        # a normal vdwbox has 0's inside and 1's outside
        vdwbox = VDWBox(vdw_array, box_mesh_size, box_offset)
        # now: 2's inside and 1's outside
        vdwbox.box = vdwbox.box + 1
        inside_pos = np.nonzero(vdwbox.box == 2)
        # now it is correct: 0's inside and 1's outside
        vdwbox.box[inside_pos] = 0

        # this also floods the box
        vdwbox.calculate_sas()
        vdwbox.prepare_for_geometric_matching(0)

        return vdwbox.box

    def get_num_of_overlap_atoms_with_given_structure(self,
                                                      other_structure_object,
                                                      energy_cutoff=1.,
                                                      vdw_radii=2.):
        '''
        Calculate the number of atoms in this structure object that are
        overlapping with the given structure object.

        :param other_structure_object: structure which might overlap with
           this one
        :type  other_structure_object: :class:`Structure`
        :param energy_cutoff: cutoff for the Lennard-Jones potential to decide
           if there is an overlap or not
        :type  energy_cutoff: float
        :param vdw_radii: default atomic radius for PDB structures
        :type  vdw_radii: float

        :returns: Number of atoms in this structure that overlap with
           the given structure object
        :rtype: int

        '''
        self_coords = self.get_all_atom_coords()
        other_coords = other_structure_object.get_all_atom_coords()

        if self.what_am_i == 'pqr':
            self_radii = self.get_info_2()
        else:
            self_radii = [vdw_radii for i in range(len(self_coords))]

        if other_structure_object.what_am_i == 'pqr':
            other_radii = other_structure_object.get_info_2()
        else:
            other_radii = [vdw_radii for i in range(len(other_coords))]

        # array them
        self_coords = np.array(self_coords, dtype=np.float)
        self_radii = np.array(self_radii, dtype=np.float)
        other_coords = np.array(other_coords, dtype=np.float)
        other_radii = np.array(other_radii, dtype=np.float)

        return pdb_cython.find_clashing_atoms(self_coords, self_radii,
                other_coords, other_radii, float(energy_cutoff))

    def get_contact_list(self, cutoff=5.):
        '''
        Find steric contacts between residues of a protein. A contact is found
        when the interatomic distance of at least two atoms taken from two
        different residues *i* and *j* is inferior to **cutoff** (Angstroem).

        :param cutoff: interatomic distance
        :type  cutoff: float

        :returns: Residues id's as found in the PDB, with 1 for a contact
           and 0 for no contact: ``[[i-1,j-1,0], [i,j,1], [i+1,j+1,0]]``.
        :rtype: list

        Example::

            >>> a = Structure.PDBFile('4N6W.pdb')
            >>> for contact in a.get_contact_list(cutoff=5.0):
            ...     print contact
            [1, 2, 1]
            [1, 3, 1]
            [1, 4, 0]
            [1, 5, 0]
            [1, 6, 0]
            [1, 7, 0]
            [1, 8, 1]
            [1, 9, 0]
            [1, 10, 0]
            [...]

        '''
        def residue_contact(res_i, res_j, cutoff):
            for atom_i in res_i:
                coord_i = atom_i.get_coord()
                for atom_j in res_j:
                    coord_j = atom_j.get_coord()
                    dist = np.linalg.norm(coord_i - coord_j)
                    if dist <= cutoff:
                        return True

            return False

        chain = self.structure.child_list[0].id
        n_res = len(self.structure[chain].child_list)
        contact_list = []
        for i in range(n_res):
            res_i = self.structure[chain].child_list[i]
            for j in range(i+1,n_res):
                res_j = self.structure[chain].child_list[j]

                if residue_contact(res_i, res_j, cutoff) is True:
                    contact_list.append([res_i.id[1], res_j.id[1], 1])
                else:
                    contact_list.append([res_i.id[1], res_j.id[1], 0])

        return contact_list


class PDBFile(Structure_Template):
    '''
    Container for PDB files.
    '''

    def __init__(self, pdb_path=None):
        Structure_Template.__init__(self, pdb_path)
        # read structure
        if self.structure_path is None:
            self.structure = None
        else:
            self._read_file()
        # i am a 'pdb'
        self.what_am_i = 'pdb'

    def get_pqr_structure(self, new_pqr_path=None, force_field='amber',
                          pdb2pqr_argv=None, pdb2pqr_path='pdb2pqr',
                          add_ions=True, add_chain=True):
        '''
        Call pdb2pqr to replace the B-factor and the occupancy information in
        the pdb file with the charges and the vdw radii. If this object has no
        ``structure_path`` property (i.e. it is ``None``), then an error is
        raised.

        If the pdb contains SO\ :sub:`4`:sup:`2--`, Ca\ :sup:`2+` or
        Zn\ :sup:`2+` ions, they will be added unless stated otherwise.

        :param new_pqr_path: path for the new pqr file. If ``None`` is given,
           **structure_path** will be used by changing its extension to .pqr
        :type  new_pqr_path: str
        :param force_field: forcefield from which charges and radii should be
           taken. Default is 'amber'
        :type  force_field: str
        :param pdb2pqr_argv: can contain additional arguments to pdb2pqr
           (e.g. ``['--assign-only']`` or ``['--noopt']``) (optional)
        :type  pdb2pqr_argv: list
        :param pdb2pqr_path: path to the executable of PDB2PQR
        :type  pdb2pqr_path: str
        :param add_ions: add Ca, Zn, SO\ :sub:`4` ions if they are in the pdb.
           The Ca, Zn and SO\ :sub:`4` atoms should have a residue name that
           fits their type (CA, ZN, SO4).
        :type  add_ions: bool

        :returns: A PQR structure.
        :rtype: :class:`PQRFile`
        :raises AttributeError: if **structure_path** is empty, or if
           **new_pqr_path** cannot be read or generated from the PDB path
        :raises NameError: if **new_pqr_path** refers to an existing file

        '''
        if self.structure_path is None:
            raise AttributeError("This pdb object has no structure path!")
        else:
            if not os.path.exists(self.structure_path):
                self.save_to_file(self.structure_path)

        if new_pqr_path is None:
            if self.structure_path.endswith(".pdb"):
                new_pqr_path = self.structure_path.replace(".pdb", ".pqr")
            else:
                raise AttributeError("This pdb object does not end with '.pdb' ('{0}') and no new pqr path is supplied!".format(self.structure_path))


        ## --chain keeps the chain id in the pqr-file
        cmdlist = [pdb2pqr_path, "--ff={0}".format(force_field)]

        if add_chain is True:
            cmdlist.append("--chain")

        if pdb2pqr_argv is not None:
            cmdlist.extend(pdb2pqr_argv)

        ## append filenames
        cmdlist.extend([self.structure_path, new_pqr_path])

        ## open file to redirect the output
        pdb2pqr_file = '/dev/zero'
        with open(pdb2pqr_file, 'w') as f:
            subprocess.call(cmdlist, stdout = f, stderr = f)

        time.sleep(0.2)
        if not(os.path.exists(new_pqr_path)):
            raise NameError('{0} could not be created'.format(new_pqr_path))

        pqr = PQRFile(new_pqr_path)

        ## add ions
        if add_ions is True:
            # calcium
            name_ca = 'CA'
            charge_ca = 2.0
            radius_ca = 1.76
            # zinc
            name_zn = 'ZN'
            charge_zn = 2.0
            radius_zn = 1.24
            # SO4
            name_so4 = 'SO4'
            radius_s = 1.85
            charge_s = 2.241
            # bound to sulfor
            radius_o_s = 1.40
            charge_o_s_2 = -1.294
            charge_o_s = -0.826
            ## search for ions
            pqr_chain_ids = [x.id for x in pqr.structure]
            for chain in self.structure:
                if chain.id not in pqr_chain_ids:
                    pqr.structure.add(Chain(chain.id))
                for res in chain:
                    if res.resname.strip() == name_ca:
                        for atom in res:
                            atom.set_info_1(charge_ca)
                            atom.set_info_2(radius_ca)
                        pqr.structure[chain.id].add(res)
                    if res.resname.strip() == name_zn:
                        for atom in res:
                            atom.set_info_1(charge_zn)
                            atom.set_info_2(radius_zn)
                        pqr.structure[chain.id].add(res)
                    if res.resname.strip() == name_so4:
                        count_oxygen = 0
                        for atom in res:
                            if 'O' in atom.name.strip():
                                count_oxygen += 1
                                atom.set_info_2(radius_o_s)
                                if count_oxygen > 2:
                                    atom.set_info_1(charge_o_s)
                                else:
                                    atom.set_info_1(charge_o_s_2)
                            elif 'S' in atom.name.strip():
                                atom.set_info_1(charge_s)
                                atom.set_info_2(radius_s)
                        pqr.structure[chain.id].add(res)

            ## save the file and reload the structure
            pqr.save_to_file(pqr.structure_path)
            pqr = PQRFile(pqr.structure_path)

        return pqr


    def _read_file(self):
        '''
        This method reads the data from the given file.
        '''
        alt_location_status = False
        if self.structure_path.endswith('.pdb'):
            with open(self.structure_path) as f:
                content = f.readlines()
        elif self.structure_path.endswith(".pdb.gz"):
            with gzip.open(self.structure_path) as f:
                content = f.readlines()
        else:
            raise AttributeError("Path does not end with either '.pdb' or '.pdb.gz'\n{0}".format(self.structure_path))

        check_chain_id = None
        check_res_id = None
        check_atom_id = None

        # first 'model' -> 0
        structure = Structure(0)
        model_counter = 0
        model = Model(model_counter)
        structure.add(model)
        for line in content:
            ## debug!!!
            if len(line) > 81:
                raise AttributeError("Wrong formatted pdb files, line is longer than 81 characters!\n{0}".format(line))
            if (line[0:4] != 'ATOM' and line[0:6] != 'HETATM' and
                    line[0:3] != 'TER' and line[0:3] != 'END' and
                    line[0:5] != 'MODEL' and line[0:6] != 'ANISOU' and
                    line[0:6] != 'ENDMDL'):
                self.remarks.append(line)
            if line[0:5] == 'MODEL':
                model_counter = +1
            if model_counter == 2:
                raise ValueError('This method is not supposed to work with different models! Please remove them manually!')
            if line[0:4] == 'ATOM' or line[0:6] == 'HETATM':
                line = line.strip('\n')
                atom_id = int(line[6:11])
                atom_name = line[12:16].strip()
                full_atom_name = line[12:16]
                alt_location = line[16:17]
                res_id = int(line[22:26])
                residue_name = line[17:20]
                chain_id = line[21:22]
                x_coord = float(line[30:38])
                y_coord = float(line[38:46])
                z_coord = float(line[47:54])
                atom_coord = np.array([x_coord, y_coord, z_coord])
                try:
                    # can not convert empty string ('  ') to float
                    # and because this information is not so important
                    # we can use this hack to keep it working
                    info_1 = float(line[54:61])
                except:
                    info_1 = 0
                try:
                    # can not convert empty string ('  ') to float
                    # and because this information is not so important
                    # we can use this hack to keep it working
                    info_2 = float(line[61:70])
                except:
                    info_2 = 0
                element = line[76:79]
                charge = line[79:80]

                if line[0:4] == 'ATOM':
                    het_flag = ' '
                elif line[0:6] == 'HETATM':
                    het_flag = 'H'

                # biopython artifact
                res_id = (het_flag, res_id, ' ')

                # if this pdb contains alternative locations, it will
                # print out a warning, but it only takes the A-position
                if alt_location != ' ':
                    alt_location_status = True
                    if alt_location != 'A':
                        continue

                if(check_atom_id is None and check_chain_id is None
                   and check_res_id is None):
                    check_atom_id = atom_id
                    check_res_id = res_id
                    check_chain_id = chain_id

                    new_chain = Chain(chain_id)
                    # register in structure
                    model.add(new_chain)

                    new_res = Residue(res_id, residue_name, '')
                    # register in chain
                    new_chain.add(new_res)
                    # name, coord, bfactor, occupancy, altloc, fullname,
                    # serial, number, element (None), charge(None)
                    new_atom = Atom(atom_name, atom_coord, info_1, info_2,
                                    ' ', full_atom_name, atom_id)
                    # set element
                    new_atom.set_element(element)
                    new_atom.set_charge(charge)
                    # register in res
                    new_res.add(new_atom)

                if chain_id != new_chain.id:
                    new_chain = Chain(chain_id)
                    # register in structure
                    model.add(new_chain)

                if res_id != new_res.id:
                    new_res = Residue(res_id, residue_name, '')
                    # register in chain
                    new_chain.add(new_res)

                if (atom_name.strip() != new_atom.id
                    or atom_id != new_atom.serial_number):
                    # name, coord, bfactor, occupancy, altloc, fullname,
                    # serial, number, element (None), charge(None)
                    new_atom = Atom(atom_name, atom_coord, info_1, info_2,
                                    alt_location, full_atom_name, atom_id)

                    # alt location change:
                    if alt_location != ' ':
                        new_res.flag_disordered()

                    # set element
                    new_atom.set_element(element)
                    new_atom.set_charge(charge)
                    # register in res
                    new_res.add(new_atom)

        # set structure
        self.structure = structure[0]

        if alt_location_status is True:
            print('Warning: This structure contains disorderd atoms! Using only A-Positions!')

    def save_to_file(self, path):
        '''
        Write the structure to a file.

        Note: We do not work with multiple models!

        :param path: path for the new pdb file.
        :type  path: str

        :raises AttributeError: if **path** does not end with '.pdb' or
           '.pdb.gz'

        .. note:

            If **path** points to a file, it will be overwritten without
            confirmation message.

        '''
        if path.endswith(".pdb"):
            with open(path, 'w') as f:
                self._save_to_file(f)
        elif path.endswith(".pdb.gz"):
            with gzip.open(path, "w") as f:
                self._save_to_file(f)
        else:
           raise AttributeError("Path does not end with '.pdb' or '.pdb.gz'\n{0}".format(path))

    def _save_to_file(self,f):
        for line in self.remarks:
            f.write(line)

        # iterate over chains
        for chain in self.structure:
            # iterate over residues
            for res in chain:
                # iterate over atoms
                for atom in res:
                    line = self._get_atom_line(chain, res, atom)
                    f.write(line)
            f.write('TER\n')

        if len(self.structure) == 1:
            chain = self.structure.child_list[0].id
            all_coords = self.get_all_atom_coords()
            if len(self.structure[chain]) == len(all_coords):
                # write connect
                for i, coord in enumerate(all_coords):
                    j = i + 1
                    if i == 0:
                        f.write('CONECT{0:5}{1:5}{2:5}{2:59}\n'.format(j,j+1,' ',''))
                    elif 0 < i < len(all_coords)-1:
                        f.write('CONECT{0:5}{1:5}{2:5}{2:59}\n'.format(j,j-1,j+1,''))
                    else:
                        f.write('CONECT{0:5}{1:5}{2:5}{2:59}\n'.format(j,j-1,' ',''))

        f.write('END\n')


    def _get_atom_line(self, chain, res, atom):
        '''
        Format a PDB line.

        :param chain: chain
        :type  chain: str
        :param res: residue
        :type  res: str
        :param atom: atom
        :type  atom: str

        :returns: Formatted PDB record.
        :rtype: str
        '''
        if res.get_id()[0] == ' ':
            record_type = "ATOM  "
        else:
            record_type = "HETATM"
        atom_number = atom.get_serial_number()
        name = atom.get_fullname()
        altloc = atom.get_altloc()
        resname = res.get_resname()
        chain_id = chain.get_id()
        resseq = res.get_id()[1]
        icode = res.get_id()[2]
        x, y, z = atom.get_coord()
        # occupancy
        info_1 = atom.get_info_1()
        # bfactor
        info_2 = atom.get_info_2()
        # atm there is no use for these parameters
        segid = ' '
        element = atom.get_element()
        charge = atom.get_charge()

        pdb_str = "{0}{1:5} {2:3}{3}{4:3} {5}{6:4}{7}   {8:8.3f}{9:8.3f}{10:8.3f}{11:6.2f}{12:6.2f}       {13:3}{14:3}{15}\n"
        #pdb_str = "{0}{1:5} {2:3}{3}{4:3} {5}{6:4}{7}   {8:8.3f}{9:8.3f}{10:8.3f}{11:6.2f}{12:6.2f}       {13:4}{14:2}{15}\n"


        return pdb_str.format(record_type, atom_number, name, altloc, resname,
                              chain_id, resseq, icode, x, y, z, info_1, info_2,
                              segid, element, charge)


class PQRFile(Structure_Template):
    '''
    Container for PQR files.
    '''

    def __init__(self, pqr_path = None):
        Structure_Template.__init__(self, pqr_path)
        # read structure
        if self.structure_path is None:
            self.structure = None
        else:
            self._read_structure()
        # i am a 'pqr'
        self.what_am_i = 'pqr'

    def _read_structure(self):
        '''
        This method reads the data from the given file.
        '''
        if self.structure_path.endswith(".pqr"):
            with open(self.structure_path) as f:
                content = f.readlines()
        elif self.structure_path.endswith(".pqr.gz"):
            with gzip.open(self.structure_path) as f:
                content = f.readlines()
        else:
            raise AttributeError("Path does not end with either '.pqr' or '.pqr.gz'\n{0}".format(self.structure_path))

        check_chain_id = None
        check_res_id = None
        check_atom_id = None

        # first 'model' -> 0
        structure = Structure(0)
        model = Model(0)
        structure.add(model)

        for line in content:
            if len(line) > 81 and line[0:6] != 'REMARK':
                raise AttributeError("Wrong formatted pqr files, line is longer than 81 characters!\n{0}".format(line))
            if (line[0:4] != 'ATOM' and line[0:6] != 'HETATM' and
                    line[0:3] != 'TER' and line[0:3] != 'END' and
                    line[0:5] != 'MODEL' and line[0:6] != 'ANISOU' and
                    line[0:6] != 'ENDMDL'):
                self.remarks.append(line)

            if line[0:4] == 'ATOM' or line[0:6] == 'HETATM':
                atom_id = int(line[6:11])
                atom_name = line[12:16].strip()
                full_atom_name = line[12:16]
                res_id = int(line[22:26])
                iCode = line[26]
                residue_name = line[17:20].strip()
                chain_id = line[21:22]
                x_coord = float(line[30:38])
                y_coord = float(line[38:46])
                z_coord = float(line[47:54])
                atom_coord = np.array([x_coord, y_coord, z_coord])
                # here info_1 and info_2 are necessary, because it is a pqr
                # file, so if they are missing, it is an error!
                info_1 = float(line[54:62])
                info_2 = float(line[62:70])


                if line[0:4] == 'ATOM':
                    het_flag = ' '
                elif line[0:6] == 'HETATM':
                    het_flag = 'H'

                # biopython artifact
                res_id = (het_flag, res_id, iCode)

                if(check_atom_id is None and check_chain_id is None
                   and check_res_id is None):
                    check_atom_id = atom_id
                    check_res_id = res_id
                    check_chain_id = chain_id

                    new_chain = Chain(chain_id)
                    # register in structure
                    model.add(new_chain)

                    new_res = Residue(res_id, residue_name, '')
                    # register in chain
                    new_chain.add(new_res)

                    new_atom = Atom(atom_name, atom_coord, info_1, info_2,
                                    iCode, full_atom_name, atom_id)
                    # register in res
                    new_res.add(new_atom)

                if chain_id != new_chain.id:
                    new_chain = Chain(chain_id)
                    # register in structure
                    model.add(new_chain)

                if res_id != new_res.id:
                    new_res = Residue(res_id, residue_name, '')
                    # register in chain
                    new_chain.add(new_res)

                if (atom_name.strip() != new_atom.id
                    or atom_id != new_atom.serial_number):
                    new_atom = Atom(atom_name, atom_coord, info_1, info_2,
                                    iCode, full_atom_name, atom_id)
                    # register in res
                    new_res.add(new_atom)

        # set structure
        self.structure = structure[0]

    def save_to_file(self, path):
        '''
        Write the structure to a file.

        Note: We do not work with multiple models!

        :param path: path for the new pqr file.
        :type  path: str

        :raises AttributeError: if **path** does not end with '.pqr' or
           '.pqr.gz'

        .. note:

            If **path** points to a file, it will be overwritten without
            confirmation message.

        '''
        if path.endswith(".pqr"):
            with open(path, 'w') as f:
                self._save_to_file(f)
        elif path.endswith(".pqr.gz"):
            with gzip.open(path, "w") as f:
                self._save_to_file(f)
        else:
           raise AttributeError("Path does not end with '.pqr' or '.pqr.gz'\n{0}".format(path))

    def _save_to_file(self,f):
        for line in self.remarks:
            f.write(line)
        # iterate over chains
        for chain in self.structure:
            # iterate over residues
            for res in chain:
                # iterate over atoms
                for atom in res:
                    line = self._get_atom_line(chain, res, atom)
                    f.write(line)
            f.write('TER\n')


    def _get_atom_line(self, chain, res, atom):
        '''
        Formats the line for writing.
        '''
        #if res.get_id()[0] == ' ':
        #    record_type = "ATOM  "
        #else:
        #    record_type = "HETATM"
        # everything is an ATOM
        record_type = "ATOM  "
        atom_number = atom.get_serial_number()
        name = atom.get_fullname()
        altloc = atom.get_altloc()
        resname = res.get_resname()
        chain_id = chain.get_id()
        resseq = res.get_id()[1]
        icode = res.get_id()[2]
        x, y, z = atom.get_coord()
        # charge
        info_1 = atom.get_info_1()
        # radius
        info_2 = atom.get_info_2()

        pqr_str = "{0}{1:5} {2:3}{3}{4:3} {5}{6:4}{7}   {8:8.3f}{9:8.3f}{10:8.3f}{11:8.4f}{12:7.4f}\n"

        return pqr_str.format(record_type, atom_number, name, altloc, resname,
                              chain_id, resseq, icode, x, y, z, info_1, info_2)


    def snap_crg_to_box(self, box_mesh_size, box_dim, box_offset, warning=True):
        '''
        This method snaps the charge for each atom of this pqr structure to a
        given dxbox. The new array contains the charge in the unit of Coulomb.

        If any coordinate lies outside the box an error will be printed to the
        standard output.

        Note: with uncharged ligands the dipole effect can be masked
        if the mesh size is too large.

        :param box_mesh_size: mesh size
        :type  box_mesh_size: :class:`np.array[3]`
        :param box_dim: box dimesions
        :type  box_dim: :class:`np.array[3]`
        :param box_offset: box origin
        :type  box_offset: :class:`np.array[3]`
        :param warning: print an error if the structure does not fit
           completely into the given box dimensions
        :type  warning: bool

        :returns: OpenDX box with Coulomb charges at the center of every atom.
        :rtype: :class:`np.array[:,:,:]`
        '''
        esp_box = np.zeros(box_dim)
        mesh_size = float(box_mesh_size[0])
        if not (box_mesh_size[0] == box_mesh_size[1] == box_mesh_size[2]):
            print ('WARNING: mesh size different depending on the direction, '
                   'using only the x-axis to snap vdw radii to grid space')
        
        # list which stores the coordinates of all atoms
        coordinates = []

        # list which contains charges
        esp_charge_list = []

        # we need a dxbox, because it lets us transform the coordinates
        dxbox = DXBox(esp_box, box_mesh_size, box_offset)

        # error warning
        error_warning = False

        #eC = 1.60217646e-19

        for atom in self.structure.get_atoms():
            box_coord = list(dxbox.transform_real_to_box_space(list(atom.get_coord())))
            # we do not want to leave the box!
            if(0 <= box_coord[0] < box_dim[0] and 0 <= box_coord[1] < box_dim[1]
               and 0 <= box_coord[2] < box_dim[2]):
                # check if coordinates already in the list
                if box_coord not in coordinates:
                    coordinates.append(box_coord)

                    esp_charge_list.append(atom.get_info_1())
                else:
                    position_index = coordinates.index(box_coord)
                    # add new charge to the old one
                    esp_charge_list[position_index] += atom.get_info_1()

            else:
                error_warning = True

        if warning is True and error_warning is True:
            print('{0} ({1}) has left the box!'.format(self.what_am_i,
                                                       self.structure_path))

        coordinates = np.array(coordinates)

        # set charge
        esp_box[(coordinates[:][:, 0], coordinates[:][:, 1], coordinates[:][:, 2])] = esp_charge_list

        return esp_box



class LatFile(Structure_Template):
    '''
    Container for Lattice files.
    '''

    def __init__(self, seq = None, fold = None, lattice = 'fcc',
            mesh_size = 3.8):
        Structure_Template.__init__(self, None)

        # make structure
        if seq is not None and fold is not None:
            self._make_structure(seq, fold, lattice, mesh_size)

        self.what_am_i = 'lat'

    def _make_structure(self, seq, fold, lattice, mesh_size):
        '''
        This method reads the data from the given file.
        '''
        if lattice.lower() == 'fcc':
            coord_list = self._convert_fcc_fold(fold, mesh_size)
        elif lattice.lower() == 'sc':
            coord_list = self._convert_sc_fold(fold, mesh_size)
        else:
            raise ValueError('Unkown lattice type: {0}'.format(lattice))

        check_chain_id = None
        check_res_id = None
        check_atom_id = None

        # first 'model' -> 0
        structure = Structure(0)
        model_counter = 0
        model = Model(model_counter)
        structure.add(model)
        counter = 0
        for i,(coord, char) in enumerate(zip(coord_list, seq)):
            atom_id = i+1
            atom_name = 'CA'
            full_atom_name = ' CA '
            alt_location = ' '
            res_id = i+1
            residue_name = self.one2threeletter[char]
            chain_id = 'A'
            x_coord, y_coord, z_coord = coord
            atom_coord = np.array([x_coord, y_coord, z_coord])
            info_1 = 1
            info_2 = 0
            element = ' C '
            charge = ' '
            het_flag = 'h'

            # biopython artifact
            res_id = (het_flag, res_id, ' ')

            if(check_atom_id is None and check_chain_id is None
               and check_res_id is None):
                check_atom_id = atom_id
                check_res_id = res_id
                check_chain_id = chain_id

                new_chain = Chain(chain_id)
                # register in structure
                model.add(new_chain)

                new_res = Residue(res_id, residue_name, '')
                # register in chain
                new_chain.add(new_res)
                # name, coord, bfactor, occupancy, altloc, fullname,
                # serial, number, element (None), charge(None)
                new_atom = Atom(atom_name, atom_coord, info_1, info_2,
                                ' ', full_atom_name, atom_id)
                # set element
                new_atom.set_element(element)
                new_atom.set_charge(charge)
                # register in res
                new_res.add(new_atom)

            if chain_id != new_chain.id:
                new_chain = Chain(chain_id)
                # register in structure
                model.add(new_chain)

            if res_id != new_res.id:
                new_res = Residue(res_id, residue_name, '')
                # register in chain
                new_chain.add(new_res)

            if (atom_name.strip() != new_atom.id
                or atom_id != new_atom.serial_number):
                # name, coord, bfactor, occupancy, altloc, fullname,
                # serial, number, element (None), charge(None)
                new_atom = Atom(atom_name, atom_coord, info_1, info_2,
                                alt_location, full_atom_name, atom_id)

                # alt location change:
                if alt_location != ' ':
                    new_res.flag_disordered()

                # set element
                new_atom.set_element(element)
                new_atom.set_charge(charge)
                # register in res
                new_res.add(new_atom)

        # set structure
        self.structure = structure[0]


    def save_to_file(self, path):
        '''
        Write the structure to a file. Try to calculate the
        correct connections. The Assumption of a lattice is made on the
        number of residues and the number of overall atom coordinates.
        If both numbers are equal, it is very probable that this is a lattice
        protein. This only works, if there is just one chain!

        Note: We do not work with multiple models!

        :param path: path for the new lattice file
        :type  path: str

        .. note:

            If **path** points to a file, it will be overwritten without
            confirmation message.

        '''
        with open(path, 'w') as f:
            for line in self.remarks:
                f.write(line)
            # iterate over chains
            for chain in self.structure:
                # iterate over residues
                for res in chain:
                    # iterate over atoms
                    for atom in res:
                        line = self._get_atom_line(chain, res, atom)
                        if len(line) > 0:
                            f.write(line)
                f.write('TER\n')

            if len(self.structure) == 1:
                chain = self.structure.child_list[0].id
                all_coords = self.get_all_atom_coords()
                if len(self.structure[chain]) == len(all_coords):
                    # write connect
                    for i, coord in enumerate(all_coords):
                        j = i + 1
                        if i == 0:
                            f.write('CONECT{0:5}{1:5}{2:5}{2:59}\n'.format(j,j+1,' ',''))
                        elif 0 < i < len(all_coords)-1:
                            f.write('CONECT{0:5}{1:5}{2:5}{2:59}\n'.format(j,j-1,j+1,''))
                        else:
                            f.write('CONECT{0:5}{1:5}{2:5}{2:59}\n'.format(j,j-1,' ',''))

            f.write('END\n')


    def _get_atom_line(self, chain, res, atom):
        '''
        Formats the line for writing.
        '''
        if res.get_id()[0] == ' ':
            record_type = "ATOM  "
        else:
            record_type = "HETATM"
        atom_number = atom.get_serial_number()
        name = atom.get_fullname()
        altloc = atom.get_altloc()
        resname = res.get_resname()
        chain_id = chain.get_id()
        resseq = res.get_id()[1]
        icode = res.get_id()[2]
        x, y, z = atom.get_coord()
        # occupancy
        info_1 = atom.get_info_1()
        # bfactor
        info_2 = atom.get_info_2()
        # atm there is no use for these parameters
        segid = ' '
        element = atom.get_element()
        charge = atom.get_charge()

        pdb_str = "{0}{1:5} {2:3}{3}{4:3} {5}{6:4}{7}   {8:8.3f}{9:8.3f}{10:8.3f}{11:6.2f}{12:6.2f}       {13:3}{14:3}{15}\n"
        #pdb_str = "{0}{1:5} {2:3}{3}{4:3} {5}{6:4}{7}   {8:8.3f}{9:8.3f}{10:8.3f}{11:6.2f}{12:6.2f}       {13:4}{14:2}{15}\n"


        return pdb_str.format(record_type, atom_number, name, altloc, resname,
                              chain_id, resseq, icode, x, y, z, info_1, info_2,
                              segid, element, charge)

    def _convert_fcc_fold(self, fold, mesh_size):
        coord_list = []
        seed = [0,0,0]
        coord_list.append(seed)
        fold_not_alphabet = {'B':'F', 'F':'B', 'D':'U', 'U':'D', 'L':'R', 'R':'L'}
        for i in range(0, len(fold), 2):
            # check move
            if fold[i] == fold[i+1] or fold[i] == fold_not_alphabet[fold[i+1]]:
                raise ValueError("False fold move in '{0}' at position {0},{1}!".format(fold, i, i+1))
            char_list = [fold[i],fold[i+1]]
            new_coord = np.zeros(3)
            for char in char_list:
                if char == 'U':
                    new_coord = new_coord + np.array([0,0,1])
                elif char == 'D':
                    new_coord = new_coord + np.array([0,0,-1])
                elif char == 'F':
                    new_coord = new_coord + np.array([1,0,0])
                elif char == 'B':
                    new_coord = new_coord + np.array([-1,0,0])
                elif char == 'R':
                    new_coord = new_coord + np.array([0,-1,0])
                elif char == 'L':
                    new_coord = new_coord + np.array([0,1,0])
                else:
                    raise ValueError('Unkown fold instruction: {0}!'.format(char))

            # 1/sqrt(2) -> ffc
            new_coord = list(np.array(coord_list[-1])
                             + mesh_size / np.sqrt(2) * new_coord)
            if new_coord not in coord_list:
                coord_list.append(new_coord)
            else:
                raise ValueError("Collision detected in '{0}', at position {1}!".format(fold,i))

        return coord_list

    def _convert_sc_fold(self, fold, mesh_size):
        coord_list = []
        seed = [0,0,0]
        coord_list.append(seed)
        for i,char in enumerate(fold):
            if char == 'U':
                new_coord = np.array([0,0,1])
            elif char == 'D':
                new_coord = np.array([0,0,-1])
            elif char == 'F':
                new_coord = np.array([1,0,0])
            elif char == 'B':
                new_coord = np.array([-1,0,0])
            elif char == 'R':
                new_coord = np.array([0,-1,0])
            elif char == 'L':
                new_coord = np.array([0,1,0])
            else:
                raise ValueError('Unkown fold instruction: {0}!'.format(char))

            new_coord = list(np.array(coord_list[-1]) + mesh_size * new_coord)

            if new_coord not in coord_list:
                coord_list.append(new_coord)
            else:
                raise ValueError("Collision detected in '{0}', at position {1}!".format(fold,i))

        return coord_list


class Atom(BioAtom):
    '''
    Atom object.

    The Atom object stores atom name (both with and without spaces),
    coordinates, B factor, occupancy, alternative location specifier
    and (optionally) anisotropic B factor and standard deviations of
    B factor and positions.

    :param name: atom name (ex: 'CA'). Note that spaces are normally stripped
    :type  name: str
    :param coord: atomic coordinates (Angstroms)
    :type  coord: :class:`np.array[3]`
    :param bfactor: isotropic B factor
    :type  bfactor: float
    :param occupancy: occupancy (0.0-1.0)
    :type  occupancy: float
    :param altloc: alternative location specifier for disordered atoms
    :type  altloc: str
    :param fullname: full atom name, including spaces, e.g. ``' CA '``.
       Normally these spaces are stripped from the atom name.
    :type  fullname: str
    :param element: atom element, ex: ``'C'`` for carbon, ``'HG'`` for mercury
    :type  fullname: str
    '''

    def __init__(self, name, coord, occupancy, bfactor, altloc, fullname,
                 serial_number, element = None, charge = None):
        BioAtom.__init__(self, name, coord, bfactor, occupancy, altloc,
                         fullname,
                         serial_number,
                         element)
        self.level = "A"
        # Reference to the residue
        self.parent = None
        # the atomic data
        self.name = name      # eg. CA, spaces are removed from atom name
        self.fullname = fullname  # e.g. " CA ", spaces included
        self.coord = coord
        self.bfactor = bfactor
        self.occupancy = occupancy
        self.altloc = altloc
        self.full_id = None   # (structure id, model id, chain id, residue id, atom id)
        self.id = name        # id of atom is the atom name (e.g. "CA")
        # my addition
        if self.altloc != ' ':
            self.disordered_flag = 1
        else:
            self.disordered_flag = 0
        self.anisou_array = None
        self.siguij_array = None
        self.sigatm_array = None
        self.serial_number = serial_number
        # Dictionary that keeps addictional properties
        self.xtra = {}
        self.element = element
        self.charge = charge

    def _assign_element(self, element): 
        '''
        Not used!.
        '''
        return None

    def set_info_1(self, info_1):
        '''
        :setter: Set slot information 1 (pdb: occupancy, pqr: charge)
        :type: float
        '''
        self.set_occupancy(info_1)

    def set_info_2(self, info_2):
        '''
        :setter: Set slot information 2 (pdb: B-factor, pqr: atomic radius)
        :type: float
        '''
        self.set_bfactor(info_2)

    def get_info_1(self):
        '''
        :getter: Get slot information 1 (pdb: occupancy, pqr: charge)
        :type: float
        '''
        return self.get_occupancy()

    def get_info_2(self):
        '''
        :getter: Get slot information 2 (pdb: B-factor, pqr: atomic radius)
        :type: float
        '''
        return self.get_bfactor()

    def set_element(self, element):
        '''
        :setter: Set atom element
        :type: str
        '''
        self.element = element
        
    def get_element(self):
        '''
        :getter: Get atom element
        :type: str
        '''
        return self.element
        
    def set_charge(self, charge):
        '''
        :setter: Set atom charge
        :type: str
        '''
        self.charge = charge

    def get_charge(self):
        '''
        :getter: Get atom charge
        :type: str
        '''
        return self.charge
