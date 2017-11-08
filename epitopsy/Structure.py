'''
Created on Jan 11, 2012

@author: chris


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
    This is a template class from which PDBFile and PQRFile will be derived,
    because they share a lot of functions.

    This class is not supposed to work with multiple models (for example NMR
    structures)!

    List of class one has to implement for each derivative:
        * read structure
        * write structure
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

    def __init__(self, structure_path = None):
        '''
        Constructor
        '''
        if structure_path is not None:
            self.structure_path = structure_path
        else:
            self.structure_path = structure_path

        self.structure = None
        # 'pdb', 'pqr', ?
        self.what_am_i = None
        self.remarks = []

    def get_all_atom_coords(self):
        '''
        Returns:
            A List of atom coordinates.
        '''
        atoms_coord = []
        for atom in self.structure.get_atoms():
            atoms_coord.append(atom.get_coord())
        return atoms_coord

    def get_coords_from_atom_list(self, atom_list):
        '''
        Args:
            atom_list -> List of atom objects.

        Returns:
            A list of atom coordinates for the specified atom objects.
        '''
        atom_coord = []
        for atom in atom_list:
            atom_coord.append(atom.get_coord())
        return atom_coord

    def get_hetero_atoms(self):
        '''
        Returns:
            A list of all hetero atoms.
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
        Returns:
            A list of all hetero atom coordinates
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
        Returns:
            A list of all amino acid atoms.
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
        Returns:
            A list of all amino acid atom coordinates.
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
        Returns:
            A list of all atom coordinates
        '''
        all_atoms = []
        for atom in self.structure.get_atoms():
            all_atoms.append(atom)
        return all_atoms

    def get_info_1(self, atoms=None):
        '''
        Returns:
            A list of all atom information 1 values. For a pdb this is
            the occupancy and for a pqr it is the charge.
        '''
        if atoms is None:
            atoms = self.get_all_atoms()

        info_1_list = []
        for atom in atoms:
            info_1_list.append(atom.get_info_1())
        return info_1_list

    def get_info_2(self, atoms=None):
        '''
        Returns:
            A list of all atom information 2 values. For a pdb this is
            the temperature factor and for a pqr it is the radius.
        '''
        if atoms is None:
            atoms = self.get_all_atoms()

        info_2_list = []
        for atom in atoms:
            info_2_list.append(atom.get_info_2())
        return info_2_list

    def get_chain_ids(self):
        '''
        Returns:
            A list of all chain ids in this structure.
        '''
        chain_list = []
        for chain in self.structure:
            chain_list.append(chain.id)

        return chain_list

    def get_first_res_id(self):
        '''
        Returns:
            The integer number of the first amino acid.
        '''
        chain_id = self.structure.child_dict.keys()[0]
        return self.structure[chain_id].child_list[0].id[1]

    def get_atoms_of_type(self, atom_type):
        '''
        Args:
            atom_type -> Type of atom to return (e.g. 'CA')
        Returns:
            All atom objects of a certain type.
        '''
        all_atoms_of_type = []
        for atom in self.structure.get_atoms():
            if atom.get_name().upper() == atom_type.upper():
                all_atoms_of_type.append(atom)
        return all_atoms_of_type


    def transform(self, T):
        '''
        Transform the pdb structure with the given matrix.

        Args:
            T -> [3,3] numpy matrix  to transform the coordinates by
                matrix multiplication.

        Returns:
            None.
        '''
        for atom in self.structure.get_atoms():
            atomcoord = atom.get_coord()
            newcoord = np.dot(T, atomcoord)
            atom.set_coord(newcoord)

    def translate(self, transVector):
        '''
        Method to translate protein structure by the given vector.

        Args:
            transvector -> Numpy Array holding translation distances for each
                dimension.

        Returns:
            None.
        '''
        for atom in self.structure.get_atoms():
            atomcoord = atom.get_coord()
            newcoord = atomcoord + np.array(transVector)
            atom.set_coord(newcoord)

    def translate_x(self, dist):
        '''
        Method to translate protein structure in x direction.

        Args:
            dist -> Amount of displacement in x direction.

        Returns:
            None.
        '''
        self.translate([dist, 0, 0])

    def translate_y(self, dist):
        '''
        Method to translate protein structure in y direction.

        Args:
            dist -> Amount of displacement in y direction.

        Returns:
            None.
        '''
        self.translate([0, dist, 0])

    def translate_z(self, dist):
        '''
        Method to translate protein structure in z direction.

        Args:
            dist -> Amount of displacement in z direction.

        Returns:
            None.
        '''
        self.translate([0, 0, dist])

    def translate_origin_and_rotate(self, phi, theta, psi):
        '''
        This methods centers the structure at the origin, rotates it with
        angle_x around the x axis (angle_y around y axis, etc.) and moves
        it back to where it was.

        Args:
            phi -> Euler angle for rotation.
            theta -> Euler angle for rotation.
            psi -> Euler angle for rotation.

        Returns:
            None.
        '''
        transvector = self.determine_geometric_center()
        self.translate(-transvector)
        self.pRotate(phi, theta, psi)
        self.translate(transvector)

    def move_to_new_position(self, new_coord):
        '''
        This method moves the geometric center of the structure to the
        supplied coordinates.

        Args:
            new_coord -> list/numpy array of the new coordinates

        Returns:
            None.
        '''
        old_coord = self.determine_geometric_center()
        self.translate(np.array(new_coord) - np.array(old_coord))

    def rotate_and_move_to_new_position(self, phi, theta, psi, new_coord):
        '''
        This method centers the structure at (0,0,0), rotates it and the
        moves it to a new position.

        Args:
            phi -> Euler angle for rotation.
            theta -> Euler angle for rotation.
            psi -> Euler angle for rotation.
            new_coord -> New coordination for the center of geometry.

        Returns:
            None.
        '''
        self.translate(-self.determineGeometricCenter())
        self.pRotate(phi, theta, psi)
        self.translate(np.array(new_coord))

    def rotate(self, angle, axis):
        '''
        Method to rotate protein structure.
        For the rotation I use the Rodrigues' rotation formula.

        Args:
            degree -> Angle by which to rotate
            axis -> Axis around which to rotate

        Returns:
            None.
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
        Args:
            degree -> Angle for rotation around the x axis.

        Returns:
            None.
        '''
        self.rotate(degree, [1, 0, 0])

    def rotateY(self, degree):
        '''
        Args:
            degree -> Angle for rotation around the y axis.

        Returns:
            None.
        '''
        self.rotate(degree, [0, 1, 0])

    def rotateZ(self, degree):
        '''
        Args:
            degree -> Angle for rotation around the z axis.

        Returns:
            None.
        '''
        self.rotate(degree, [0, 0, 1])

    def pRotate(self, phi, theta, psi):
        '''
        Apply euler angle rotation to the structure.

        Args:
            phi -> Euler angle for rotation.
            theta -> Euler angle for rotation.
            psi -> Euler angle for rotation.

        Returns:
            None.
        '''
        euler_rotation = MathTools.calculate_rotation_matrix(phi, theta, psi)
        self.transform(euler_rotation)

    def rotate_by_matrix(self, rot_matrix):
        '''
        Args:
            rot_matrix -> [3,3] numpy matrix to transform the coordinates by
                matrix multiplication.

        Returns:
            None.
        '''
        self.transform(rot_matrix)

    def determineCenterOfMass(self):
        '''
        Method to determine the center of mass for the protein structure.

        returns:
            None.
        '''
        print('not implemented ... missing atom weights!')

    def determine_geometric_center(self):
        '''
        Returns:
            A vector pointing to the geometric center of the structure.
        '''
        all_coords = np.array(self.get_all_atom_coords())
        return np.mean(all_coords, 0)

    def determine_center_of_extremes_of_atoms(self, atoms):
        '''
        Args:
            atoms -> List of atom objects.

        Returns:
            A vector pointing to the geometric center of the coordination
            extremes of the given atom coordinates.
        '''
        centerVector = np.zeros(3)
        extremes = self.determine_coordinate_extremes(atoms)
        centerVector[0] = (extremes[0][1] + extremes[0][0]) / 2.0
        centerVector[1] = (extremes[1][1] + extremes[1][0]) / 2.0
        centerVector[2] = (extremes[2][1] + extremes[2][0]) / 2.0
        return centerVector

    def determine_center_of_extremes(self):
        '''
        Returns:
            A vector pointing to the geometric center of the coordination
            extremes of this structure.
        '''
        centerVector = np.zeros(3)
        extremes = self.determine_coordinate_extremes()
        centerVector[0] = (extremes[0][1] + extremes[0][0]) / 2.0
        centerVector[1] = (extremes[1][1] + extremes[1][0]) / 2.0
        centerVector[2] = (extremes[2][1] + extremes[2][0]) / 2.0
        return centerVector

    def determine_max_diameter(self, atoms = None):
        '''
        Args:
            atoms -> Optional a list of atom objects, otherwise it uses all
                atoms of this structure and calculates the maximum diameter.

        Returns:
            A float number of the maximum diameter.
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
        Determine the geometric center and calculate the minimal radius that
        encapsulates all atoms.

        Args:
            atoms -> Optional a list of atom objects, otherwise it uses all
                atoms of this structure and calculates the radius.

        Returns:
            A float number of the radius.
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
        Translate geometric center to (0.0/0.0/0.0).

        Returns:
            None.
        '''
        center_vector = self.determine_geometric_center()
        self.translate(-center_vector)

    def determine_coordinate_extremes(self, atoms = None):
        '''
        Args:
            atoms -> Optional a list of atom objects, otherwise it uses all
                atoms of this structure and calculates the extreme coordinates
                in each  direction.
        Returns:
            Extreme values in each direction as a 3*2 array.
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
        :returntype: :class:`numpy.ndarray`[3,3]
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
        :type  base: :class:`numpy.ndarray`[3,3]
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
        This method calculates the radius of gyration. It is the maximum
        distance of an atom to the geometrical center.

        Returns:
            A float number as the radius of gyration.
        '''
        all_coords = np.array(self.get_all_atom_coords())
        geo_center = np.mean(all_coords,0)
        r_g = 0

        diff = all_coords - geo_center
        r = np.sqrt(diff[:,0]**2 + diff[:,1]**2 + diff[:,2]**2)

        r_g = np.amax(r)

        return r_g

    def clone(self, chain_id_list = None,
              res_id_list = None, res_type_list = None, atom_types_list = None):
        '''
        This method returns a clone of this structure. Through the list
        parameters specific items can be selected. If one supplies one letter
        codes for the residues they will be translated to three letter codes!

        Args:
            chain_id_list -> List of chains to copy to the clone.
            res_id_list -> List of residues in each chain to copy to the
                clone.
            res_type_list -> Types of residues to copy to the new clone.
            atom_types_list -> Types of atoms to copy to the new clone.

        Returns:
            A new PDBFile / PQRFile / LatFile object.
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

    def get_residue_id_list(self, chain_id = None):
        '''
        Args:
            chain_id -> Optionally, if None, it uses the all available chains.
        Returns:
            A list with all residue ID's of the structure.
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
        Returns:
            A dictionary of the first chain, which contains the residue id as
            key and the corresponding amino acid as value.
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
            print('Encountered the following non amino acids: {0}'.format(non_amino_acid_list))

        return res_map

    def contains_chain_break(self, chain_id = None):
        '''
        Args:
            chain_id -> Optionally, if None, it uses the all available chains.

        Returns:
            Either True (chain break) or False (no chain break).
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
                                    return True
                        res_id_list.append(res_id)
        else:
            res_id_list = []
            for i,res in enumerate(self.structure[chain_id]):
                if res.id[0] != "H":
                    res_id = res.get_id()[1]
                    if i > 0:
                        if res_id_list[-1] + 1 != res_id:
                                return True
                    res_id_list.append(res_id)

        ## no chain break
        return False


    def get_res_id_array_mapping(self):
        '''
        Returns:
            A dictionary with residue ids as keys and values, which can be
            used as an index in an array.
        '''
        id_dict = {}
        chain = self.structure.child_list[0].id
        for i,res in enumerate(self.structure[chain]):
            id_dict[res.id[1]] = i

        return id_dict

    def get_residue_names_from_res_id_list(self, res_id_list, chain_id = None):
        '''
        Args:
            res_id_list -> List of residue ids.
            chain_id -> Optionally, if None, it uses the all available chains.

        Returns:
            A list with the names of the residue ids in the
            given list.
        '''
        if chain_id is None:
            # if no chain id has been supplied take the first one
            chain_id = self.structure.child_dict.keys()[0]

        res_name = []
        for res in self.structure[chain_id].child_list:
            if res.get_id()[1] in res_id_list:
                res_name.append(res.resname)

        return res_name

    def get_residues(self, chain_id = None, res_id_list = None):
        '''
        This method returns a list with residue objects of the residue ids in
        the given list, if None is given, it returns all residues.

        Args:
            chain_id -> Optionally, if None, it uses the all available chains.
            res_id_list -> List of residue ids from which one wants the
                objects.

        Returns:
            A list of residues objects matching the given criteria.
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
                res_list.append(res.resname)

        return res_list

    def get_atoms_by_id(self, atom_id_list):
        '''
        This function returns the Atoms,  by their corresponding number from
        the pdb-file.

        Args:
            atom_id_list -> List of atom id numbers.

        Returns:
            A list of atom objects to the matching atom ids.
        '''
        atom_list = []
        for atom in self.get_all_atoms():
            if atom.serial_number in atom_id_list:
                atom_list.append(atom)
        return atom_list

    def get_atoms_close_to_reference(self, reference_point, max_radius,
            min_radius = 0, atoms = None):
        '''
        This function returns all atoms of this structure, which lie in the
        range [min_radius, max_radius] to the reference point.

        Args:
            reference_point -> Numpy array of the reference.
            max_radius -> Maximal distance to include the atoms.
            min_radius -> Minimal distance of the atoms to the the reference.
            atoms -> Optional list of atoms, if None is given it uses all
                atoms from the protein.

        Returns:
            A list of atom objects close to the given reference point.
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
        Finds Atoms of chain1 which are within max_distance of chain2

        Args:
            chain1 -> First chain id.
            chain2 -> Second chain id.
            max_distance -> Maximal distance to include atoms in the
                calculation.

        Returns:
            A list of atom objects from chain1 which distance to chain2 is
            smaller than max_distance.
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
        Returns:
            A dictionary which contains the chain and the sequence:
            'A' : 'RG...CC'
            'B' : 'PW...FV'
            Non standard amino acids will not be returned!!!
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
                        warning = True, vdw_radii = 2,
                        increase_vdw_by = 0):
        '''
        This method snaps a structure to a given dxbox. If the structure is
        a pqr it uses the supplied vdw radii otherwise it uses the variable
        'vdw_radii' for each atom.
        
        If any coordinate lies outside the box an error will be printed to the
        standard output.
        
        :param box_mesh_size: mesh size
        :type  box_mesh_size: tuple(float,float,float)
        :param box_dim: box dimesions
        :type  box_dim: tuple(int,int,int)
        :param box_offset: box origin
        :type  box_offset: tuple(float,float,float)
        :param warning: print an error if the structure does not fit
           completely into the given box dimensions
        :type  warning: bool
        :param vdw_radii: default vdw radius if not a PQR file (Angstroms)
        :type  vdw_radii: float
        :param increase_vdw_by: increase all vdw radii by a constant (Angstroms)
        :type  increase_vdw_by: float
        :returns: OpenDX box with 0's outside the molecule and 1's inside
        :returntype: :class:`np.ndarray[:,:,:]`
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
                                                         atom_types = None):
        '''
        This method tries to fit the given pdb onto this one. The method
        returns a dictionary, which contains the 'rmsd', 'rotation' matrix
        and the 'translation'.
        The parameter 'atom_types' can be used to supply a list of atoms, which
        should be fit onto each other, if 'None' is supplied, it will try to
        fit all atoms onto each other.
        In case the number of atoms does not match, it raises an Error!

        I guess the units are Angstroem.

        Args:
            pdb_to_rotate -> Structure derivative object to superimpose onto
                this object.
            atom_types -> If None, it tries to fit all atoms, but it can also
                be used to fit specific types (e.g. ['CA'], ['CA','N'])

        Returns:
            A dictionary which contains the following keys:
            * 'rmsd' : root mean square deviation
            * 'rotation' : rotation matrix
            * 'translation' : translation vector
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
                                        atom_types = None):
        '''
        This method superimposes this structure onto the given structure!
        If the number of atoms differ it raises an error!

        Args:
            pdb_to_superimpose ->  Structure derivative object to superimpose onto
                this object.
            atom_types -> If None, it tries to fit all atoms, but it can also
                be used to fit specific types (e.g. ['CA'], ['CA','N'])

        Returns:
            None.
        '''
        super_instructions = self.get_rmsd_rotation_translation_from_superposition(pdb_to_superimpose, atom_types)
        trans_vector = super_instructions['translation']
        rot_matrix = super_instructions['rotation']

        for atom in pdb_to_superimpose.structure.get_atoms():
            atom.transform(rot_matrix, trans_vector)

    def superimpose_self_onto_given_pdb(self, pdb_to_superimpose,
                                        atom_types = None):
        '''
        This method superimposes this structure onto the given structure!
        If the number of atoms differ it raises an error!

        Args:
            pdb_to_superimpose ->  Structure derivative object to superimpose
                this object onto.
            atom_types -> If None, it tries to fit all atoms, but it can also
                be used to fit specific types (e.g. ['CA'], ['CA','N'])

        Returns:
            None.
        '''
        super_instructions = pdb_to_superimpose.get_rmsd_rotation_translation_from_superposition(self, atom_types)
        trans_vector = super_instructions['translation']
        rot_matrix = super_instructions['rotation']

        for atom in self.structure.get_atoms():
            atom.transform(rot_matrix, trans_vector)

    def get_dxbox_dim(self, box_mesh_size, extend=None, cubic_box=False,
                      nlev=4):
        '''
        Calculate the dimensions of an APBS box large enough to contain the
        protein. The box center is taken as the protein geometric center.
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

        Args:
            box_mesh_size -> [m,m,m]
            box_dim -> [x,y,z]
            box_center -> [x_c,y_c,z_c]

        Returns:
            A list box_offset: [x_o,y_o,z_o].
        '''
        box_offset = np.array(box_center) - (np.array(box_mesh_size) *
                                             np.array(box_dim) // 2)
        return box_offset

    def get_hydrophobic_potential(self, box_mesh_size, box_dim, box_offset):
        '''Calculate the hydrophobic potential.
        This method uses a simplified model for the hydrophobic potential.
        The 'charges' are taken from the Kyte and Doolittle hydrophobicity
        scale. For each residue the center is calculated and the potential
        is modelled as:
            phi = sum ( hydrophobic_charge * exp( - distance ) )

        Args:
            box_mesh_size: meshsize of the grid
            box_dim: dimension of the grid
            box_offset: offset of the grid

        Returns:
            Numpy array with the hydrophoic potential for each grid
            point.
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
            increase_vdw_by = 0):
        '''Get the van der Waals hull of the protein.

        Args:
            box_mesh_size: meshsize of the grid
            box_dim: dimension of the grid
            box_offset: offset of the grid
            vdw_radii: atom radii, if this is not a pqr
            increase_vdw_by: extend each radius by this value

        Returns:
            Numpy array with 1's at the hull grid points and 0's everywhere
            else.
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
            energy_cutoff = 1., vdw_radii = 2.):
        '''Calculate the number of atoms in this structure object that are
        overlapping with the given structure object.

        Args:
            other_structure_object: Structure which might have an overlapp with
                this one
            energy_cutoff: cutoff for the lennard jones potential to decide
                if there is an overlapp or not
            vdw_radii: If this is pdb object there are no radii information
                available

        Returns:
            Integer: counts of atoms in this structure that overlapp with
            the given structure object.
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

    def get_contact_list(self, cutoff = 5.):
        '''
        Args:
            cutoff -> if any distance between atoms from two residues is less
                than the cutoff, it is a contact

        Returns:
            A list:
            [[...], [i,j,0], [...]] in this case 0 means no contact.
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

    def __init__(self, pdb_path = None):
        Structure_Template.__init__(self, pdb_path)
        # read structure
        if self.structure_path is None:
            self.structure = None
        else:
            self._read_file()
        # i am a 'pdb'
        self.what_am_i = 'pdb'

    def get_pqr_structure(self, new_pqr_path = None, force_field = "amber",
            pdb2pqr_argv = None, pdb2pqr_path = "pdb2pqr",
            add_ions = True, add_chain = True):
        '''Call pdb2pqr to replace the b-factor and the occupancy information
        in the pdb file with the charges and the vdw radii. If this pdb Object
        has no structure_path property (i.e. it is None), then an error is
        raised.

        If the pdb contains CA,ZN or SO4 ions, they will be added, if not
        stated otherwise.

        Args:
            new_pqr_path -> Path for the new pqr file. If None is given, it
            replaces the ".pdb" with ".pqr" at the end.
            force_field -> Forcefield from which charges and radii should be
                taken. Default is amber.
            pdb2pqr_argv -> Can contain additional arguments to pdb2pqr as a
                list (e.g. ['--assign-only'], oder ['--noopt']). If multiple
                additional arguments are given, they also have to be given as
                a list (e.g. ['--assign-only', '--noopt']).
            pdb2pqr_path -> Path to the executable of PDB2PQR.
            add_ions -> Add CA, ZN, SO4 ions if they are in the pdb. The CA,ZN
                and SO4 atoms should have a residue name that fits their type
                (CA, ZN, SO4).

        Returns:
            An PQRFile object.
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
        This method writes the structure to a given path.
        If the structure is a lattice, it will try to calculate the
        correct connections. The Assumption of a lattice is made on the
        number of residues and the number of overall atom coordinates.
        If both numbers are equal, it is very probable that this is a lattice
        protein. This only works, if there is just one chain!

        Notice: We do not work with multiple models!
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
        icode = ' '
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

                    new_atom = Atom(atom_name, atom_coord, info_1, info_2,
                                    ' ', full_atom_name, atom_id)
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
                                    ' ', full_atom_name, atom_id)
                    # register in res
                    new_res.add(new_atom)

        # set structure
        self.structure = structure[0]

    def save_to_file(self, path):
        '''
        This method writes the structure to a given path.

        Notice: We do not work with multiple models!
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
        icode = ' '
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
        
        Be carefull!!! If you use a neutral probe and you mesh size is to
        large the dipole effect is not visible!
        
        :param box_mesh_size: mesh size
        :type  box_mesh_size: tuple(float,float,float)
        :param box_dim: box dimesions
        :type  box_dim: tuple(int,int,int)
        :param box_offset: box origin
        :type  box_offset: tuple(float,float,float)
        :param warning: print an error if the structure does not fit
           completely into the given box dimensions
        :type  warning: bool
        :returns: OpenDX box with Coulomb charges at the center of every atom
        :returntype: :class:`np.ndarray[:,:,:]`
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
        This method writes the structure to a given path.
        If the structure is a lattice, it will try to calculate the
        correct connections. The Assumption of a lattice is made on the
        number of residues and the number of overall atom coordinates.
        If both numbers are equal, it is very probable that this is a lattice
        protein. This only works, if there is just one chain!

        Notice: We do not work with multiple models!
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
        icode = ' '
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
    def __init__(self, name, coord, occupancy, bfactor, altloc, fullname,
                 serial_number, element = None, charge = None):
        """
        Atom object.

        The Atom object stores atom name (both with and without spaces),
        coordinates, B factor, occupancy, alternative location specifier
        and (optionally) anisotropic B factor and standard deviations of
        B factor and positions.

        @param name: atom name (eg. "CA"). Note that spaces are normally stripped.
        @type name: string

        @param coord: atomic coordinates (x,y,z)
        @type coord: Numeric array (Float0, size 3)

        @param bfactor: isotropic B factor
        @type bfactor: number

        @param occupancy: occupancy (0.0-1.0)
        @type occupancy: number

        @param altloc: alternative location specifier for disordered atoms
        @type altloc: string

        @param fullname: full atom name, including spaces, e.g. " CA ". Normally
        these spaces are stripped from the atom name.
        @type fullname: string

        @param element: atom element, e.g. "C" for Carbon, "HG" for mercury,
        @type fullname: uppercase string (or None if unknown)

        """
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
        This is just a more abstract method which I will use instead of
        'get_bfactor' and 'get_occupancy'.
        pdb: occupancy
        pqr: charge
        '''
        self.set_occupancy(info_1)

    def set_info_2(self, info_2):
        '''
        This is just a more abstract method which I will use instead of
        'get_bfactor' and 'get_occupancy'.
        pdb: b_factor
        pqr: vdw radius
        '''
        self.set_bfactor(info_2)

    def get_info_1(self):
        '''
        This is just a more abstract method which I will use instead of
        'get_bfactor' and 'get_occupancy'.
        pdb: occupancy
        pqr: charge
        '''
        return self.get_occupancy()

    def get_info_2(self):
        '''
        This is just a more abstract method which I will use instead of
        'get_bfactor' and 'get_occupancy'.
        pdb: bfactor
        pqr: radius
        '''
        return self.get_bfactor()

    def set_element(self, element):
        '''
        Set element.
        '''
        self.element = element
        
    def get_element(self):
        '''
        Return element.
        '''
        return self.element
        
    def set_charge(self, charge):
        self.charge = charge

    def get_charge(self):
        return self.charge
