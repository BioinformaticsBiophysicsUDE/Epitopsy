
:mod:`Structure` --- Additional classes for PDB files
=====================================================

.. module:: Structure
   :synopsis: Additional classes for PDB files.
.. moduleauthor:: Christoph Wilms <christoph.wilms@uni-due.de>
.. moduleauthor:: Thomas Hamelryck <thamelry@binf.ku.dk>
.. sectionauthor:: Jean-Noel Grad <jean-noel.grad@uni-due.de>


This module provides additional classes for the handling of PDB files.
Following classes are part or derivatives of the original biopython software
(`homepage <http://biopython.org/wiki/Biopython>`_):
:class:`entity`, :class:`Structure`, :class:`Model`,
:class:`Chain`, :class:`Residue` and :class:`atom`.
They are therefore covered by their original Biopython license
(:ref:`read <BioLicense>`, or
:download:`download <../_static/licenses/biopython.txt>`).

.. _Structure-syntax:

Module Syntax
-------------

PDB files should be loaded this way:

    >>> a = Structure.PDBFile("/path/4N6W.pdb")
    Warning: This structure contains disorderd atoms! Using only A-Positions!
    >>> # further operations...

.. _contents-of-module-Structure:

Module Contents
---------------

Main file class
^^^^^^^^^^^^^^^

.. class:: Structure_Template(object)

    This is a template class from which PDBFile and PQRFile will be derived,
    because they share a lot of functions.

    This class is not supposed to work with multiple models (for example NMR
    structures)!

    List of class one has to implement for each derivative:

    * read structure
    * write structure

    .. method:: get_all_atom_coords()

        :returns: a List of atom coordinates

        Example::

            >>> a = Structure.PDBFile("/path/4N6W.pdb")
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


    .. method:: get_coords_from_atom_list(atom_list)

        :param atom_list: list of atom objects

        :returns: a list of atom coordinates for the specified atom objects

        Example::

            >>> a = Structure.PDBFile("/path/4N6W.pdb")
            >>> b = a.get_hetero_atoms()[1:3] # retrieving two atom objects
            >>> print b
            [<Atom FE>, <Atom CAC>]
            >>> print a.get_coords_from_atom_list(b) # getting coordinates
            [array([-28.139,  20.588,  11.849]), array([-30.158,  21.541,   5.145])]
            >>> print a.get_hetero_atoms_coords()[1:3] # just checking
            [array([-28.139,  20.588,  11.849]), array([-30.158,  21.541,   5.145])]

    .. method:: get_hetero_atoms()

    Navigate through all HETATM entries of the PDB file and return :class:`atom` objects.

        :returns: a list of all hetero atoms (:class:`atom`)

        Example::

            >>> a = Structure.PDBFile("/path/4N6W.pdb")
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

    .. method:: get_hetero_atoms_coords()

        :returns: a list of all hetero atom coordinates

        Example::

            >>> a = Structure.PDBFile("/path/4N6W.pdb")
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

    .. method:: get_amino_atoms()

        :returns: a list of all amino acid atoms (:class:`atom`)

        Example::

            >>> a = Structure.PDBFile("/path/4N6W.pdb")
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

    .. method:: get_amino_atoms_coords()

        :returns: a list of all amino acid atom coordinates.

        Example::

            >>> a = Structure.PDBFile("/path/4N6W.pdb")
            >>> for ATOM in a.get_amino_atoms_coords():
            ...     print ATOM
            ... 
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

    .. method:: get_all_atoms()

        :returns: a list of all atom coordinates

        Example::

            >>> a = Structure.PDBFile("/path/4N6W.pdb")
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

    .. method:: get_info_1(atoms=None)

        :returns: a list of all :class:`atom` information 1 values. For a pdb this is
            the occupancy and for a pqr it is the charge

    .. method:: get_info_2(atoms=None)

        :returns: a list of all :class:`atom` information 2 values. For a pdb this is
            the temperature factor and for a pqr it is the radius

    .. method:: get_chain_ids()

        :returns: a list of all chain ids in this structure

        Example::

            >>> a = Structure.PDBFile("/path/4N6W.pdb")
            >>> a.get_chain_ids()
            ['A']

    .. method:: get_first_res_id()

        :returns: the integer number of the first amino acid

        Example::

            >>> a = Structure.PDBFile("/path/4N6W.pdb")
            >>> a.get_first_res_id()
            1

    .. method:: get_atoms_of_type(atom_type)

        :param atom_type: type of atom to return (e.g. 'CA')
        :type atom_type: str

        :returns: a list of all :class:`atom` objects of a certain type

        Example::

            >>> a = Structure.PDBFile("/path/4N6W.pdb")
            >>> a.get_atoms_of_type("FE")
            [<Atom FE>, <Atom FE>]

    .. method:: transform(T)

        Transform the pdb structure with the given matrix.

        :param T: [3,3] numpy matrix to transform the coordinates by
            matrix multiplication

        :returns: ``None``

    .. method:: translate(transVector)

        Method to translate protein structure by the given vector.

        :param transvector: Numpy array holding translation distances for each
            dimension

        :returns: ``None``

    .. method:: translate_x(dist)

        Translate protein structure in the x direction.

        :param dist: amount of displacement in the x direction (Angstrom)
        :type dist: float

        :returns: ``None``

        Example::

            >>> a = Structure.PDBFile("/path/4N6W.pdb")
            >>> print a.get_all_atom_coords()[0] # before translation
            [-42.316  38.638  20.425]
            >>> a.translate_x(2.0)
            >>> print a.get_all_atom_coords()[0] # after translation
            [-40.316  38.638  20.425]

    .. method:: translate_y(dist)

        Translate protein structure in the y direction.

        :param dist: amount of displacement in the y direction (Angstrom)
        :type dist: float

        :returns: ``None``

    .. method:: translate_z(dist)

        Translate protein structure in the z direction.

        :param dist: amount of displacement in the z direction (Angstrom)
        :type dist: float

        :returns: ``None``

    .. method:: translate_origin_and_rotate(phi, theta, psi)

        Center the structure at the origin, rotate it with *angle_x*
        around the x axis (*angle_y* around y axis, etc.) and move
        it back to where it was. ??? confusion *angle_x* with *phi* ???

        :param phi: euler angle for rotation
        :type phi: float
        :param theta: euler angle for rotation
        :type theta: float
        :param psi: euler angle for rotation
        :type psi: float

        :returns: ``None``

    .. method:: move_to_new_position(new_coord)

        Move the geometric center of the structure to the
        supplied coordinates.

        :param new_coord: list/numpy array of the new coordinates

        :returns: ``None``

    .. method:: rotate_and_move_to_new_position(phi, theta, psi, new_coord)

        Center the structure at the origin, rotate it and
        move it to the new position *new_coord*.

        :param phi: euler angle for rotation
        :type phi: float
        :param theta: euler angle for rotation
        :type theta: float
        :param psi: euler angle for rotation
        :type psi: float
        :param new_coord: new coordination for the center of geometry

        :returns: ``None``

    .. method:: rotate(angle, axis)

        Rotate protein structure using the Rodrigues' rotation formula.

        :param degree: angle by which to rotate (in degrees)
        :type degree: float
        :param axis: axis around which to rotate (x = [1,0,0], 
           y = [0,1,0], z = [0,0,1])
        :type axis: array

        :returns: ``None``

        Example::

            >>> a = Structure.PDBFile("/path/4N6W.pdb")
            >>> print a.get_all_atom_coords()[0] # before rotation
            [-42.316,  38.638,  20.425]
            >>> a.rotate(+30.0, [1,0,0]) # x axis
            >>> print a.get_all_atom_coords()[0] # after rotation
            [-42.316  23.24898955  37.00756887]

    .. method:: rotateX(degree)

        :param degree: angle for rotation around the x axis (in degrees)
        :type  degree: float

        :returns: ``None``

        Example::

            >>> a = Structure.PDBFile("/path/4N6W.pdb")
            >>> print a.get_all_atom_coords()[0] # before rotation
            [-42.316,  38.638,  20.425]
            >>> a.rotateX(+30.0)
            >>> print a.get_all_atom_coords()[0] # after rotation
            [-42.316  23.24898955  37.00756887]

    .. method:: rotateY(degree)

        :param degree: angle for rotation around the y axis (in degrees)
        :type degree: float

        :returns: ``None``

    .. method:: rotateZ(degree)

        :param degree: angle for rotation around the z axis (in degrees)
        :type  degree: float

        :returns: ``None``

    .. method:: pRotate(phi, theta, psi)

        Apply euler angle rotation to the structure.

        :param phi: euler angle for rotation
        :type  phi: float
        :param theta: euler angle for rotation
        :type  theta: float
        :param psi: euler angle for rotation
        :type  psi: float

        :returns: ``None``

    .. method:: rotate_by_matrix(rot_matrix)

        :param rot_matrix: [3,3] numpy matrix to transform the coordinates by
            matrix multiplication

        :returns: ``None``

    .. method:: determineCenterOfMass()

        Method to determine the center of mass for the protein structure.

        :returns: ``None``

        .. warning::

            Chris: This method is empty.

    .. method:: determine_geometric_center()

        :returns: (np.array) a vector pointing to the geometric center of the
            whole structure

        Example::

            >>> a = Structure.PDBFile("/path/4N6W.pdb")
            >>> print a.determine_geometric_center()
            [-27.25387547  25.13667925  13.83051887]

    .. method:: determine_center_of_extremes_of_atoms(atoms)

        :param atoms: a list of atoms
        :type  atoms: :class:`atom`

        :returns: (np.array) a vector pointing to the geometric center of the
            coordination extremes of the given atom coordinates

    .. method:: determine_center_of_extremes()

        :returns: (np.array) a vector pointing to the geometric center of the
            coordination extremes of the whole structure

        Example::

            >>> a = Structure.PDBFile("/path/4N6W.pdb")
            >>> print a.determine_center_of_extremes()
            [-26.8235  24.865  13.304]

    .. method:: determine_max_diameter(atoms = None)

        :param atoms: a list of atoms, or all
            atoms of the structure if ``None``
        :type  atoms: :class:`atom`

        :returns: (float) the maximum diameter in Angstroem

        Example::

            >>> a = Structure.PDBFile("/path/4N6W.pdb")
            >>> print a.determine_max_diameter()
            52.4878883362

    .. method:: determine_radius(atoms = None)

        Determine the geometric center and calculate the minimal radius that
        encapsulates all atoms.

        :param atoms: a list of atoms (optional), or all
            atoms of the structure if ``None``
        :type atoms: :class:`Structure.atom` object

        :returns: the radius (Angstrom)
        :rtype: float

        .. warning::

            JN: In the source code, 
            ``center = self.determine_goometric_center_of_atoms(atoms)``
            should be replaced with
            ``center = self.determine_geometric_center()``
            but then, ValueError is raised

    .. method:: center()

        Translate the geometric center to the origin.

        :returns: ``None``

        Strictly equivalent to::

            >>> a.translate(-a.determine_geometric_center())

    .. method:: determine_coordinate_extremes(atoms = None)

        :param atoms: a list of :class:`atom` objects, otherwise it uses all
          atoms of this structure and calculates the extreme coordinates
          in each direction (optional)

        :returns: extreme values in each direction as a 3*2 array

    .. method:: get_radius_of_gyration()

        This method calculates the radius of gyration. It is the maximum
        distance of an atom to the geometrical center.

        :returns: radius of gyration
        :rtype: float

        Example::

            >>> a = Structure.PDBFile("/path/4N6W.pdb")
            >>> a.get_radius_of_gyration()
            28.848387249760471

    .. method:: clone(chain_id_list = None, res_id_list = None, res_type_list = None, atom_types_list = None)

        Return a clone of self. Through the list parameters specific items can
        be selected, namely the list of residues or certain types of residues
        (mutually exclusive). One-letter codes for the residues will be
        translated to three-letter codes.

        :param chain_id_list: list of chains to copy to the clone
        :param res_id_list: list of residues in each chain to copy to the clone
        :param res_type_list: types of residues to copy to the new clone
        :param atom_types_list: types of atoms to copy to the new clone

        :returns: a new :class:`PDBFile` / :class:`PQRFile` / :class:`LatFile` object
        :raises AttributeError: if both **res_id_list** and **res_type_list**
           were used, or if **self.what_am_i** is empty.

        Example::

            >>> a = Structure.PDBFile("/path/4N6W.pdb")
            >>> b = a.clone()
            >>> print b
            <epitopsy.Structure.PDBFile object at 0x2088a90>

    .. method:: get_residue_id_list(chain_id = None)

        Display all residue numbers as found in the PDB file.

        :param chain_id: which chain in the PDB file should be used
           (optional), if ``None``, it uses all the available chains
        :type chain_id: str

        :returns: a list with all residue ID's of the structure

        Example::

            >>> a = Structure.PDBFile("/path/4N6W.pdb")
            >>> print a.get_residue_id_list()
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, [...], 444]

    .. method:: get_res_id_aa_dict(chain_id)

        Display all residues from chain *chain_id* in a dictionary with
        residue number as key and amino acid one-letter code as value. The
        advantage of a dictionary over a list is that gaps in the sequence
        numbering are preserved. All non-amino acids are ignored.

        :param chain_id: which chain in the PDB file should be used
        :type chain_id: str

        :returns: a dictionary with residues number and one-letter code

        Example::

            >>> a = Structure.PDBFile("/path/4N6W.pdb")
            >>> a.get_res_id_aa_dict('A')
            Encountered the following non amino acids: ['FE', 'FLC', 'HOH']
            {1: 'M', 2: 'S', 3: 'L', 4: 'S', 5: 'N', 6: 'S', 7: 'S', 8: 'K', 9: 'V', 10: 'S', [...],  187: 'E'}

        .. note::

            JN: In the source code, it would be nice to replace
            ``.format(non_amino_acid_list)`` at line 741 by
            ``.format(sorted(set(non_amino_acid_list)))``.
            The complete list brings nothing.

    .. method:: contains_chain_break(chain_id = None)

        Tell if a chain break exist in *chain_id*, or in the whole structure if
        omitted. HETATM are skipped.

        :param chain_id: the chain id where to look for a break (optional),
           uses all available chains if ``None``

        :returns: either ``True`` (chain break) or ``False`` (no chain break).
        :rtype: bool

        .. note::

            JN: In the source code, ``return True`` could be replaced with
            ``return "Break between {0} and {1}".format(res_id_list[-1],res_id)``
            (a string always evaluate to ``True``) to indicate the user where the
            break is. This string could be displayed in the terminal or stored
            in a log file.

    .. method:: get_res_id_array_mapping()

        Remove gaps in the residue sequence and return the new mapping
        in a dictionary, with the old residue id's as key and the new
        id's as value. The new mapping starts at zero and is suitable
        for use as an index for an array.

        :returns: a dictionary with old residue id's as key and newly mapped
           id's as value

        Example::

            >>> a = Structure.PDBFile("/path/4N6W.pdb")
            >>> a.get_res_id_array_mapping()
            {1: 0, 2: 1, 3: 2, 4: 3, 5: 4, 6: 5, 7: 6, 8: 7, 9: 8, 10: 9, ..., 444: 326}

    .. method:: get_residue_names_from_res_id_list(res_id_list, chain_id = None)

        Display the three-letter code of the residues given in *res_id_list*
        from chain *chain_id* (if ``None``, take the first chain).

        :param res_id_list: residue id's from which one wants the names
        :type res_id_list: list
        :param chain_id: if ``None``, it uses the all available chains (optional)
        :type chain_id: str

        :returns: a list with the three-letter codes of the residues matching
           the given criteria

        Example::

            >>> a = Structure.PDBFile("/path/4N6W.pdb")
            >>> a.get_residue_names_from_res_id_list([1,2]) # residue id starts at 1 in the PDB file
            ['MET', 'SER']

    .. method:: get_residues(chain_id = None, res_id_list = None)

        Display all residues from *res_id_list* contained in chain *chain_id*
        as :class:``Structure.Residue``. If *res_id_list* is ``None``, display
        all residues from chain *chain_id*. If *chain_id* is ``None``, return
        an error excepted when there is only one chain in the PDB file.
        
        This method returns a list with residue objects of the residue ids in
        the given list, if ``None`` is given, it returns all residues.

        :param chain_id: if ``None``, it uses the all available chains (optional)
        :type chain_id: str
        :param res_id_list: residue id's from which one wants the objects
        :type res_id_list: list

        :returns: a list of residues matching the given criteria
        :rtype: :class:`Structure.Residue` object

        Example::

            >>> a = Structure.PDBFile("/path/4N6W.pdb")
            >>> a.get_residues('A')[0:2] # residue id starts at 0 in a list
            [<Residue MET het=  resseq=1 icode= >, <Residue SER het=  resseq=2 icode= >]

        .. warning::

            JN: this method has a bug, please consider replacing
            ``res_list.append(res.resname)`` by ``res_list.append(res)``
            to get the objects instead of just the names!

    .. method:: get_atoms_by_id(atom_id_list)

        This function returns atoms by their corresponding number from
        the pdb-file.

        :param atom_id_list: atom id numbers
        :type atom_id_list: list

        :returns: a list of atoms matching the given criteria
        :rtype: :class:`Structure.atom` object

        Example::

            >>> a = Structure.PDBFile("/path/4N6W.pdb")
            >>> a.get_atoms_by_id([1,2,3])
            [<Atom N>, <Atom CA>, <Atom C>]

    .. method:: get_atoms_close_to_reference(reference_point, max_radius, min_radius = 0, atoms = None)

        Find all atoms lying in the range [**min_radius**, **max_radius**] to
        **reference_point**.

        :param reference_point: coordinates of the reference point
        :type reference_point: :class:`numpy.array`
        :param max_radius: maximal distance from the reference (Angstroem)
        :type max_radius: float
        :param min_radius: minimal distance from the reference (Angstroem)
        :type min_radius: float
        :param atoms: list of the atoms to which the results should be
           restricted (optional), if ``None`` it uses all
           atoms from the protein
        :type atoms: :class:`Structure.atom`

        :returns: a list of atoms close to the given reference point
        :rtype: :class:`Structure.atom` object

        Example::

            >>> a = Structure.PDBFile("/path/4N6W.pdb")
            >>> # retrieve the iron catalytic center
            >>> fe_list = a.get_atoms_of_type("FE")               # fe_list = [<Atom FE>, <Atom FE>]
            >>> coord = a.get_coords_from_atom_list([fe_list[1]]) # coord = numpy.array([-28.139,20.588,11.849])
            >>> # get all 6 atoms chelating Fe, excepted Fe
            >>> a.get_atoms_close_to_reference(coord, 2.5, 0.1)
            [<Atom NE2>, <Atom NE2>, <Atom OD2>, <Atom OD1>, <Atom OG1>, <Atom O>]

        Biopython equivalent::

            >>> from Bio.PDB import PDBParser
            >>> from Bio.PDB import NeighborSearch
            >>> p = PDBParser()
            >>> a = p.get_structure('4N6W.pdb', '/path/4N6W.pdb')
            >>> fe_list = []
            >>> atom_list = []
            >>> for residue in struct[0]['A']:
            ...     for atom in residue:
            ...         if atom.name == "FE":
            ...             fe_list.append(atom) # fe_list = [<Atom FE>, <Atom FE>]
            ...         elif atom.name[0] <> "H":
            ...             atom_list.append(atom)
            >>> ns = NeighborSearch(atom_list)
            >>> rd = 2.5
            >>> coord = fe_list[1].get_coord() # coord = numpy.array([-28.139,20.588,11.849])
            >>> print sorted(ns.search(coord,rd,'A'))
            [<Atom NE2>, <Atom NE2>, <Atom OD2>, <Atom OD1>, <Atom OG1>, <Atom O>]

    .. method:: find_chain_contacts(chain1, chain2, max_distance)

        Find atoms of **chain1** which are within **max_distance** of **chain2**.

        :param chain1: first chain id
        :type chain1: str
        :param chain2: second chain id
        :type chain2: str
        :param max_distance: maximal distance to include atoms in the calculation
        :type max_distance: float

        :returns: a list of atoms from **chain1** for which distance to
            **chain2** is smaller than **max_distance**
        :rtype: :class:`Structure.atom` object

    .. method:: get_sequence()

        Display the one-letter code sequence of each chain of the PDB file in
        a dictionary. Non standard amino acids will not be returned.

        :returns: a dictionary with the amino acid sequence of each chain

        Example::

            >>> a = Structure.PDBFile("/path/4N6W.pdb")
            >>> a.get_sequence()
            {'A': 'MSL[...]LSE'}

    .. method:: snap_vdw_to_box(box_mesh_size, box_dim, box_offset, warning = True, vdw_radii = 2, increase_vdw_by = 0)

        Snap a structure to a given dxbox. If the structure is a pqr, it uses
        the supplied vdw radii otherwise it uses  **vdw_radii** for each
        atom. If any coordinate lies outside the box an error will be printed
        to the standard output, unless **warning** is set to ``False``.

        :param box_mesh_size: mesh size [m,m,m]
        :type box_mesh_size: list
        :param box_dim: [x,y,z]
        :typebox_dim: list
        :param box_offset: [x_o,y_o,z_o]
        :type box_offset: list
        :param warning: if ``True`` print an error when the structure does not
           fit completely into the given box dimensions
        :type warning: bool
        :param vdw_radii: if this is a pdb file there are no other radii
            available (in Angstroem)
        :type vdw_radii: float
        :param increase_vdw_by: can be used to blow up the radii of each atom
           (in Angstroem)
        :type increase_vdw_by: float

        :returns: Numpy array with 0's outside and 1's inside the protein.

    .. method:: get_rmsd_rotation_translation_from_superposition(pdb_to_rotate, atom_types = None)

        Fit the :class:`Struture.Structure` object **pdb_to_rotate** onto self
        Return a dictionary containing a *rmsd* value, a *rotation* matrix
        and a *translation* value.
        The parameter **atom_types** restricts the fitting on a particular set
        of atoms, if ``None`` is supplied, all atoms will be used.

        I guess the units are Angstroem.

        :param pdb_to_rotate: object to superimpose onto this object
        :type pdb_to_rotate: :class:`Struture.Structure`
        :param atom_types: restrict the fitting to specific types (e.g. ['CA'],
           ['CA','N']), or fit all atoms if ``None`` (optional)
        :type atom_types: list

        :returns: a dictionary with following keys:
            * 'rmsd' : root mean square deviation
            * 'rotation' : rotation matrix
            * 'translation' : translation vector
        :raises ValueError: if there is an atom missmatch between the two pdb's
           or if at least one of the atom lists is empty

    .. method:: superimpose_given_pdb_onto_self(pdb_to_superimpose, atom_types = None)

        Fit the :class:`Struture.Structure` object **pdb_to_superimpose** onto
        self. The atomic coordinates are updated in the process.

        :param pdb_to_superimpose: object to superimpose onto this object
        :type pdb_to_superimpose: :class:`Struture.Structure`
        :param atom_types: restrict the fitting to specific types (e.g. ['CA'],
           ['CA','N']), or fit all atoms if ``None`` (optional)
        :type atom_types: list

        :returns: ``None``
        :raises ValueError: if there is an atom missmatch between the two pdb's
           or if at least one of the atom lists is empty

    .. method:: superimpose_self_onto_given_pdb(pdb_to_superimpose, atom_types = None)

        Fit self onto the :class:`Struture.Structure` object
        **pdb_to_superimpose** and update self's atomic coordinates.

        :param pdb_to_superimpose: object on which to superimpose self
        :type pdb_to_superimpose: :class:`Struture.Structure`
        :param atom_types: restrict the fitting to specific types (e.g. ['CA'],
           ['CA','N']), or fit all atoms if ``None`` (optional)
        :type atom_types: list

        :returns: ``None``

    .. method:: get_dxbox_dim(box_mesh_size, extend = None, cubic_box = True, nlev = 4)

        Return the dimensions of a DXbox. The edges are calculated using the
        protein maximal diameter in each direction, **extend** if given,
        and the grid resolution **box_mesh_size**, using the formula:

        :math:`a[i] = \frac{\displaystyle protein\_diameter[i] +
        2 \cdot extend[i]}{\displaystyle box\_mesh\_size[i]}`

        with *i* = {x,y,z} the coordinates. If **cubic_box** is ``True``,
        all edges have the same length. The lengths are rounded up to the
        next value of *n* calculated with **nlev** according to the formula:

        :math:`n = c \cdot 2^{nlev + 1} + 1`

        with `nlev <http://www.poissonboltzmann.org/apbs/user-guide/running-apbs/input-files/elec-input-file-section/elec-keywords/nlev>`_
        the multilevel hierarchy of the calculation. See `dime
        <http://www.poissonboltzmann.org/apbs/user-guide/running-apbs/input-files/elec-input-file-section/elec-keywords/dime>`_
        for an explanation of this parameter.

        The center of the box is the geometric center of
        the protein if not otherwise specified.

        .. note:

            Chris: The calculation is copied from :class:`InFile`. If there
            have been changes this result might be wrong!

        :param box_mesh_size: resolution of the grid (Angstroems)
        :type  box_mesh_size: np.array
        :param extend: extension of the box dimensions (Angstroems) (optional)
        :type  extend: float
        :param cubix_box: use a cubic box if ``True`` (optional)
        :type  cubic_box: bool
        :param nlev: depth of the multilevel hierarchy used by the multigrid
            solver (optional)
        :type  nlev: int

        :returns: (np.array) dimensions of the DXbox

    .. method:: get_dxbox_offset(box_mesh_size, box_dim, box_center)

        Returns the offset for the given dimensions.

        :param box_mesh_size: dimensions of the mesh (Angstroems)
        :type  box_mesh_size: array
        :param box_dim: dimensions of the box (Angstroems)
        :type  box_dim: array
        :param box_center: center of the box (Angstroems)
        :type  box_center: array

        :returns: a list box_offset: [x_o,y_o,z_o].

    .. method:: get_hydrophobic_potential(box_mesh_size, box_dim, box_offset)

        Calculate the hydrophobic potential of a protein.
        This method uses a simplified model for the hydrophobic potential.
        The charges are taken from the Kyte and Doolittle hydrophobicity
        scale. For each residue the center is calculated and the potential
        is modelled as:

            :math:`\phi = \sum \left ( hydrophobic\_charge \times e^{( - distance )} \right )`

        .. seealso::

            Kyte, Doolittle, *A simple method for displaying the hydropathic
            character of a protein*, *J. Mol. Biol.* **1982**, *157*, 105-132.

        :param box_mesh_size: meshsize of the grid
        :param box_dim: dimension of the grid
        :param box_offset: offset of the grid

        :returns: Numpy array with the hydrophobic potential for each grid point.

        .. function:: get_residue_center(res)

            Subroutine of :func:`Structure_Template.get_hydrophobic_potential`.

    .. method:: get_vdw_hull(box_mesh_size, box_dim, box_offset, vdw_radii = 2, increase_vdw_by = 0)

        Get the van der Waals hull of the protein.

        :param box_mesh_size: meshsize of the grid
        :param box_dim: dimension of the grid
        :param box_offset: offset of the grid
        :param vdw_radii: atom radii, if this is not a pqr
        :param increase_vdw_by: extend each radius by this value

        :returns: Numpy array with 1's at the hull grid points and 0's everywhere else.

    .. method:: get_num_of_overlap_atoms_with_given_structure(other_structure_object, energy_cutoff = 1., vdw_radii = 2.)

        Calculate the number of atoms in this structure object that are
        overlapping with the given structure object.

        :param other_structure_object: structure which might have an overlapp
           with this one
        :param energy_cutoff: cutoff for the lennard jones potential to decide
           if there is an overlapp or not
        :param vdw_radii: if this is pdb object there are no radii information
           available

        :returns: counts of atoms in this structure that overlapp with
          the given structure object

    .. method:: get_contact_list(cutoff = 5.)

        Find steric contacts between residues of a protein. A contact is found
        when the interatomic distance of at least two atoms taken from two
        different residues *i* and *j* is inferior to **cutoff** (Angstroem).

        :param cutoff: interatomic distance
        :type cutoff: float

        :returns: a list of residues id's as found in the PDB, with 1 for a
           contact and 0 for no contact: ``[[i-1,j-1,0], [i,j,1], [i+1,j+1,0]]``

        Example::

            >>> a = Structure.PDBFile("/path/4N6W.pdb")
            >>> for contact in a.get_contact_list(cutoff = 5.0):
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

        .. note::

            Perhaps we could replace the lists ``[i,j,0]`` by tuples ``(i,j,0)``

        .. function:: residue_contact(res_i, res_j, cutoff)

            Subroutine of :func:`Structure_Template.get_contact_list`.

Subclasses
^^^^^^^^^^

PDB object
""""""""""

.. class:: PDBFile(Structure_Template)

    This class is used for all reading/writing operations on PDB files.

    .. warning::

        where are all the ``self`` gone?

    .. method:: get_pqr_structure(new_pqr_path = None, force_field = "amber", pdb2pqr_argv = ``None``, pdb2pqr_path = "pdb2pqr", add_ions = True, add_chain = True)

        Call pdb2pqr to replace the b-factor and the occupancy information in
        the pdb file with the charges and the vdw radii. If this object has no
        ``structure_path`` property (i.e. it is ``None``), then an error is
        raised.

        If the pdb contains Ca\ :sup:`2+`, Zn\ :sup:`2+` or
        SO\ :sub:`4`:sup:`2--`,
        these ions will be added unless stated otherwise.

        :param new_pqr_path: path for the new pqr file. If ``None`` is given,
            ``structure_path`` will be used by changing its extension to ".pqr"
        :type new_pqr_path: str
        :param force_field: forcefield from which charges and radii should be
            taken. Default is "amber"
        :type force_field: str
        :param pdb2pqr_argv: can contain additional arguments to pdb2pqr
           (e.g. ['--assign-only'], oder ['--noopt']) (optional)
        :type pdb2pqr_argv: list
        :param pdb2pqr_path: path to the executable of PDB2PQR
        :type pdb2pqr_path: str
        :param add_ions: add Ca, Zn, SO\ :sub:`4` ions if they are in the pdb.
            The Ca, Zn and SO\ :sub:`4` atoms should have a residue name that
            fits their type (CA, ZN, SO4).
        :type add_ions: bool

        :returns: a :class:`PQRFile` object
        :raises AttributeError: if **structure_path** is empty, or if
           **new_pqr_path** cannot be read or generated from the PDB path.
        :raises NameError: if **new_pqr_path** refers to a file already existing

    .. method:: _read_file()

        This method reads the data from the file given in ``structure_path``.
        Subroutine of the constructor.

    .. method:: save_to_file(path)

        This method writes the structure to a given path.
        If the structure is a lattice, it will try to calculate the correct
        connections. The Assumption of a lattice is made on the number of
        residues and the number of overall atom coordinates.
        If both numbers are equal, it is very probable that this is a lattice
        protein. This only works, if there is just one chain!

        Notice: We do not work with multiple models!

        :param path: path for the new pdb file.
        :type path: str

        :returns: ``None``
        :raises AttributeError: if **path** does not end with '.pdb' or '.pdb.gz'
        
        .. note:

            If **path** points to a file, it will be overwritten without
            confirmation message.

    .. method:: _save_to_file(f)

        Subroutine of :func:`save_to_file`.

    .. method:: _get_atom_line(chain, res, atom)

        Retrieve one line from the stored structure.
        Subroutine of :func:`_save_to_file`.

        :param chain: chain
        :type chain: str
        :param res: residue
        :type res: str
        :param atom: atom
        :type atom: str

        :returns: a formatted string

PQR object
""""""""""

.. class:: PQRFile(Structure_Template)

    This class is used for all reading/writing operations on PQR files. It
    shares many similarities with :class:`PDBFile`.

    .. method:: _read_structure()

        This method reads the data from the file given in ``structure_path``.
        Subroutine of the constructor.

    .. method:: save_to_file(path)

        This method writes the structure to a given path.
        If the structure is a lattice, it will try to calculate the correct
        connections. The Assumption of a lattice is made on the number of
        residues and the number of overall atom coordinates.
        If both numbers are equal, it is very probable that this is a lattice
        protein. This only works, if there is just one chain!

        Notice: We do not work with multiple models!

        :param path: path for the new pqr file.
        :type path: str

        :returns: ``None``
        :raises AttributeError: if **path** does not end with '.pqr' or '.pqr.gz'
        
        .. note:

            If **path** points to a file, it will be overwritten without
            confirmation message.

    .. method:: _save_to_file(,f)

        Subroutine of :func:`save_to_file`.

    .. method:: _get_atom_line(chain, res, atom)

        Retrieve one line from the stored structure.
        Subroutine of :func:`_save_to_file`.

        :param chain: chain
        :type chain: str
        :param res: residue
        :type res: str
        :param atom: atom
        :type atom: str

        :returns: a formatted string

    .. method:: snap_esp_to_dxbox(dxbox, warning = True)

        Snap the charge for each atom of this pqr structure to a
        given dxbox. The new array contains the charge in Coulombs.

        If any coordinate lies outside the box an error will be printed to the
        standard output.

        :param dxbox: a dxbox object
        :param warning: either ``True`` or ``False``

        :returns: a numpy array with the charges in units of e at the center
          of each atom

        .. warning::

            If you use a neutral probe and your mesh size is to
            large the dipole effect is not visible!

Lattice object
""""""""""""""

.. class:: latFile(Structure_Template)

    This class is used for all reading/writing operations on lattice files.

    .. method:: _make_structure(seq, fold, lattice, mesh_size)

        This method reads the data from the given file.

    .. method:: save_to_file(path)

        This method writes the structure to a given path.
        If the structure is a lattice, it will try to calculate the
        correct connections. The Assumption of a lattice is made on the
        number of residues and the number of overall atom coordinates.
        If both numbers are equal, it is very probable that this is a lattice
        protein. This only works, if there is just one chain!

        Notice: We do not work with multiple models!

        :param path: path for the new lattice file.
        :type path: str

        :returns: ``None``
        
        .. note:

            If **path** points to a file, it will be overwritten without
            confirmation message.

    .. method:: _get_atom_line(chain, res, atom)

        Formats the line for writing.

    .. method:: _convert_fcc_fold(fold, mesh_size)

        Docstring missing.

    .. method:: _convert_sc_fold(fold, mesh_size)

        Docstring missing.

Main entity class
^^^^^^^^^^^^^^^^^

.. class:: entity(object)

    Docstring missing.

    .. method:: __len__()

        Return the number of children.

    .. method:: __getitem__(id)

        Return the child with given id.

    .. method:: __delitem__(id)

        Remove a child.

    .. method:: __iter__()

        Iterate over children.

    .. method:: get_level()

        :returns: level in hierarchy:

            * A - atom
            * R - residue
            * C - chain
            * M - model
            * S - structure

    .. method:: set_parent(entity)

        Set the parent Entity object.

    .. method:: detach_parent()

        Detach the parent.

    .. method:: detach_child(id)

        Remove a child.

    .. method:: add(entity)

        Add a child to the Entity.

    .. method:: get_iterator()

        :returns: iterator over children

    .. method:: get_list()

        :returns: a copy of the list of children

    .. method:: has_id(id)

        :returns: ``True`` if a child with given **id** exists

    .. method:: get_parent()

        :returns: the parent Entity object

    .. method:: get_id()

        :returns: the id

    .. method:: get_full_id()

        Return the full id.

        The full id is a tuple containing all id's starting from
        the top object (Structure) down to the current object. A full id for
        a Residue object e.g. is something like:

        ("1abc", 0, "A", (" ", 10, "A"))

        This corresponds to:

        * Structure with id "1abc"
        * Model with id 0
        * Chain with id "A"
        * Residue with id (" ", 10, "A")

        The Residue id indicates that the residue is not a hetero-residue
        (or a water) beacuse it has a blank hetero field, that its sequence
        identifier is 10 and its insertion code "A".


Subclasses
^^^^^^^^^^

Structure
"""""""""

.. class:: Structure(Entity)

    The Structure class contains a collection of Model instances.

    .. method:: __repr__()

        Docstring missing.

    .. method:: _sort(m1, m2)

        Sort models.

        This sorting function sorts the Model instances in the Structure instance.
        The sorting is done based on the model id, which is a simple int that
        reflects the order of the models in the PDB file.

        :param m1: model instance
        :param m2: model instance

    .. method:: get_chains()

        Docstring missing.

    .. method:: get_residues()

        Docstring missing.

    .. method:: get_atoms()

        Docstring missing.

Model
"""""

.. class:: Model(Entity)

    The object representing a model in a structure. In a structure
    derived from an X-ray crystallography experiment, only a single
    model will be present (with some exceptions). NMR structures
    normally contain many different models.

    .. attribute:: id

        identifiant

    .. attribute:: serial_num

        serial number

    .. method:: _sort(c1, c2)

        Sort the Chains instances in the Model instance.

        Chain instances are sorted alphabetically according to their
        chain id. Blank chains come last, as they often consist of waters.

        :param c1: Chain object
        :param c2: Chain object

    .. method:: __repr__()

        Docstring missing.

    .. method:: get_residues()

        Docstring missing.

    .. method:: get_atoms()

        Docstring missing.

Chain
"""""

.. class:: Chain(Entity)

    Docstring missing.

    .. method:: _sort(r1, r2)

        Sort function for residues in a chain

        Residues are first sorted according to their hetatm records.
        Protein and nucleic acid residues first, hetatm residues next,
        and waters last. Within each group, the residues are sorted according
        to their resseq's (sequence identifiers). Finally, residues with the
        same resseq's are sorted according to icode.

        :param r1: Residue object
        :param r2: Residue object

    .. method:: _translate_id(id)

        A residue id is normally a tuple (hetero flag, sequence identifier,
        insertion code). Since for most residues the hetero flag and the
        insertion code are blank (i.e. " "), you can just use the sequence
        identifier to index a residue in a chain. The _translate_id method
        translates the sequence identifier to the (" ", sequence identifier,
        " ") tuple.

        :param id: residue resseq
        :type id: int

    .. method:: __getitem__(id)

        Return the residue with given id.

        The id of a residue is (hetero flag, sequence identifier, insertion code).
        If id is an int, it is translated to (" ", id, " ") by the _translate_id
        method.

        :param id: residue resseq
        :type id: int or string?

    .. method:: __delitem__(id)

        :param id: residue resseq
        :type id: int or string?

    .. method:: __repr__()

        Docstring missing.

    .. method:: get_unpacked_list()

        Return a list of undisordered residues.

        Some Residue objects hide several disordered residues
        (DisorderedResidue objects). This method unpacks them,
        ie. it returns a list of simple Residue objects.

    .. method:: has_id(id)

        The id of a residue is (hetero flag, sequence identifier, insertion code).
        If id is an int, it is translated to (" ", id, " ") by the _translate_id
        method.

        :param id: residue resseq
        :type id: int or string?

        :returns: 1 if a residue with given id is present

    .. method:: get_atoms()

        Docstring missing.

Residue
"""""""

.. class:: Residue(Entity)

    Represents a residue. A Residue object stores atoms.

    .. method:: __repr__()

        Docstring missing.

    .. method:: _sort(a1, a2)

        Sort the Atom objects.

        Atoms are sorted alphabetically according to their name,
        but N, CA, C, O always come first.

        Arguments:
        :param a1: Atom object
        :param a2: Atom object

    .. method:: add(atom)

        Add an Atom object.

        Checks for adding duplicate atoms, and raises a
        PDBConstructionException if so.

    .. method:: sort()

        Docstring missing.

    .. method:: flag_disordered()

        Set the disordered flag.

    .. method:: is_disordered()

        Return 1 if the residue contains disordered atoms.

    .. method:: get_resname()

        Docstring missing.

    .. method:: get_unpacked_list()

        Returns the list of all atoms, unpack DisorderedAtoms.

    .. method:: get_segid()

        Docstring missing.

Atom
""""

.. class:: atom(object)

    Atom object.

    The Atom object stores atom name (both with and without spaces),
    coordinates, B factor, occupancy, alternative location specifier
    and (optionally) anisotropic B factor and standard deviations of
    B factor and positions.

    .. attribute:: name

        (str) atom name (eg. "CA"). Note that spaces are normally stripped.

    .. attribute:: coord

        (Numeric array (Float0, size 3)) atomic coordinates (x,y,z).

    .. attribute:: bfactor

        (float) isotropic B factor.

    .. attribute:: occupancy

        (float) occupancy (0.0-1.0).

    .. attribute:: altloc

        (str) alternative location specifier for disordered atoms

    .. attribute:: fullname

        (str) full atom name, including spaces, e.g. " CA ".
        Normally these spaces are stripped from the atom name.

    .. attribute:: element

        (uppercase string or ``None`` if unknown) atom element,
        e.g. "C" for Carbon, "HG" for mercury

    .. method:: __repr__()

        Print Atom object as <Atom atom_name>.

    .. method:: __sub__(other)

        Calculate distance between two atoms.

        :param other: the other atom

    .. method:: set_serial_number(n)

        Docstring missing.

    .. method:: set_bfactor(bfactor)

        Docstring missing.

    .. method:: set_coord(coord)

        Docstring missing.

    .. method:: set_altloc(altloc)

        Docstring missing.

    .. method:: set_occupancy(occupancy)

        Docstring missing.

    .. method:: set_sigatm(sigatm_array)

        Set standard deviation of atomic parameters.

        The standard deviation of atomic parameters consists
        of 3 positional, 1 B factor and 1 occupancy standard
        deviation.

        :param sigatm_array: standard deviations of atomic parameters
        :type sigatm_array: numeric array (length 5)

    .. method:: set_siguij(siguij_array)

        Set standard deviations of anisotropic temperature factors.

        :param siguij_array: standard deviations of anisotropic temperature factors
        :type siguij_array: numeric array (length 6)

    .. method:: set_anisou(anisou_array)

        Set anisotropic B factor.

        :param anisou_array: anisotropic B factor
        :type anisou_array: numeric array (length 6)

    .. method:: set_element(element)

        Set element.

    .. method:: set_info_1(info_1)

        This is just a more abstract method which I will use instead of
        'get_bfactor' and 'get_occupancy'.

        * pdb: occupancy
        * pqr: charge

    .. method:: set_info_2(info_2)

        This is just a more abstract method which I will use instead of
        'get_bfactor' and 'get_occupancy'.

        * pdb: b_factor
        * pqr: vdw radius

    .. method:: flag_disorder()

        Set the disordered flag to 1.

        The disordered flag indicates whether the atom is disordered or not.

    .. method:: is_disordered()

        Return the disordered flag (1 if disordered, 0 otherwise).

    .. method:: set_parent(parent)

        Set the parent residue.

        :param parent: Residue object

    .. method:: detach_parent()

        Remove reference to parent.

    .. method:: get_sigatm()

        :returns: standard deviation of atomic parameters

    .. method:: get_siguij()

        :returns: standard deviations of anisotropic temperature factors

    .. method:: get_anisou()

        :returns: anisotropic B factor

    .. method:: get_parent()

        :returns: parent residue

    .. method:: get_serial_number()

        Docstring missing.

    .. method:: get_name()

        :returns: atom name

    .. method:: get_id()

        :returns: the id of the atom (which is its atom name)

    .. method:: get_full_id()

        Return the full id of the atom.

        The full id of an atom is the tuple
        (structure id, model id, chain id, residue id, atom name, altloc)

    .. method:: get_coord()

        :returns: atomic coordinates

    .. method:: get_bfactor()

        :returns: B factor

    .. method:: get_occupancy()

        :returns: occupancy

    .. method:: get_fullname()

        :returns: the atom name, including leading and trailing spaces

    .. method:: get_altloc()

        :returns: alternative location specifier

    .. method:: get_level()

        Docstring missing.

    .. method:: transform(rot, tran)

        Apply rotation and translation to the atomic coordinates.

        :param rot: a right multiplying rotation matrix
        :type rot: 3x3 numeric array

        :param tran: the translation vector
        :type tran: size 3 Numeric array

        Example::

                >>> rotation = rotmat(pi, Vector(1,0,0))
                >>> translation = array((0,0,1), 'f')
                >>> atom.transform(rotation, translation)

    .. method:: get_info_1()

        This is just a more abstract method which I will use instead of
        'get_bfactor' and 'get_occupancy'.

        * pdb: occupancy
        * pqr: charge

    .. method:: get_info_2()

        This is just a more abstract method which I will use instead of
        'get_bfactor' and 'get_occupancy'.

        * pdb: bfactor
        * pqr: radius

    .. method:: get_element()

        :returns: element

    .. method:: set_charge(charge)

        Docstring missing.

    .. method:: get_charge()

        Docstring missing.

