
:mod:`Structure` --- YYYYYYY
======================================================

.. module:: Structure
   :synopsis: YYYYYYYY.
.. moduleauthor:: Christoph Wilms <christoph.wilms@uni-due.de>
.. sectionauthor:: Jean-Noel Grad <jean-noel.grad@uni-due.de>


This module provides YYYYY.


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

.. _Structure-syntax:

Module Syntax
-------------

Empty.

.. _contents-of-module-Structure:

Module Contents
---------------

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

    .. method:: get_coords_from_atom_list(atom_list)

        :param atom_list: list of atom objects

        :returns: a list of atom coordinates for the specified atom objects

    .. method:: get_hetero_atoms()

        :returns: a list of all hetero atoms

    .. method:: get_hetero_atoms_coords()

        :returns: a list of all hetero atom coordinates

    .. method:: get_amino_atoms()

        :returns: a list of all amino acid atoms

    .. method:: get_amino_atoms_coords()

        :returns: a list of all amino acid atom coordinates.

    .. method:: get_all_atoms()

        :returns: a list of all atom coordinates

    .. method:: get_info_1(atoms=None)

        :returns: a list of all atom information 1 values. For a pdb this is
            the occupancy and for a pqr it is the charge

    .. method:: get_info_2(atoms=None)

        :returns: a list of all atom information 2 values. For a pdb this is
            the temperature factor and for a pqr it is the radius

    .. method:: get_chain_ids()

        :returns: a list of all chain ids in this structure

    .. method:: get_first_res_id()

        :returns: the integer number of the first amino acid

    .. method:: get_atoms_of_type(atom_type)

        :param atom_type: type of atom to return (e.g. 'CA')

        :returns: all atom objects of a certain type

    .. method:: transform(T)

        Transform the pdb structure with the given matrix.

        :param T: [3,3] numpy matrix  to transform the coordinates by
            matrix multiplication

        :returns: ``None``

    .. method:: translate(transVector)

        Method to translate protein structure by the given vector.

        :param transvector: Numpy array holding translation distances for each
            dimension

        :returns: ``None``

    .. method:: translate_x(dist)

        Method to translate protein structure in x direction.

        :param dist: amount of displacement in x direction.

        :returns: ``None``

    .. method:: translate_y(dist)

        Method to translate protein structure in y direction.

        :param dist: amount of displacement in y direction

        :returns: ``None``

    .. method:: translate_z(dist)

        Method to translate protein structure in z direction.

        :param dist: amount of displacement in z direction

        :returns: ``None``

    .. method:: translate_origin_and_rotate(phi, theta, psi)

        This methods centers the structure at the origin, rotates it with
        angle_x around the x axis (angle_y around y axis, etc.) and moves
        it back to where it was.

        :param phi: euler angle for rotation
        :param theta: euler angle for rotation
        :param psi: euler angle for rotation

        :returns: ``None``

    .. method:: move_to_new_position(new_coord)

        This method moves the geometric center of the structure to the
        supplied coordinates.

        :param new_coord: list/numpy array of the new coordinates

        :returns: ``None``

    .. method:: rotate_and_move_to_new_position(phi, theta, psi, new_coord)

        This method centers the structure at (0,0,0), rotates it and the
        moves it to a new position.

        :param phi: euler angle for rotation
        :param theta: euler angle for rotation
        :param psi: euler angle for rotation
        :param new_coord: new coordination for the center of geometry

        :returns: ``None``

    .. method:: rotate(angle, axis)

        Method to rotate protein structure.
        For the rotation I use the Rodrigues' rotation formula.

        :param degree: angle by which to rotate
        :param axis: axis around which to rotate

        :returns: ``None``

    .. method:: rotateX(degree)

        :param degree: angle for rotation around the x axis

        :returns: ``None``

    .. method:: rotateY(degree)

        :param degree: angle for rotation around the y axis

        :returns: ``None``

    .. method:: rotateZ(degree)

        :param degree: angle for rotation around the z axis

        :returns: ``None``

    .. method:: pRotate(phi, theta, psi)

        Apply euler angle rotation to the structure.

        :param phi: euler angle for rotation
        :param theta: euler angle for rotation
        :param psi: euler angle for rotation

        :returns: ``None``

    .. method:: rotate_by_matrix(rot_matrix)

        :param rot_matrix: [3,3] numpy matrix to transform the coordinates by
            matrix multiplication

        :returns: ``None``

    .. method:: determineCenterOfMass()

        Method to determine the center of mass for the protein structure.

        :returns: ``None``

    .. method:: determine_geometric_center()

        :returns: a vector pointing to the geometric center of the structure

    .. method:: determine_center_of_extremes_of_atoms(atoms)

        :param atoms: list of atom objects

        :returns: a vector pointing to the geometric center of the coordination
            extremes of the given atom coordinates

    .. method:: determine_center_of_extremes()

        :returns: a vector pointing to the geometric center of the coordination
            extremes of this structure

    .. method:: determine_max_diameter(atoms = None)

        :param atoms: a list of atom objects, otherwise it uses all
          atoms of this structure and calculates the maximum diameter (optional)

        :returns: a float number of the maximum diameter

    .. method:: determine_radius(atoms = None)

        Determine the geometric center and calculate the minimal radius that
        encapsulates all atoms.

        :param atoms: optional a list of atom objects, otherwise it uses all
            atoms of this structure and calculates the radius

        :returns: a float number of the radius

    .. method:: center()

        Translate geometric center to (0.0/0.0/0.0).

        :returns: ``None``

    .. method:: determine_coordinate_extremes(atoms = None)

        :param atoms: optional a list of atom objects, otherwise it uses all
          atoms of this structure and calculates the extreme coordinates
          in each  direction.

        :returns: extreme values in each direction as a 3*2 array

    .. method:: get_radius_of_gyration()

        This method calculates the radius of gyration. It is the maximum
        distance of an atom to the geometrical center.

        :returns: a float number as the radius of gyration

    .. method:: clone(chain_id_list = None, res_id_list = None, res_type_list = None, atom_types_list = None)

        This method returns a clone of this structure. Through the list
        parameters specific items can be selected. If one supplies one letter
        codes for the residues they will be translated to three letter codes!

        :param chain_id_list: list of chains to copy to the clone
        :param res_id_list: list of residues in each chain to copy to the clone
        :param res_type_list: types of residues to copy to the new clone
        :param atom_types_list: types of atoms to copy to the new clone

        :returns: a new PDBFile / PQRFile / LatFile object.

    .. method:: get_residue_id_list(chain_id = None)

        :param chain_id: optionally, if ``None``, it uses the all available chains

        :returns: a list with all residue ID's of the structure

    .. method:: get_res_id_aa_dict(chain_id)

        :returns: a dictionary of the first chain, which contains the residue id as
            key and the corresponding amino acid as value

    .. method:: contains_chain_break(chain_id = None)

        :param chain_id: if ``None``, it uses the all available chains (optional)

        :returns: either ``True (chain break) or ``False`` (no chain break).

    .. method:: get_res_id_array_mapping()

        :returns: a dictionary with residue ids as keys and values, which can be
            used as an index in an array

    .. method:: get_residue_names_from_res_id_list(res_id_list, chain_id = None)

        :param res_id_list: list of residue ids
        :param chain_id: optionally, if ``None``, it uses the all available chains

        :returns: a list with the names of the residue ids in the given list.

    .. method:: get_residues(chain_id = None, res_id_list = None)

        This method returns a list with residue objects of the residue ids in
        the given list, if ``None`` is given, it returns all residues.

        :param chain_id: optionally, if ``None``, it uses the all available chains
        :param res_id_list: list of residue ids from which one wants the objects

        :returns: a list of residues objects matching the given criteria

    .. method:: get_atoms_by_id(atom_id_list)

        This function returns the Atoms,  by their corresponding number from
        the pdb-file.

        :param atom_id_list: list of atom id numbers

        :returns: a list of atom objects to the matching atom ids

    .. method:: get_atoms_close_to_reference(reference_point, max_radius, min_radius = 0, atoms = None)

        This function returns all atoms of this structure, which lie in the
        range [min_radius, max_radius] to the reference point.

        :param reference_point: Numpy array of the reference
        :param max_radius: maximal distance to include the atoms
        :param min_radius: minimal distance of the atoms to the the reference
        :param atoms: optional list of atoms, if ``None`` is given it uses all
            atoms from the protein

        :returns: a list of atom objects close to the given reference point

    .. method:: find_chain_contacts(chain1, chain2, max_distance)

        Finds Atoms of chain1 which are within max_distance of chain2

        :param chain1: first chain id
        :param chain2: second chain id
        :param max_distance: Maximal distance to include atoms in the calculation

        :returns: a list of atom objects from chain1 which distance to chain2 is
            smaller than max_distance

    .. method:: get_sequence()

        :returns: a dictionary which contains the chain and the sequence:
            'A' : 'RG...CC'
            'B' : 'PW...FV'
            Non standard amino acids will not be returned!!!

    .. method:: snap_vdw_to_box(box_mesh_size, box_dim, box_offset, warning = True, vdw_radii = 2, increase_vdw_by = 0)

        This method snaps a structure to a given dxbox. If the structure is
        a pqr it uses the supplied vdw radii otherwise it uses the variable
        'vdw_radii' for each atom.

        If any coordinate lies outside the box an error will be printed to the
        standard output.

        :param box_mesh_size: mesh size [m,m,m]
        :param box_dim: [x,y,z]
        :param box_offset: [x_o,y_o,z_o]
        :param warning: boolean, print an error, if the structure does not fit
            completely into the given box dimensions
        :param vdw_radii: if this is a pdb file there are no other radii
            available (in Angstroem)
        :param increase_vdw_by: can be used to blow up the radii of each atom
            (in Angstroem)

        :returns: Numpy array with 0's outside and 1's inside the protein.

    .. method:: get_rmsd_rotation_translation_from_superposition(pdb_to_rotate, atom_types = None)

        This method tries to fit the given pdb onto this one. The method
        returns a dictionary, which contains the 'rmsd', 'rotation' matrix
        and the 'translation'.
        The parameter 'atom_types' can be used to supply a list of atoms, which
        should be fit onto each other, if 'None' is supplied, it will try to
        fit all atoms onto each other.
        In case the number of atoms does not match, it raises an Error!

        I guess the units are Angstroem.

        :param pdb_to_rotate: Structure derivative object to superimpose onto
            this object.
        :param atom_types: if ``None``, it tries to fit all atoms, but it can also
            be used to fit specific types (e.g. ['CA'], ['CA','N'])

        :returns: a dictionary which contains the following keys:

            * 'rmsd' : root mean square deviation
            * 'rotation' : rotation matrix
            * 'translation' : translation vector

    .. method:: superimpose_given_pdb_onto_self(pdb_to_superimpose, atom_types = None)

        This method superimposes this structure onto the given structure!
        If the number of atoms differ it raises an error!

        :param pdb_to_superimpose: structure derivative object to superimpose onto
            this object
        :param atom_types: if ``None``, it tries to fit all atoms, but it can also
            be used to fit specific types (e.g. ['CA'], ['CA','N'])

        :returns: ``None``

    .. method:: superimpose_self_onto_given_pdb(pdb_to_superimpose, atom_types = None)

        This method superimposes this structure onto the given structure!
        If the number of atoms differ it raises an error!

        :param pdb_to_superimpose: structure derivative object to superimpose
            this object onto
        :param atom_types: if ``None``, it tries to fit all atoms, but it can also
            be used to fit specific types (e.g. ['CA'], ['CA','N'])

        :returns: ``None``

    .. method:: get_dxbox_dim(box_mesh_size, extend = None, cubic_box = True, nlev = 4)

        This method returns the dimensions of a dxbox. The calculation is
        copied from the InFile class. If there have been changes this result
        might be wrong!
        The center of the box is the geometric center of the protein if not
        otherwise specified.

    .. method:: get_dxbox_offset(box_mesh_size, box_dim, box_center)

        Returns the offset for the given dimensions.

        :param box_mesh_size: [m,m,m]
        :param box_dim: [x,y,z]
        :param box_center: [x_c,y_c,z_c]

        :returns: a list box_offset: [x_o,y_o,z_o].

    .. method:: get_hydrophobic_potential(box_mesh_size, box_dim, box_offset)

        Calculate the hydrophobic potential.
        This method uses a simplified model for the hydrophobic potential.
        The 'charges' are taken from the Kyte and Doolittle hydrophobicity
        scale. For each residue the center is calculated and the potential
        is modelled as:

            :math:`\phi = \sum \left ( hydrophobic\_charge \times e^{( - distance )} \right )

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

        :param other_structure_object: structure which might have an overlapp with
          this one
        :param energy_cutoff: cutoff for the lennard jones potential to decide
          if there is an overlapp or not
        :param vdw_radii: if this is pdb object there are no radii information
          available

        :returns: counts of atoms in this structure that overlapp with
          the given structure object
        :rtype: int

    .. method:: get_contact_list(cutoff = 5.)

        :param cutoff: if any distance between atoms from two residues is less
            than the cutoff, it is a contact

        :returns: a list: [[...], [i,j,0], [...]] where 0 means no contact.

    .. function:: residue_contact(res_i, res_j, cutoff)

        Subroutine of :func:`Structure_Template.get_contact_list`.

.. class:: PDBFile(Structure_Template)

    Docstring missing.


    .. method:: get_pqr_structure(new_pqr_path = None, force_field = "amber", pdb2pqr_argv = None, pdb2pqr_path = "pdb2pqr", add_ions = True, add_chain = True)

        Call pdb2pqr to replace the b-factor and the occupancy information in
        the pdb file with the charges and the vdw radii. If this pdb Object has
        no structure_path property (i.e. it is None), then an error is raised.

        If the pdb contains CA, ZN or SO4 ions, they will be added, if not
        stated otherwise.

        :param new_pqr_path: path for the new pqr file. If ``None`` is given, it
            replaces the ".pdb" with ".pqr" at the end
        :param force_field: forcefield from which charges and radii should be
            taken. Default is amber
        :param pdb2pqr_argv: can contain additional arguments to pdb2pqr as a
            list (e.g. ['--assign-only'], oder ['--noopt']). If multiple
            additional arguments are given, they also have to be given as
            a list (e.g. ['--assign-only', '--noopt'])
        :param pdb2pqr_path: path to the executable of PDB2PQR
        :param add_ions: add CA, ZN, SO4 ions if they are in the pdb. The CA,ZN
            and SO4 atoms should have a residue name that fits their type
            (CA, ZN, SO4).

        :returns: a PQRFile object

    .. method:: _read_file()

        This method reads the data from the given file.

    .. method:: save_to_file(path)

        This method writes the structure to a given path.
        If the structure is a lattice, it will try to calculate the
        correct connections. The Assumption of a lattice is made on the
        number of residues and the number of overall atom coordinates.
        If both numbers are equal, it is very probable that this is a lattice
        protein. This only works, if there is just one chain!

        Notice: We do not work with multiple models!

    .. method:: _save_to_file(f)

        Docstring missing.

    .. method:: _get_atom_line(chain, res, atom)

        Formats the line for writing.

.. class:: PQRFile(Structure_Template)

    Docstring missing.


    .. method:: _read_structure()

        This method reads the data from the given file.

    .. method:: save_to_file(path)

        This method writes the structure to a given path.

        Notice: We do not work with multiple models!

    .. method:: _save_to_file(,f)

        Docstring missing.

    .. method:: _get_atom_line(chain, res, atom)

        Formats the line for writing.

    .. method:: snap_esp_to_dxbox(dxbox, warning = True)

        This method snaps the charge for each atom of this pqr structure to a
        given dxbox. The new array contains the charge in the unit of Coulomb.

        If any coordinate lies outside the box an error will be printed to the
        standard output.

        Be carefull!!! If you use a neutral probe and you mesh size is to
        large the dipole effect is not visible!

        :param dxbox: a dxbox object
        :param warning: either ``True`` or ``False``

        :returns: a numpy array with the charges in units of e at the center
          of each atom

.. class:: latFile(Structure_Template)

    Docstring missing.


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

    .. method:: _get_atom_line(chain, res, atom)

        Formats the line for writing.

    .. method:: _convert_fcc_fold(fold, mesh_size)

        Docstring missing.

    .. method:: _convert_sc_fold(fold, mesh_size)

        Docstring missing.

.. class:: entity(object)

    This class is more or less copied from biopython.

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

