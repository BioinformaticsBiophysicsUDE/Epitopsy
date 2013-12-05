
:mod:`APBS` --- APBS automation
===============================

.. module:: APBS
   :synopsis: APBS automation.
.. moduleauthor:: Christoph Wilms <christoph.wilms@uni-due.de>
.. sectionauthor:: Jean-Noel Grad <jean-noel.grad@uni-due.de>


This module helps automatize APBS calculations.


.. _APBS-syntax:

Module Syntax
-------------

Empty.

.. _contents-of-module-APBS:

Module Contents
---------------

.. class:: APBSWrapper

    Docstring missing.

    .. method:: runPDB2PQR(pdb_path, pqr_path, force_field = 'amber', pdb2pqr_argv = None)

        Call pdb2pqr to replace the b-factor and the occupancy information
        in the pdb file with the charges and the vdw radii.

        :param pdb_path: path to the pdb file
        :param pqr_path: path for the new pqr file
        :param ph: pH at for which the charges should be calculated. Default is
            ``None``, which is close to a pH of 7
        :param force_field: forcefield from which charges and radii should be
            taken. Default is "amber"
        :param pdb2pqr_argv: Can contain additional arguments to pdb2pqr as a
            list (e.g. ['--assign-only'], oder ['--noopt']). If multiple
            additional arguments are given, they also have to be given as
            a list (e.g. ['--assign-only', '--noopt']).

        :returns: ``None``

    .. method:: runAPBS(apbsInParameters, inFilename)

        :param apbsInParameters: contains all neccessary parameters
        :param inFilename: name of the inputfile for apbs

        :returns: ``None``

    .. method:: get_dxbox(pdb_path, mesh_size, **kwds)

        This method starts pdb2pqr and apbs and returns the specified boxes.
        Be careful with absolute paths!

        :param pdb_path: path to the structure
        :param mesh_size: grid dimensions
        
        :param pqr_path: if ``None`` is supplied, it uses pdb2pqr to
            calculate one (optional)
        :param box_dim: if box dimensions and a box center are supplied, it
            will not use the pdb for the construction of the
            box. Furthermore it is not possible to extend a
            box, if the box properties are supplied (optional)
        :param box_center: determine the center of the box. If ``None`` is
            supplied, it calculates the center of the pqr file (optional)
        :param extend: Extend grid by the given number of Angstroms (optional)
        :param box_type: specify the returned objects: 'esp', 'vdw', 'smol' (optional)
        :param cubic_box: if the box is generated from the pqr file, it can
            be made to be cubic (optional)
        :param close_boundaries: if ``True`` it uses another algorithm for the
            boundary conditions of the box which
            yields better results but is slower (optional)
        :param temperatur: the temperature in Kelvin (optional)
        :param ph: the pH (optional)

        :returns: a dictionary with DXBox objects, which keys are the requested
            box types

    .. method:: get_binding_energy(complex_pqr_path, ligand_pqr_path, fixed_pqr_path, box_mesh_size, extend = None, **kwds)

        Binding energy is in kJ/mol.
        Obtain the change in total polar solvation energy by:

            :math:`dG\_bind = comp\_solv - fixed\_solv - ligand\_solv`

        It gives the total polar solvation energy which includes both the
        reaction field and coulombic contributions.
        Positive values are favorable.

        **kwds is used to set the options of the infile for the apbs calculation.

        :param complex_pqr_path: path of the complex of ligand and fixed protein
        :param ligand_pqr_path: path of the isolated ligand
        :param fixed_pqr_path: path of the isolated fixed protein
        :param box_mesh_size: specify the mesh size of the grid in Angstroem
        :param extend: increase the box dimensions

        :returns: binding energy
        :rtype: float

    .. method:: get_binding_energy_long(complex_pqr_path, ligand_pqr_path, fixed_pqr_path, box_mesh_size, extend = None, **kwds)

        Binding energy is in kJ/mol:

            :math:`solvation\_energy = (complex\_solv - complex\_ref) - (ligand\_solv - ligand\_ref) - (fixed\_solv - fixed\_ref)`

            :math:`coulomb\_energy = frac{complex\_coulomb - ligand\_coulomb - fixed\_coulomb}{pdie}`

            :math:`binding\_energy = solvation\_energy + coulomb\_energy`

        Positive values are favorable.

        **kwds is used to set the options of the infile for the apbs calculation.

        :param complex_pqr_path: path of the complex of ligand and fixed protein
        :param ligand_pqr_path: path of the isolated ligand
        :param fixed_pqr_path: path of the isolated fixed protein
        :param box_mesh_size: specify the mesh size of the grid in Angstroem
        :param extend: increase the box dimensions

        :returns: binding energy
        :rtype: float

    .. method:: get_dissociation_energy(complex_pqr_path, ligand_pqr_path, fixed_pqr_path, box_mesh_size, extend = None, **kwds)

        Dissociation energy is in kJ/mol.
        Obtain the change in total polar solvation energy by:

            :math:`dG\_diss = - \left ( comp\_solv - fixed\_solv - ligand\_solv \right )`

        It gives the total polar solvation energy which includes both the
        reaction field and coulombic contributions.

        **kwds is used to set the options of the infile for the apbs
        calculation.

        :param complex_pqr_path: path of the complex of ligand and fixed protein
        :param ligand_pqr_path: path of the isolated ligand
        :param fixed_pqr_path: path of the isolated fixed protein
        :param box_mesh_size: specify the mesh size of the grid in Angstroem
        :param extend: increase the box dimensions

        :returns: dissociation energy
        :rtype: float

.. function:: get_coulomb_energy(protein_pqr_path)

    :param protein_pqr_path: path to the pqr of the protein

    :returns: energy in kJ/mol

.. class:: InFile

    This class holds parameters to run APBS calculations.
    Only the most basic parameters are covered here,  intended for the
    generation of basic grids holding information of electrostatic
    potential and/or the van-der-Waals surface

    .. attribute:: pqr_path

        pqr path, from which the potential energy should be
        calculated or the path of the complex for which the binding
        energy should be calculated

    .. attribute:: calculation_type

        potential or binding_energy

    .. attribute:: box_mesh_size

        mesh size

    .. attribute:: box_dim

        special dimension can be supplied (although they might be
        fixed to work properly with apbs), if ``None`` is given it
        will be calculated from the size of the protein.

    .. attribute:: box_center

        center of the box, if ``None`` is given if will be set to the
        geometric center of the protein

    .. attribute:: extend

        can only be used when no box_dim is supplied, extends the box
        size by the given amount (in Angstroem)

    .. attribute:: cubic_box

        determine wheter it is a cubic box or not, in the case of
        a cubic box the largest dimension is used for all
        dimensions

    .. attribute:: box_type

        * esp : electrostatic potential
        * vdw : van der waals based solvent accessibility
        * smol : solvent accessibility
        * ndens : total mobile ion number density in units of M
        * qdens : total mobile charge density in units of e_c M

    .. attribute:: ligand_pqr_path

        needs to be supplied if the calculation is a
        binding energy calculation

    .. attribute:: fixed_pqr_path

        needs to be supplied if the calculation is a
        binding energy calculation

    .. method::  setGridCenter(center)

        Docstring missing.

    .. method::  setGridSize(size)

        Docstring missing.

    .. method::  setMeshSize(meshSize)

        Docstring missing.

    .. method::  generateFromPDB(pdb, padding, cubicBox)

        Docstring missing.

    .. method::  generateFromPDB2(pdb, padding, minDiameter, cubicBox)

        Docstring missing.

    .. method:: calculateValidDimension(c)

        Due to a multilevel approach APBS requires the grid to be of certain sizes.
        (See APBS manual for more information.)

        Self method ensures, that chosen grid dimensions meet these requirements.
        Current grid dimensions will be enlarged accordingly.

        :param c: test grid dimension

        :returns: integer number that has the correct dimension

    .. method:: fixGridSize()

        Docstring missing.

    .. method:: write(file_path)

        Docstring missing.

    .. method:: write_potential(file_path)

        This is the function that writes an infile for the calculation of
        potential grids.

        :param file_path: path to the new infile

        :returns: ``None``

    .. method:: write_binding_energy(file_path)

        This is the function that writes an infile for binding energy
        calculations.

        :param file_path: path to the new infile

        :returns: ``None``

    .. method:: write_binding_energy_long(file_path)

        This is the function that writes an infile for binding energy
        calculations.

        :param file_path: path to the new infile

        :returns: ``None``

    .. method:: checkInFile()

        Checks the information in the infile, java leftover ... pointless!

    .. method:: set_options(apbs_input_dict)

        This method accepts a dictionary and sets the elements from the
        dictionary (keys). If there is an element which could not be set, it
        will raise an error.
        The given values of the dictionary are not checked for validation.
        If some options are not given the default values will be used.

