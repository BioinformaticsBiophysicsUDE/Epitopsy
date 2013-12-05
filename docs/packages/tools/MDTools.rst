
:mod:`MD_Tools` --- Molecular Dynamics tools
============================================

.. module:: MD_Tools
   :synopsis: Molecular Dynamics tools.
.. moduleauthor:: Christoph Wilms <christoph.wilms@uni-due.de>
.. sectionauthor:: Jean-Noel Grad <jean-noel.grad@uni-due.de>


This module provides some tools for Molecular Dynamics simulations.


.. _MD_Tools-syntax:

Module Syntax
-------------

Empty.

.. _contents-of-module-MD_Tools:

Module Contents
---------------

.. function:: dump_gromacs_snapshots(trajectory_file, tpr_file, delta_t, new_pdb_path)

    This function dumps a pdb snapshot at the given time interval. The data
    will be stored in the new pdb path.


.. function:: calculate_rmsd_from_pdbmodel(pdb_model_path, pdb_frame_dir, ref_pdb_path, atom_types, delta_t)

    Docstring missing.

.. function:: calculate_rmsf_from_pdbmodel(pdb_model_path, pdb_frame_dir, fit_atom_types = None, rmsf_type = None)

    Docstring missing.

.. function:: calculate_center_of_mass(res)

    Docstring missing.

.. function:: calculate_rmsd_from_trajectory(traj_file, tpr_file, delta_t, ref_pdb_path, atom_types = None)

    This function returns two lists, which contain the time and the
    rmsd to the given reference structure path. 
    The function extracts snapshots from the trajectory file and stores them 
    in /dev/shm.


.. function:: calculate_rmsf_from_trajectory(traj_file, tpr_file, delta_t, fit_atom_types = None, rmsf_type = None)

    This function returns two lists, one with a list of indieces (e.g. res 
    numbers) and one with the calculated rmsfs.
    The unit of the rmsf is Angstrom.

    :param fit_atom_types: atom types that are used to fit the trajectory 
          onto the reference structure
    :param rmsf_type: 'res', 'CA', None; types that are used for the
          calculation of the rmsf. 'res' means rmsf per residue, 'CA' 
          calculates the rmsf for each c alpha atom and None for the
          geometric center of the protein.

