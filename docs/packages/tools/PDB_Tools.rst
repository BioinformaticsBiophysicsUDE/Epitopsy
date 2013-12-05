
:mod:`PDB_Tools` --- YYYYYYY
======================================================

.. module:: PDB_Tools
   :synopsis: YYYYYYYY.
.. moduleauthor:: Christoph Wilms <christoph.wilms@uni-due.de>
.. sectionauthor:: Jean-Noel Grad <jean-noel.grad@uni-due.de>


This module provides YYYYY.


.. _PDB_Tools-syntax:

Module Syntax
-------------

Empty.

.. _contents-of-module-PDB_Tools:

Module Contents
---------------

.. function:: extract_pdbs_from_pdb_models(pdb_model_ensemble_path, new_dir, name_template, ref_pdb_path = None)

    This function reads all models from 'pdb_model_ensemble_path' and writes
    each model as a new pdb to the supplied directory. The new name of the 
    pdbs is: 'given_name' + '_x.pdb' where x is the model id.
    If a reference pdb is supplied, the function will fit the 'CA' atoms.


.. function:: get_one_atom_pdb(pdb_path, atom_coord)

    Write a PDB file which consists of only one C-alpha atom.
    
    :param pdb_path: path of the PDB file to create.
    :type pdb_path: str
    :param atom_coord: coordinates of the C-alpha atom.
    
    :returns: a PDB object.

