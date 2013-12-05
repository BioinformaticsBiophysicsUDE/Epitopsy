
:mod:`CreateHull` --- YYYYYYY
======================================================

.. module:: CreateHull
   :synopsis: YYYYYYYY.
.. moduleauthor:: Christoph Wilms <christoph.wilms@uni-due.de>
.. sectionauthor:: Jean-Noel Grad <jean-noel.grad@uni-due.de>


This module provides YYYYY.


.. _CreateHull-syntax:

Module Syntax
-------------

Usage::

    python CreateHull.py --pdb=<path> --out=<path>

Optional arguments:

    * --extend=<number_in_A>, default=6 (Nikos suggestion)
    * --atom-coord=yes, default=no
    * --clean-up=yes, default=no

.. _contents-of-module-CreateHull:

Module Contents
---------------

.. function:: Create_Hull(pdb_name, out_file_name, mesh_size, extend = 6, clean_up = False)

    This class calculates the electrostatic hull for one given protein 
    structure and saves it under the given path. 
    The hull can be extended by the command '--extend=<number>' to account
    for structural changes if it is compared to some other profile.
    If the option '--atom-coord=yes' is given, the hull is saved with 
    atomic coordinates instead of the array indices.

.. class:: CreateHull2(object)

    This class calculates the electrostatic hull for one given protein 
    structure and saves it under the given path. 
    The hull can be extended by the command '--extend=<number>' to account
    for structural changes if it is compared to some other profile.
    If the option '--atom-coord=yes' is given, the hull is saved with 
    atomic coordinates instead of the array indices.

    .. method:: run_apbs(pdb, extend, mesh_size)

        Docstring missing.

