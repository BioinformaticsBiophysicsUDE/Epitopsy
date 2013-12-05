
:mod:`CalculateESPProfiles` --- YYYYYYY
======================================================

.. module:: CalculateESPProfiles
   :synopsis: YYYYYYYY.
.. moduleauthor:: Christoph Wilms <christoph.wilms@uni-due.de>
.. sectionauthor:: Jean-Noel Grad <jean-noel.grad@uni-due.de>


This module provides YYYYY.


.. _CalculateESPProfiles-syntax:

Module Syntax
-------------

Usage::

    python CalculateESPProfiles --list=<path> --out=<path>

Optional arguments:

    * --extend=

        * biggest :math:`\rightarrow` extends the hull of
          the biggest structure until all structure
          fit into this hull
        * all :math:`\rightarrow` default! Finds a hull
          that covers all  structures and extends it
          by one mesh unit

    * --hull-coordinates=<path>
    * --meshsize=<number>
    * --clean-up=no

.. _contents-of-module-CalculateESPProfiles:

Module Contents
---------------

.. class:: CalculateESPProfiles

    This class reads pdb-names from a given list and calculates the 
    electrostatic potential for each structure. Then it extends the surface
    until all structures are included. Therefore the structures should be 
    already aligned!
    Once the optimal shell is found it extracts the electrostatic potential
    at these points for all structures and writes them to a file.
    The template for the optimal shell is the largest structure. 

    With the option '--no-clean-up' all '.pqr' and  '.in' files are kept, 
    otherwise they will be deleted.

    .. method:: read_pdb_list(pdb_list_name)

        Reads pdb-filenames from a given list and loads them as pdb structures.

    .. method:: extend_all(big_box, template_box, extend_by = 1)

        This method works like this:
        `template_box` is 0 on all points, where a protein is and so all
        points are already inside the 'hull'. Then the hull is extended by a
        given unit and thats it.
        The hull sits much more closely around all structures.

    .. method:: extend_biggest(big_box, template_box, pdb_list, max_diameter_id)

        Flip `template_box` 1 :math:`\rightarrow` 0, 0 :math:`\rightarrow` 1, then
        multiplication with the  biggest hull should give all 0, if not, then the
        hull needs to be extended. 

