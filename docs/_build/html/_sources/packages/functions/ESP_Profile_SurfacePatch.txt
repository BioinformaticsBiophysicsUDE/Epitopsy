
:mod:`ESP_Profile_SurfacePatch` --- YYYYYYY
======================================================

.. module:: ESP_Profile_SurfacePatch
   :synopsis: YYYYYYYY.
.. moduleauthor:: Christoph Wilms <christoph.wilms@uni-due.de>
.. sectionauthor:: Jean-Noel Grad <jean-noel.grad@uni-due.de>


This module provides YYYYY.


.. _ESP_Profile_SurfacePatch-syntax:

Module Syntax
-------------

Usage:
    python ESP_Profile_SurfacePatch.py --pdb=<path> --out=<path> --coord=<[x,y,z] or [x  y  z]> --extend=<number> --sphere-radius=<number>

.. note::

    the best way to supply the '--coord' is to use 'str(list(coord))'
    in a python call (real space)
    
Additional options:

    * --meshsize=<meshsize_in_x_y_z>

.. _contents-of-module-ESP_Profile_SurfacePatch:

Module Contents
---------------

.. class:: ESP_Profile_SurfacePatch(object)

    This class calculates the ESP-Profie for a patch on the surface of the 
    given protein. The patch is calculated like this:

    #. the vdw surface of the protein is extended by the supplied value
    #. a sphere with the given value is created around the supplied
       coordinates
    #. all points on the protein surface which lie inside the sphere
       make up the surface patch
    
    This patch is then written to a file for further analysis.
    
    The input pdb should only contain the proteins, from which the esp hull
    is of interest. If one wants to know the esp-profile of a binding site
    for a ligand, that ligand has to be removed!!!
    
    Parameters of 'extend' and 'sphere radius' are in Angstroem.


