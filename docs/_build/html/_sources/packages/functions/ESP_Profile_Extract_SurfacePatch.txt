
:mod:`ESP_Profile_Extract_SurfacePatch` --- YYYYYYY
======================================================

.. module:: ESP_Profile_Extract_SurfacePatch
   :synopsis: YYYYYYYY.
.. moduleauthor:: Christoph Wilms <christoph.wilms@uni-due.de>
.. sectionauthor:: Jean-Noel Grad <jean-noel.grad@uni-due.de>


This module provides YYYYY.


.. _ESP_Profile_Extract_SurfacePatch-syntax:

Module Syntax
-------------

Usage::

    python ESP_Profile_Extract_SurfacePatch.py --pdb=<path> --coord=<[x,y,z]> --sphere-radius=<number>
    
Additional options:
    * --dx=True, if this if True, it looks for the files
      '*_vdw-PE0.dx' and '*_esp-PE0.dx', where '*' is the pdb_name.
    * --meshsize=<meshsize_in_x_y_z>

.. _contents-of-module-ESP_Profile_Extract_SurfacePatch:

Module Contents
---------------

.. class:: ESP_Profile_Extract_SurfacePatch(object)

    This class extracts the surface patch, which lies inside the volume of a 
    sphere around the given coordinates, in form a '.dx' file, which then
    contains the esp-values on the vdw surface.
    
    The extracted '.dx' file will be saved with the extension: 
    '*_extract_sp.dx', where '*' stands for the pdb_name without the'.pdb'.

