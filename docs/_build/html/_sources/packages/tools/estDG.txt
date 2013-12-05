
:mod:`estDG` --- YYYYYYY
======================================================

.. module:: estDG
   :synopsis: YYYYYYYY.
.. moduleauthor:: Christoph Wilms <christoph.wilms@uni-due.de>
.. sectionauthor:: Jean-Noel Grad <jean-noel.grad@uni-due.de>


This module provides YYYYY.


.. _estDG-syntax:

Module Syntax
-------------

Empty.

.. _contents-of-module-estDG:

Module Contents
---------------

.. function:: findVolLASPro(counter_box)

    Docstring missing.

.. function:: calcDG(energy_box,counter_box,Temp=310,conc=None)

    Calculates the binding free energy of protein to ligand in kJ/mol and 
    the k_D.

    :param energy_box: DXBox containing the energies
    :param counter_box: DXBox containing the number of possible rotations at
      each grid point
    :param Temp: temperature, default 310
    :param conc: concentration

    :returns: a list containing the binding free energy of protein
      to ligand and the k_D.
    
.. function:: calcConc(void,LASpoints,ProteinPoints,grid=[0.5,0.5,0.5])

    Calculate the concentration of protein and ligand for a given LAS layer
    and the additional volume of the protein and a given volume which
    constitutes the rest of the water in which all is soluted.

    All has to be given in grid points, the grid defines afterwards the volume
    of each gridpoint.


.. function:: calVoidVol(conc,LASpoints,ProteinPoints,grid=[0.5,0.5,0.5])

    Calculate a void volume, from a given LAS-layer, a protein volume (both in
    grid points) and a concentration in which the protein should be soluted.

