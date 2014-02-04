
:mod:`StartMD` --- Gromacs MD simulations
=========================================

.. module:: StartMD
   :synopsis: Gromacs MD simulations.
.. moduleauthor:: Christoph Wilms <christoph.wilms@uni-due.de>
.. sectionauthor:: Jean-Noel Grad <jean-noel.grad@uni-due.de>


This module provides an environment to start Gromacs MD simulations.

.. highlight:: bash

.. _StartMD-syntax:

Module Syntax
-------------

This script is called from the shell::

    python StartMD <pdb_path> <simulation_time> <temperature> <new_dir_template> <seeds>

* **pdb_path** refers to the PDB file,
* **simulation_time** is the simulation duration in nanoseconds,
* **temperature** is the temperature in Kelvins,
* **new_dir_template** is a directory (what for?),
* **seeds** are the MD seeds.


.. highlight:: python

.. _contents-of-module-StartMD:

Module Contents
---------------

There is no function of class in this module.
