
:mod:`run_dca` --- Subroutines for Direct Coupling Analysis
===========================================================

.. module:: run_dca
   :synopsis: Subroutines for Direct Coupling Analysis.
.. moduleauthor:: Christoph Wilms <christoph.wilms@uni-due.de>
.. sectionauthor:: Jean-Noel Grad <jean-noel.grad@uni-due.de>


This module provides subroutines for Direct Coupling Analysis.


.. _run_dca-syntax:

Module Syntax
-------------

This script is called from cmd::

    $> python run_dca.py algorithm=<algo> in=<input> out=<output> x=<threshold>

* **algo** can be either *dca_new*, *dca_old* or *mi*,
* **input** is a path to the aligned structure,
* **output** is the path for the output,
* **threshold** is the identity threshold.

.. _contents-of-module-run_dca:

Module Contents
---------------

.. function:: run_dca_old(aln_path, di_path, identity_threshold)

    Docstring missing.

.. function:: run_dca_new(aln_path, di_path, identity_threshold)

    Docstring missing.

.. function:: run_mi(aln_path, di_path)

    Docstring missing.

