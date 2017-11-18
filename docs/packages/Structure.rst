
:mod:`Structure` --- Additional classes for PDB files
=====================================================

.. module:: Structure
   :synopsis: Additional classes for PDB files.
.. moduleauthor:: Christoph Wilms <christoph.wilms@uni-due.de>
.. moduleauthor:: Thomas Hamelryck <thamelry@binf.ku.dk>
.. sectionauthor:: Jean-Noel Grad <jean-noel.grad@uni-due.de>


This module provides additional classes for the handling of PDB and PQR files.
Following classes are part of or derivatives of the original Biopython software
(`homepage <http://biopython.org/wiki/Biopython>`_,
`documentation <http://biopython.org/wiki/Documentation>`_,
`download <https://github.com/biopython/biopython>`_,
`license <../_static/licenses/biopython.txt>`_):
:class:`entity`, :class:`Structure`, :class:`Model`,
:class:`Chain`, :class:`Residue` and :class:`atom`.
They are still covered by their original Biopython license in Epitopsy.

.. _Structure-syntax:

Module Syntax
-------------

PDB files should be loaded this way:

    >>> a = Structure.PDBFile('4N6W.pdb')
    Warning: This structure contains disorderd atoms! Using only A-Positions!
    >>> # further operations...

.. _contents-of-module-Structure:

Module Contents
---------------

.. automodule:: epitopsy.Structure
    :members:
    :undoc-members:
    :show-inheritance:

