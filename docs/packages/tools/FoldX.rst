
:mod:`FoldX` --- FoldX
======================

.. module:: FoldX
   :synopsis: YYYYYYYY.
.. moduleauthor:: Christoph Wilms <christoph.wilms@uni-due.de>
.. sectionauthor:: Jean-Noel Grad <jean-noel.grad@uni-due.de>


This module provides YYYYY.


.. _FoldX-syntax:

Module Syntax
-------------

Empty.

.. _contents-of-module-FoldX:

Module Contents
---------------

.. class:: FoldX()

   The :class:`FoldX` class supports the following methods and attributes:

    .. attribute:: pdb_name

        refers to a structure that should be mutated and analyzed
        without '.pdb'!

    .. attribute:: sequence_parent

        the sequence of the pdb (pdb_name)

    .. attribute:: sequence_child

        a mutated version of sequence_parent

    .. attribute:: repair_flag

        indicates, if the pdb should be repaired at first

    .. attribute:: new_name

        the name of the scored structure

    .. attribute:: Temp

        the temperature of the experiment in K

    .. attribute:: pH

        the pH of the experiment

    .. method:: _make_FoldX_run_Repair()

        Docstring missing.

    .. method:: _make_FoldX_run_Build()

        Docstring missing.

    .. method:: _make_FoldX_mut_List()

        Docstring missing.

    .. method:: _make_FoldX_pdb_List(pdb_name)

        :param pdb_name: refers to a structure that should be mutated
           and analyzed without '.pdb'!

    .. method:: _extract_FoldX_ddG()

        Docstring missing.

    .. method:: extract_FoldX_ddG()

        Docstring missing.

    .. method:: run_FoldX()

        Docstring missing.

