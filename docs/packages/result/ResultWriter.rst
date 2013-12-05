
:mod:`ResultWriter` --- YYYYYYY
======================================================

.. module:: ResultWriter
   :synopsis: YYYYYYYY.
.. moduleauthor:: Christoph Wilms <christoph.wilms@uni-due.de>
.. sectionauthor:: Jean-Noel Grad <jean-noel.grad@uni-due.de>


This module provides YYYYY.


.. _ResultWriter-syntax:

Module Syntax
-------------

.. warning:: Imports at the bottom ... circular imports ... what a mess!

.. _contents-of-module-ResultWriter:

Module Contents
---------------

.. class:: ResultWriter

    Again I copied Niko's class, but I had to change the function print as this 
    is a keyword in python:

        print(...) :math:`\rightarrow` printResult(...)

    .. method:: write(filename, rs, append = False)

        Writes a given ResultSet to specified file, overwriting old results.

    .. method:: printResult(r)

        Prints Result as formatted String to given PrintWriter Object.

    .. method:: append(filename, rs)

        Docstring missing.

    .. method:: createResultFile(path, psiAxis = False)

        Docstring missing.

