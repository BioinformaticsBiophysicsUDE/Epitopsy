
:mod:`ScoreHull` --- YYYYYYY
======================================================

.. module:: ScoreHull
   :synopsis: YYYYYYYY.
.. moduleauthor:: Christoph Wilms <christoph.wilms@uni-due.de>
.. sectionauthor:: Jean-Noel Grad <jean-noel.grad@uni-due.de>


This module provides YYYYY.


.. _ScoreHull-syntax:

Module Syntax
-------------

Usage::

    python ScoreHull.py --hullpdb=<path> --hull=<path> --target=<path> 

Optional arguments:
    * --meshsize=<number>, default=1
    * --out=<path>, default=no
    * --atom-coord=yes, default=no
    * --clean-up=yes, default=no 
    * --score=diff-vec, default=diff-vec
    * ?abs-vec
    * ?=no

.. _contents-of-module-ScoreHull:

Module Contents
---------------

.. class:: ScoreHull(object)

    This class scores the given profile for the given pdb structure.
    To compare structures it is necessary to provide the target structure, 
    so that the infile with the apbs parameters can be constructed from 
    this structure, so that the offset is the same. Furthermore it is
    necessary to provide the '--extend=<number>' parameter, that has been
    used with the target structure. If none has been used, there is no need
    to use it here either!
    If the option '--atom-coord=yes' is given, the hull is saved as 
    atomic coordinates instead of the array indices.
    The hull of the current protein can be saved if neccessary by
    '--out=<path>'.

    .. method:: score_diff_vec(vec1, vec2): # works best diff_vec = vec1 - vec2 score_diff = np.sqrt(np.dot(diff_vec, diff_vec)) if self.print_result is True

        Docstring missing.

    .. method:: score_abs_vec(vec1, vec2)

        Docstring missing.

    .. method:: get_score()

        Docstring missing.

