
:mod:`Analysis` --- Hypervolume and Pareto analysis tools
=========================================================

.. module:: Analysis
   :synopsis: Hypervolume and Pareto analysis tools.
.. moduleauthor:: Christoph Wilms <christoph.wilms@uni-due.de>
.. moduleauthor:: Aaron Garrett <aaron.lee.garrett@gmail.com>
.. sectionauthor:: Jean-Noel Grad <jean-noel.grad@uni-due.de>


This module provides tools for the calculation of Hypervolume by Slicing
Objectives (HSO) of sets and for the Pareto analysis os sets. Part of this
module belongs to the package **inspyred 1.0** (:math:`\copyright` 2012
Inspired Intelligence Initiative), which can be found at `inspyred`_.
Following functions are part or derivatives of the original package:
:class:`hypervolume_max`, :class:`hypervolume_min`. They are therefore 
covered by the GNU General Public License version 3.0 (GPLv3), available
online at `opensource.org`_.

.. _inspyred: https://pypi.python.org/pypi/inspyred
.. _opensource.org: http://opensource.org/licenses/gpl-3.0.html


.. _Analysis-syntax:

Module Syntax
-------------

Empty.

.. _contents-of-module-Analysis:

Module Contents
---------------

.. function:: hypervolume(pareto_set, optimization_type, reference_point = None)

    Calculates the hypervolume by slicing objectives (HSO).

    :param pareto_set: the list or lists of objective values comprising the Pareto front
    :param optimization_type: 'min' for minimization or 'max' for maximization
    :param reference_point: the reference point to be used (default None)

.. function:: hypervolume_max(pareto_set, reference_point = None)

    Calculates the hypervolume by slicing objectives (HSO).
    
    This function calculates the hypervolume (or S-measure) of a nondominated
    set using the Hypervolume by Slicing Objectives (HSO) procedure of While, 
    *et al.*.

    .. seealso::

        While, *et al.*, *Heuristics for Optimising the Calculation of
        Hypervolume for Multi-objective Optimisation Problems*, IEEE CEC 2005
        (`original paper <http://www.lania.mx/~ccoello/EMOO/while05a.pdf.gz>`_)

    :param pareto_set: the list or lists of objective values comprising the
       Pareto front
    :param reference_point: the reference point may be specified or it may be
        left as the default value ``None``. In that case, the reference point
        is calculated to be the maximum value in the set for all objectives
        (the ideal point). This function  assumes that objectives are to be
        maximized.

    .. function:: dominates(p, q, k = None)

        Subroutine of :func:`hypervolume_max`.

    .. function:: insert(p, k, pl)

        Subroutine of :func:`hypervolume_max`.

    .. function:: slice(pl, k, ref)

        Subroutine of :func:`hypervolume_max`.

.. function:: hypervolume_min(pareto_set, reference_point = None)

    See :func:`hypervolume_max`.

    .. function:: dominates(p, q, k = None)

        Subroutine of :func:`hypervolume_min`.

    .. function:: insert(p, k, pl)

        Subroutine of :func:`hypervolume_min`.

    .. function:: slice(pl, k, ref)

        Subroutine of :func:`hypervolume_min`.

