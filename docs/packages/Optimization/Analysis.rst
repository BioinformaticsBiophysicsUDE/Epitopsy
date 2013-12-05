
:mod:`Analysis` --- YYYYYYY
======================================================

.. module:: Analysis
   :synopsis: YYYYYYYY.
.. moduleauthor:: Christoph Wilms <christoph.wilms@uni-due.de>
.. sectionauthor:: Jean-Noel Grad <jean-noel.grad@uni-due.de>


This module provides YYYYY.


.. _Analysis-syntax:

Module Syntax
-------------

Empty.

.. _contents-of-module-Analysis:

Module Contents
---------------

.. function:: hypervolume(pareto_set, optimization_type, reference_point = None)

    Calculates the hypervolume by slicing objectives (HSO).
    
    The calculation has been 'stolen' from inspyred, but it is open source
    and so is this code :-)

    :param pareto_set: the list or lists of objective values comprising the Pareto front
    :param optimization_type: 'min' for minimization or 'max' for maximization
    :param reference_point: the reference point to be used (default None)

.. function:: hypervolume_max(pareto_set, reference_point = None)

    Calculates the hypervolume by slicing objectives (HSO).
    
    This function calculates the hypervolume (or S-measure) of a nondominated
    set using the Hypervolume by Slicing Objectives (HSO) procedure of `While, et al. 
    (IEEE CEC 2005) <http://www.lania.mx/~ccoello/EMOO/while05a.pdf.gz>`_.
    The *pareto_set* should be a list of lists of objective values.
    The *reference_point* may be specified or it may be left as the default 
    value of None. In that case, the reference point is calculated to be the
    maximum value in the set for all objectives (the ideal point). This function 
    assumes that objectives are to be maximized.

    :param pareto_set: the list or lists of objective values comprising the Pareto front
    :param reference_point: the reference point may be specified or it may be left as
        the default value ``None``. In that case, the reference point is calculated
        to be the maximum value in the set for all objectives (the ideal point).
        This function  assumes that objectives are to be maximized.

.. function:: dominates(p, q, k = None)

    Docstring missing.

.. function:: insert(p, k, pl)

    Docstring missing.

.. function:: slice(pl, k, ref)

    Docstring missing.

.. function:: hypervolume_min(pareto_set, reference_point = None)

    
    Calculates the hypervolume by slicing objectives (HSO).
    
    This function calculates the hypervolume (or S-measure) of a nondominated
    set using the Hypervolume by Slicing Objectives (HSO) procedure of `While, et al. 
    (IEEE CEC 2005) <http://www.lania.mx/~ccoello/EMOO/while05a.pdf.gz>`_.
    The *pareto_set* should be a list of lists of objective values.
    The *reference_point* may be specified or it may be left as the default 
    value of None. In that case, the reference point is calculated to be the
    maximum value in the set for all objectives (the ideal point). This function 
    assumes that objectives are to be minimized.
    
    ### I have hijacked the function and now it calculates the correct hypervolume 
    for a minimization! ###

    :param pareto_set: the list or lists of objective values comprising the Pareto front
    :param reference_point: the reference point may be specified or it may be left as
        the default value ``None``. In that case, the reference point is calculated
        to be the maximum value in the set for all objectives (the ideal point).
        This function  assumes that objectives are to be maximized.

.. function:: dominates(p, q, k = None)

    Docstring missing.

.. function:: insert(p, k, pl)

    Docstring missing.

.. function:: slice(pl, k, ref)

    Docstring missing.

