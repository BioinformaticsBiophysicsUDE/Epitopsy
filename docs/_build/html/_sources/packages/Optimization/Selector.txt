
:mod:`Selector` --- YYYYYYY
======================================================

.. module:: Selector
   :synopsis: YYYYYYYY.
.. moduleauthor:: Christoph Wilms <christoph.wilms@uni-due.de>
.. sectionauthor:: Jean-Noel Grad <jean-noel.grad@uni-due.de>


This module provides YYYYY.


.. _Selector-syntax:

Module Syntax
-------------

Empty.

.. _contents-of-module-Selector:

Module Contents
---------------

.. function:: random_tournament(parents, children, pareto_front, optimization_type)

    Docstring missing.

.. function:: nsga_tournament(parents, children, pareto_front, optimization_type)

    Replaces population using the non-dominated sorting technique from
    NSGA-II. Fill new parent population according to the best front
    respectively crowding distance within a front

.. function:: nsga_random_tournament(parents, children, pareto_front, optimization_type)

    Docstring missing.

.. function:: pareto_energy_tournament(parents, children, pareto_front, optimization_type)

    Docstring missing.

.. function:: pool_pareto_distance_tournament(parents, children, pareto_front, optimization_type)

    Docstring missing.

.. function:: who_is_fitter_random(c1, c2, optimization_type)

    This method returns the fitter individual form the two candidates
    c1 and c2. In case none dominates the other, one is randomly chosen.

.. function:: who_is_fitter_max_distance_all(c1, c2, population, optimization_type)

    This method returns the fitter individual form the two candidates
    c1 and c2. In case none dominates the other, the distance to all other
    individual is calculated and the one with the greater distance is 
    chosen to increase the diversity.

.. function:: who_is_fitter_min_energy_all(c1, c2, population, optimization_type)

    This method returns the fitter individual form the two candidates
    c1 and c2. In case none dominates the other, the distance to all other
    individual is calculated and the one with the greater distance is 
    chosen to increase the diversity.

