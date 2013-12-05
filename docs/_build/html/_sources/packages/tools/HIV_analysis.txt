
:mod:`HIV_analysis` --- HIV analysis operations
===============================================

.. module:: HIV_analysis
   :synopsis: YYYYYYYY.
.. moduleauthor:: Christoph Wilms <christoph.wilms@uni-due.de>
.. sectionauthor:: Jean-Noel Grad <jean-noel.grad@uni-due.de>


This module provides HIV analysis operations.


.. _HIV_analysis-syntax:

Module Syntax
-------------

Empty.

.. _contents-of-module-HIV_analysis:

Module Contents
---------------

.. function:: get_HXB2_seq_map(ref_seq)

    :param ref_seq: HXB2 sequence with gaps

    :returns: a dictionary mapping the di columns (starting with 1, ..., N) to
      the HXB2 sequence (starting with 1, ..., N).
    
.. function:: map_dca(dca_list, seq_map, seq, seq_dist)

    :param dca_list: [i,j,di], sorted with respect to di
    :param seq_map: maps the di columns (1,...,N) to the sequence positions
      (1, ..., N)
    :param seq: sequence without gaps
    :param seq_dist: minimum distance seperating two residues to be relevant

    :returns: a list of list, containing the di columns (1,...,N) which refer
      to HXB2, the sequence positions (1,...,N) itself, the domain
      information and the di value.
    
.. function:: reduce_dca(dca_list, n_dis, include_neighbors)

    :param dca_list: [\*,di,seq_neighbor], sorted with respect to di
    :param n_dis: number of di values
    :param include_neighbors: include neighbors or not (True/False)

    :returns: a list.
    
.. function:: filter_dca(dca_list, seq_map, seq, seq_dist, n_dis)

    :param dca_list: [i,j,di], sorted with respect to di
    :param seq_map: maps the di columns (1,...,N) to the sequence positions
      (1, ..., N)
    :param seq: sequence without gaps
    :param seq_dist: minimum  distance seperating two residues to be relevant
    :param n_dis: number of di values to check

    :returns: a list of list, containing the di columns (1,...,N) which refer
      to HXB2, the sequence positions (1,...,N) itself, the domain
      information and the di value.
    
.. function:: map_domain(x)

    :param x: position in the HXB2 sequence (1,...,N)

.. function:: get_color_for_domain(domain)

    Docstring missing.

