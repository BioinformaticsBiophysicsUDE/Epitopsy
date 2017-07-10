
:mod:`analysis` --- Analysis of energy grids
============================================

.. module:: analysis
   :synopsis: Analysis of energy grids.
.. moduleauthor:: Christoph Wilms <christoph.wilms@uni-due.de>
.. sectionauthor:: Jean-Noel Grad <jean-noel.grad@uni-due.de>


This module is used by
:func:`epitopsy.EnergyGrid.calculation.scan`
to determine the favorable (resp. unfavorable) regions of space around a
protein with regards to the ligand used as probe. Similarly, favorable
(resp. unfavorable) regions of the protein surface are determined. The results
of this calculation can be used to predict the results of an Alanine scan.


.. include:: <mmlalias.txt>
.. include:: <isogrk3.txt>
.. include:: <isotech.txt>
.. |_| unicode:: 0xA0 
   :trim:
.. |kbT| replace:: k\ :sub:`B`\ T


The energies returned by :func:`calc_vol_and_surface` are calculated using
the following formula (with *t* the **energy_cutoff**):

    +----------------------------+---------------------------------+---------------------------------------------------------------+
    |            Key             |          Value (formula)        |                        Description                            |
    +============================+=================================+===============================================================+
    | normal_total_volume        | :math:`V = V_+ + V_- + V_0`     | Total number of grid points outside the protein (int)         |
    |                            |                                 |                                                               |
    +----------------------------+---------------------------------+---------------------------------------------------------------+
    | normal_neutral_volume      | :math:`V_0 = \sum_i^V         \ | Number of neither favorable or unfavorable grid               |
    |                            | [+t \leq \Delta G_i \leq -t]`   | points outside the protein (int)                              |
    +----------------------------+---------------------------------+---------------------------------------------------------------+
    | normal_fav_volume          | :math:`V_+ = \sum_i^V         \ | Number of favorable grid points outside the protein (int)     |
    |                            | [\Delta G_i < -t]`              |                                                               |
    +----------------------------+---------------------------------+---------------------------------------------------------------+
    | normal_unfav_volume        | :math:`V_- = \sum_i^V         \ | Number of unfavorable grid points outside the protein (int)   |
    |                            | [\Delta G_i > +t]`              |                                                               |
    +----------------------------+---------------------------------+---------------------------------------------------------------+
    | normal_fav_volume_score    | :math:`wV_+ = \sum_i^V        \ | Summation of the Gibbs free energy over all favorable grid    |
    |                            | \Delta G_i [\Delta G_i < -t]`   | points outside the protein (float, < 0, in units of |kbT|)    |
    +----------------------------+---------------------------------+---------------------------------------------------------------+
    | normal_unfav_volume_score  | :math:`wV_- = \sum_i^V        \ | Summation of the Gibbs free energy over all unfavorable grid  |
    |                            | \Delta G_i [\Delta G_i > +t]`   | point outside the protein (float, > 0, in units of |kbT|)     |
    +----------------------------+---------------------------------+---------------------------------------------------------------+
    | normal_total_surface       | :math:`S = A_+ + A_- + V_0`     | Total number of grid points on the protein surface (int)      |
    |                            |                                 |                                                               |
    +----------------------------+---------------------------------+---------------------------------------------------------------+
    | normal_neutral_surface     | :math:`A_0 = \sum_i^S         \ | Number of neither favorable or unfavorable grid points on     |
    |                            | [+t \leq \Delta G_i \leq -t]`   | the protein surface (int)                                     |
    +----------------------------+---------------------------------+---------------------------------------------------------------+
    | normal_fav_surface         | :math:`A_+ = \sum_i^S         \ | Number of favorable grid points on the protein surface (int)  |
    |                            | [\Delta G_i < -t]`              |                                                               |
    +----------------------------+---------------------------------+---------------------------------------------------------------+
    | normal_unfav_surface       | :math:`A_- = \sum_i^S         \ | Number of unfavorable grid points on the protein surface      |
    |                            | [\Delta G_i > +t]`              | (int)                                                         |
    +----------------------------+---------------------------------+---------------------------------------------------------------+
    | normal_fav_surface_score   | :math:`wA_+ = \sum_i^S        \ | Summation of the Gibbs free energy over all favorable grid    |
    |                            | \Delta G_i [\Delta G_i < -t]`   | points on the protein surface (float, < 0, in units of |kbT|) |
    +----------------------------+---------------------------------+---------------------------------------------------------------+
    | normal_unfav_surface_score | :math:`wA_- = \sum_i^S        \ | Summation of the Gibbs free energy over all unfavorable grid  |
    |                            | \Delta G_i [\Delta G_i > +t]`   | points on the protein surface (float, > 0, in units of |kbT|) |
    +----------------------------+---------------------------------+---------------------------------------------------------------+

    Using the Iverson bracket notation: :math:`[P] = \begin{cases} 1 & \text{when }
    P \text{ is true} \\ 0 & \text{when } P \text{ is false} \end{cases}`


    Output from :file:`result_data.txt`:

    .. code-block:: none

        # energy cutoff:                      1.0
        DG:                        -7.41489284032
        conc:                                 1.0
        normal_total_volume:             16156569
        normal_neutral_volume:           15469239
        normal_fav_volume:                 566194
        normal_unfav_volume:               121136
        normal_fav_volume_score:   -989039.489084
        normal_unfav_volume_score:  177463.490532
        normal_total_surface:               98262
        normal_neutral_surface:             86543
        normal_fav_surface:                 11719
        normal_unfav_surface:                   0
        normal_fav_surface_score:  -20950.6659235
        normal_unfav_surface_score:           0.0


.. _EnergyGrid_analysis-syntax:

Module Syntax
-------------

Empty.

.. _contents-of-module-EnergyGrid_analysis:

Module Contents
---------------

.. automodule:: epitopsy.EnergyGrid.analysis
    :members:
    :undoc-members:
    :show-inheritance:

