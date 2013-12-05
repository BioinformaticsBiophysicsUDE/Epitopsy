
:mod:`calculate_partition_function` --- YYYYYYY
======================================================

.. module:: calculate_partition_function
   :synopsis: YYYYYYYY.
.. moduleauthor:: Christoph Wilms <christoph.wilms@uni-due.de>
.. sectionauthor:: Jean-Noel Grad <jean-noel.grad@uni-due.de>


This module provides YYYYY.


.. _calculate_partition_function-syntax:

Module Syntax
-------------

Empty.

.. _contents-of-module-calculate_partition_function:

Module Contents
---------------

.. function:: calculate_interaction_energy(pdb_path, ligand_path, mesh_size, number_of_rotations=150, Temperature=310., extend=None, ph=None, use_pdb2pqr=True, cubic_box=True, center_pdb=False, box_center=[0, 0, 0], box_dim=None, box_type=["esp", "vdw"], write_top_hits=False, explicit_sampling=False, zipped=False)

    This function calculates the difference in Gibbs free energy in units
    of :math:`k_bT`.

    Probability of all states at position :math:`\vec{r}`:

        :math:`p(r) = \frac{1}{Z} \sum_i^N e^{\left ( \frac{-E_i}{k_b T} \right )}`

    Probability of all states at position :math:`\vec{r}` in solution (with :math:`E_i = 0` for all rotations):

        :math:`p(r_{solv}) = \frac{1}{Z} \sum_i^N e^{\left ( \frac{-E_i}{k_b T} \right )} = \frac{N}{Z}`

    Dissociation constant:

        :math:`K = \frac{p(r)}{p(r_{solv})} = \frac{1}{N} \sum_i^N e^{\left ( \frac{-E_i}{k_b T} \right )}`

    Gibbs free energy:

        :math:`\Delta G = - k_bT \times ln(K)`

    :param pdb_path: path to the fixed pdb file
    :param ligand_path: path to the .pqr file of the ligand. the .pdb file has
        to be in the same folder
    :param mesh_size: grid mesh size in all three dimensions [mx,my,mz]
    :param number_of_rotations: how many different orientations for the ligand?
    :param Temperature: temperature for the electrostatic potential calculations
    :param extend: extend the box dimensions by this factor
    :param ph: can be used to calculate the pqr for a specific pH, default is
        ``None``
    :param use_pdb2pqr: use pdb2pqr to calculate the pqr of pdb_path, if
        ``False`` a pqr with the same name, but ending with '.pqr' has
        to be in the same folder
    :param cubic_box: cubic box or not, hence ``True`` or ``False``
    :param center_pdb: if ``True`` the pdb will be centered at box_center,
        default is ``False``
    :param box_center: center of the apbs calculations, default is [0,0,0]
    :param box_dim: dimension of the box, default is ``None`` and it is
        calculated automatically
    :param box_type: types of boxes APBS should write to disk. We need at least
        ["esp","vdw"], "smol" is optional :math:`\rightarrow` ["esp","vdw","smol"]
    :param write_top_hits: write the top scoring conformations
    :param explicit_sampling: if ``True`` it does not use FFT correlation for
        the detection of overlapping orienations. This option does
        not work with write_top_hits. Default is ``False``
    :param zipped: if ``True`` all dxboxes will be gzipped, else the energy
        matrix and the energy matrix with the LAS surface
        are not zipped. Default is ``False``

    :returns: the energy

