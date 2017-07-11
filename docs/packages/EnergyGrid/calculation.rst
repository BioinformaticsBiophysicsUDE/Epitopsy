
:mod:`calculation` --- Grid-based protein-ligand affinity screening
===================================================================

.. module:: calculation
   :synopsis: Grid-based protein-ligand affinity screening.
.. moduleauthor:: Christoph Wilms <christoph.wilms@uni-due.de>
.. sectionauthor:: Jean-Noel Grad <jean-noel.grad@uni-due.de>


This module provides a function to probe the surface of a protein with a
ligand and return an energy grid.

.. note::

    Following dependencies must be present on the system:

    * either anfft or FFTW3 (see :doc:`../../install/fftw`),
    * PyMOL v1.6.0.0 or superior (see :doc:`../../install/pymol`),
    * APBS (see :doc:`../../install/apbs`),
    * optionally, PDB2PQR (see :doc:`../../install/pdb2pqr`)


..  /========================================================================\
    |                 Insert unicode aliases in reST files                   |
    |========================================================================|
    |                                                                        |
    | Find a special character in charmap, take its code (example: with      |
    | "U+00C5", take "00C5") and grep it from the reST parser libraries:     |
    |                                                                        |
    |   $> (cd /usr/share/docutils/parsers/rst/include; grep "00C5" *)       |
    |   isolat1.txt      |Aring|  unicode: U+000C5 LETTER A WITH RING ABOVE  |
    |   xhtml1-lat1.txt  |Aring|  unicode: U+000C5 LETTER A WITH RING ABOVE  |
    |                                                                        |
    | Include "isolat1.txt" in the reST file:                                |
    |                                                                        |
    |   .. include:: <isolat1.txt>                                           |
    |                                                                        |
    | Insert the special character somewhere:                                |
    |                                                                        |
    |   Below 0.40 |Aring|, the calculation will slow down.                  |
    \________________________________________________________________________/


.. include:: <isolat1.txt>
.. include:: <mmlalias.txt>
.. include:: <isogrk1.txt>
.. |kbT| replace:: k\ :sub:`B`\ T 
.. |_| unicode:: 0xA0 
   :trim:

.. _calculation-syntax:

Module Syntax
-------------

Following files have to be gathered in the same directory before the
calculation starts:

* :file:`protein.pdb`
* :file:`ligand.pqr`

To start the Epitopsy calculation::

    >>> from epitopsy import EnergyGrid as EG
    >>> m = 0.80 # mesh size = grid resolution
    >>> EG.electrostatics('protein.pdb', 'ligand.pqr', mesh_size=[m,m,m],
    ...                   center_pdb=True, cubic_box=False, verbose=True)
    >>> EG.scan('protein.pdb', 'ligand.pqr', number_of_rotations=150)

Following files are generated in the same directory:

* :file:`protein.in` (APBS input file)
* :file:`io.mc` (APBS log file, showing the calculation progress)
* :file:`result_data.txt` (various variables on the calculation)
* :file:`protein_esp.dx` (electrostatic potential)
* :file:`protein_vdw.dx` (van der Waals surface)
* :file:`protein_epi.dx` (energy grid)
* :file:`protein_mic.dx.gz` (number of allowed rotations)

The most important file is :file:`protein_epi.dx`. It can be displayed
as isosurfaces in PyMOL::

    >>> reinitialize
    >>> load protein.pdb
    >>> load protein_epi.dx
    >>> bg_color white
    >>> as cartoon, protein
    >>> show spheres, protein and elem CA+ZN
    >>> isosurface negative, protein_epi, -1
    >>> isosurface positive, protein_epi, +1
    >>> color firebrick, negative
    >>> color skyblue, positive

Algorithm
---------

A fixed and rigib protein (**pdb_path**) is centered in a DXBox slightly
larger than its maximal radius, optionally extended by the user (**extend**),
and filled with implicit water at a defined temperature (**Temperature**) and
pH (**ph**). The DXBox volume is discretized in a grid of resolution
**mesh_size**. A rigid ligand (**ligand_path**) is introduced at the center of
the DXBox and rotated **number_of_rotations** times around its geometrical
center. For each rotational state, a correlation function computes the
complementarity of the protein surface to the ligand surface over all grid
points.

In a first step, an instance of
:class:`epitopsy.EnergyGrid.FFT.FFT_correlation_scoring` is created
(**shape_scoring**) by computing the conjugated reciprocal lattice of the
protein Van der Waals surface (**vdwbox.box**) by FFT (Fast Fourier
Transform). For any rotational state **angle_element** of the ligand, the
molecular surface **pqr_vdw** is projected to the grid *via*
:meth:`Structure.Structure_Template.snap_vdw_to_box` and passed to
:meth:`epitopsy.EnergyGrid.FFT.FFT_correlation_scoring.get_correlation`.
Within this method, the reciprocal lattice of the ligand molecular surface is
computed by :meth:`epitopsy.EnergyGrid.FFT.FFT_correlation_scoring.do_fft`,
multiplied by **shape_scoring** to get a scalar matrix in the reciprocal space,
which is transformed to the real space by
:meth:`epitopsy.EnergyGrid.FFT.FFT_correlation_scoring.do_ifft`. Only the real
part is returned, which represent the cross-correlation of the ligand molecular
surface to the protein molecular surface. The offset is deduced by
:meth:`epitopsy.EnergyGrid.FFT.FFT_correlation_scoring.shift_fft`, which calls
:meth:`epitopsy.EnergyGrid.FFT.FFT_correlation_scoring.shift_fft`.


.. _contents-of-module-calculation:

Module Contents
---------------

.. automodule:: epitopsy.EnergyGrid.calculation
    :members:
    :undoc-members:
    :show-inheritance:

