
:mod:`calculation` --- Grid-based protein-ligand affinity screening
===================================================================

.. module:: calculation
   :synopsis: Grid-based protein-ligand affinity screening.
.. moduleauthor:: Christoph Wilms <christoph.wilms@uni-due.de>
.. sectionauthor:: Jean-Noel Grad <jean-noel.grad@uni-due.de>


This module provides a function to probe the surface of a protein with a
ligand, and return the area of space around the protein where the ligand
has the highest probability of being present (in solution).

.. note::

    In order to run :func:`calculate_interaction_energy()`, following
    dependencies must be present on the system:

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

* protein.pdb
* protein.pqr (only if you don't use PDB2PQR within Biopython)
* ligand.pqr
* ligand.mol2

To start the Epitopsy calculation::

    >>> from epitopsy.Structure import PDBFile
    >>> from epitopsy.tools.style import style
    >>> from epitopsy.functions.calculate_partition_function import calculate_interaction_energy
    >>> m = 0.65 # mesh size = grid resolution
    >>> calculate_interaction_energy(pdb_path = "sonic.pdb", # your pdb file
    ...     ligand_path = "lig.pqr",
    ...     mesh_size = [m,m,m],         # your grid resolution
    ...     number_of_rotations = 150,   # 150 is sufficient
    ...     extend = None,               # extension in every direction (in Angstroms)
    ...     use_pdb2pqr = True,          # if False, please provide a PQR file of your protein
    ...     center_pdb=True,             # centers the pdb to (0,0,0), set to True if use_pdb2pqr is True, False otherwise
    ...     box_type=["esp","vdw"],      # always use at least ["esp","vdw"], others can be added
    ...     cubic_box = True,            # cubic boxes are memory-consuming
    ...     zipped=False)                # PyMOL can't read zipped files yet

Following files are generated in the same directory:

* sonic.in (APBS input file)
* io.mc (APBS log file, showing the calculation progress)
* result_data.txt (various variables on the calculation)
* sonic_esp.dx (electrostatic potential)
* sonic_vdw.dx (van der Waals surface)
* gibbs_free_energy.dx (ligand probability of presence around the protein)
* counter_matrix.dx.gz

The most important file is :file:`gibbs_free_energy.dx`. It can be displayed
as a negative isosurface in PyMOL using the APBS Tools2 plugin. The isosurface
represents the probability of presence of the ligand around the protein. To
render the picture in PyMOL (:ref:`Figure 1 <rendering>`)::

    >>> reinitialize
    >>> load protein.pdb, protein
    >>> load gibbs_free_energy.dx, protein_gibbs
    >>> # display
    >>> cmd.bg_color('white')
    >>> cmd.hide('lines', 'protein')
    >>> cmd.show('cartoon', 'protein')
    >>> cmd.show('spheres', "protein and elem CA+ZN")
    >>> cmd.color('yellow', 'protein and not elem CA+ZN')
    >>> # APBS
    >>> cmd.isosurface('negative', 'protein_gibbs', -2)
    >>> cmd.color('gray50', 'negative')
    >>> cmd.set("transparency", 0.5, 'negative')

.. _rendering:
.. figure:: ../../_static/figures/0.25Angstroem-SHH-heparin.jpg
    :height: 1181px
    :width: 2000 px
    :scale: 40 %
    :alt: Rendering of PDB file 3M1N in cartoon representation with a grey
          cloud around its positively-charged residues (Lys and Arg).
    :align: center

    **Figure 1**: Probability of presence of heparin around Sonic Hedgehog.

    Negative potential isosurface for Sonic Hedgehog (PDB: 3M1N) probed by a
    heparine dissacharid (PDB: 1HPN, charges assigned by *ab initio*
    computation with the HF/6-31G** method) at 0.40 |Aring| grid resolution.
    Since heparin features a net negative charge of -3, the areas around Lys
    and Arg are the most probable sites for heparin docking. This picture
    clearly highlights the Cardin-Weintraub motif on the N-terminus of Sonic
    Hedgehog, as well as a secondary docking site on the globular domain.

.. note::

    PyMOL v1.5.0.1 (with plugin APBS Tools2, APBS v1.4.0, FFTW v3.0, on Ubuntu
    13.10) can crash the OS when rendering positive or negative *isosurfaces*
    (but *molecular surface* remains unaffected). The issue is fixed in the
    v1.6.0.0 release of PyMOL. You may want to upgrade your PyMOL installation
    as well as your Biopython pymol module (see :doc:`../../install/pymol`)
    before proceeding with the isosurface rendering in PyMOL.

Mathematical approach
---------------------

This module uses APBS to calculate the affinity of a ligand for a protein
molecular surface. The calculation takes into account electrostatics
interactions and desolvation energies. Contrary to MD simulations, which use
an explicit solvent model and a stochastic approach, Epitopsy uses an implicit
solvent model and covers the whole space available around the protein on
discretized grid points.

Probability of all states at position :math:`\vec{r}`:

    :math:`p(r) = \frac{1}{Z} \sum_i^N e^{\left ( \frac{-E_i}{k_B T} \right )}`

Probability of all states at position :math:`\vec{r}` in solution (with :math:`E_i = 0` for all rotations):

    :math:`p(r_{solv}) = \frac{1}{Z} \sum_i^N e^{\left ( \frac{-E_i}{k_B T} \right )} = \frac{N}{Z}`

Dissociation constant:

    :math:`K = \frac{\displaystyle p(r)}{\displaystyle p(r_{solv})} = \frac{1}{N} \sum_i^N e^{\left ( \frac{-E_i}{k_B T} \right )}`

Gibbs free energy:

    :math:`\Delta G = - k_B T ln(K)`

Use an infinite reference in solution, for which the probability of
all rotations is given as

    :math:`p(\vec{r}_{solv}) = \frac{1}{Z_K} \sum_i^N exp^{\left ( \frac{-E_i}{k_B T} \right ) } = \frac{N}{Z_K}`

The ratio of the probabilty around the protein and the one in solution
yields a equilibrium constant K

    :math:`K = \frac{\displaystyle p(\vec{r})}{\displaystyle p(\vec{r}_{solv})}`

This enable us to calculate the difference in Gibbs free energy

    :math:`\Delta G(solv \rightarrow prot) = - k_B T ln(K)`

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
:class:`FFTCorrelationScoring.FFT_correlation_scoring` is created
(**shape_scoring**) by computing the conjugated reciprocal lattice of the
protein Van der Waals surface (**vdwbox.box**) by FFT (Fast Fourier
Transform). For any rotational state **angle_element** of the ligand, the
molecular surface **pqr_vdw** is projected to the grid *via*
:meth:`Structure.Structure_Template.snap_vdw_to_box` and
passed to :func:`FFTCorrelationScoring.FFT_correlation_scoring.get_correlation`.
Within this method, the reciprocal lattice of the ligand molecular surface is
computed by :func:`FFTCorrelationScoring.FFT_correlation_scoring.do_fft`,
multiplied by **shape_scoring** to get a scalar matrix in the reciprocal space,
which is transformed to the real space by
:func:`FFTCorrelationScoring.FFT_correlation_scoring.do_ifft`. Only the real
part is returned, which represent the cross-correlation of the ligand molecular
surface to the protein molecular surface. The offset is deduced by
:func:`FFTCorrelationScoring.FFT_correlation_scoring.shift_fft`, which calls
:func:`FFTCorrelationScoring.FFT_correlation_scoring.shift_fft`.

Benchmarking
------------

A mesh size of 0.40 |Aring| usually give satisfactory results. A mesh size of
0.50 |Aring| or above will generate artifacts in the form of small bubbles
(they disappear above 5.00 |Aring|, but the picture will be crude).
At a resolution of 0.40 |Aring|, each .dx file will require between 300 and
400 Mio of free disk space (each experiment generates at least three .dx
files). The function :func:`calculate_interaction_energy()` should not be run
below 0.25 |Aring|, since the volume of data generated is proportional to the
inverse of the grid resolution at the power 3.

For high grid resolution (**mesh_size** > 1.0 |_| |Aring|), the calculation
can be executed in a matter of minutes, but for a grid resolution of 0.50 |_|
|Aring| or less, the execution time will increase rapidly (:ref:`Figure 2
<benchmarking>`, graph a), especially for the Numpy implementation of FFT.
Using anfft instead of PyFFTW may be preferred for some automated analyzes,
namely alanine test scans, since it runs 5 |_| min faster at 0.40 |_| |Aring|
(:ref:`Figure 2 <benchmarking>`, graph b). For more casual use of Epitopsy,
pyFFTW can be used instead.

The 3-dimensional DFT has a complexity of |THgr|\ (n\ :sup:`6`\ ), which can
be brought down to |THgr|\ (n\ :sup:`3`\.log(n\ :sup:`3`\ )) using a FFT
algorithm. According to our benchmarking on PDB: 3M1N, the execution time is
approximately proportional to the inverse of the grid resolution at the
power 3 in the working range 0.25 |_| |leq| |_| **mesh_size** |_| |leq| |_|
2.0 |_| |Aring| (:ref:`Figure 2 <benchmarking>`, graph c).
For **mesh_size** |_| |geq| 2.0 |_| |Aring|, the calculation is too
approximative to be exploited, and for **mesh_size**
|_| |leq| 0.20 |_| |Aring|, the gain in refinement doesn't bring
much valuable information (and the APBS calculation consumes too much memory).
The picture quality is already sufficient at 0.40 |_| |Aring|, hence the
inverse power 3 model can be used as a rule of thumb when planing the execution
time of your experiments. Don't forget the additional calculatiom required at
the beginning by APBS to generate the electrostatic potential around the
protein. This calculation can't be multi-threaded and took approximately
5 |_| min at 0.50 |_| |Aring| and 40 |_| min at 0.25 |_| |Aring| on our system.

The anfft and Numpy calculations are identical until the 7\ :sup:`th` decimal
position. In two simulations at a resolution of 0.50 |_| |Aring|, a 100 |_|
MiB :file:`gibbs_free_energy.dx` file was generated by Numpy and anfft,
with only 44 non-matching lines over 2'396'366 lines (compare the non-matching
lines in :download:`anfft.dx <../../_static/files/anfft.dx>` and
:download:`numpy.dx <../../_static/files/numpy.dx>`; files obtained with
``grep -Fxvf numpy.dx anfft.dx``).

.. _benchmarking:
.. figure:: ../../_static/figures/benchmarking.*
    :target: ../../_static/figures/benchmarking.pdf
    :width: 1200 px
    :alt: Benchmarking three FFT algorithms. Numpy is clearly the slowest one,
          anfft and pyfftw are quite close, anfft being the quicker.
    :align: center

    **Figure 2**: Implementation of three FFT algorithms.

    Legend: Numpy (blue, on 1 core), pyFFTW (green, on 4 cores), anfft (red,
    on 4 cores). Benchmarking on the Intel Core 2 Quad Processor Q9650
    (4 cores, 3.00 GHz) with PDB: 3M1N and a heparin disaccharid as probe.
    (**a**) Numpy lacks multiprocessing capabilities and runs slowly. PyFFTW
    v0.9.2 and anfft v0.2 (compiled both with FFTW v3.3.3) give better results.
    (**b**) Zoom with a power model fitted on the data. The orange curve shows
    the absolute difference in execution time.
    (**c**) Linearization of the data using a power model (only linear for the
    subset {0.25 |_| |leq| |_| **mesh_size** |_| |leq| |_| 2.0 |_| |Aring|}.


Profiling the code (:ref:`Table 1 <profiling>`) reveals that a significant
amount of time is lost to :func:`APBS.APBSWrapper.get_dxbox` and
:func:`DXFile.DXBox.write`, which both aren't paralellizable.

.. table:: **Table 1**: Code profiling of :func:`calculate_interaction_energy`
           on an Intel Core 2 Quad Processor Q9650 (4 cores, 3.00 GHz) with
           anfft, PDB: 3M1N and a heparin disaccharide as probe. Total time:
           212 |_| s, Timer unit: 1e-06 |_| s.
    :class: right-align-col
    :name: profiling

    ===================================================================== ==== ==== ======== ============ ====== =====
    Directive                                                             Line Hits Time (s) Per Hit (ms) % Time Cores
    ===================================================================== ==== ==== ======== ============ ====== =====
    :meth:`APBS.APBSWrapper.get_dxbox`                                     211    1     52.3      52354.3   24.7   1
    :meth:`Structure.Structure_Template.snap_vdw_to_box`                   303  150      3.8         25.1    1.8   1
    :func:`FFTCorrelationScoring.FFT_correlation_scoring.get_correlation`  310  150     31.0        206.8   14.6   4
    :func:`FFTCorrelationScoring.FFT_correlation_scoring.shift_fft`        312  150      3.0         20.0    1.4   4
    :func:`FFTCorrelationScoring.FFT_correlation_scoring.get_correlation`  315  150     31.0        206.7   14.6   4
    :func:`FFTCorrelationScoring.FFT_correlation_scoring.shift_fft`        317  150      3.0         20.0    1.4   4
    :func:`CalculateInteractionEnergyCython.count_rotations`               323  150      3.3         21.8    1.5   1
    :func:`np.exp`                                                         325  150     29.2        194.9   13.8   4
    :meth:`DXFile.DXBox.write`                                             417    1      3.9       3931.1    1.9   1
    :meth:`DXFile.DXBox.write`                                             424    1     44.8      44765.4   21.1   1
    :func:`calc_vol_and_surface.calc_vol_and_surface`                      433    1      3.3       3338.3    1.6   1
    ===================================================================== ==== ==== ======== ============ ====== =====


.. _contents-of-module-calculation:

Module Contents
---------------

.. automodule:: epitopsy.EnergyGrid.calculation
    :members:
    :undoc-members:
    :show-inheritance:
