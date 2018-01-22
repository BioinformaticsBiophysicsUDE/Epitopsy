
.. _EnergyGrid-index:

EnergyGrid package
==================

.. toctree::

    packages/EnergyGrid/calculation
    packages/EnergyGrid/analysis
    packages/EnergyGrid/tools
    packages/EnergyGrid/FFT

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
.. |pm| unicode:: 0xB1
.. |_| unicode:: 0xA0
   :trim:

.. _EnergyGrid-syntax:


Preparing input files
---------------------

Protein
^^^^^^^

Open the PDB file of your choice in the PyMOL software and perform the following steps:

* remove unnecessary cofactors, ligands, salts and water molecules
* for multimeric proteins, generate the biological unit through symmetry if necessary (command ``symexp sym, <protein_name>, <protein_name>, <cutoff_Angstroms>``)
* make sure solvent-exposed metallic centers (Zn(II), Ca(II), etc.) are properly water-coordinated, with properly oriented hydrogens on the water molecules
* add the ``OXT`` atom at the C-terminus if missing, unless the protein is a fragment, in which case both termini must be capped with neutral groups: ``ACE`` at the N-terminus (Build/Residue/Acetyl), ``NME`` at the C-terminus (Build/Residue/N-Methyl)
* identify residues with alternate side-chain conformations (``select alternate, not(alt '')``), select which conformation to keep, delete atoms from alternate conformations, then erase alternate information (``alter <protein_name>, alt=''``)
* instruct PyMOL to preserve the original atom order when saving PDB files (``set retain_order, on``)
* save the protein in a new PDB file

Run the PDB2PQR tool using the Amber forcefield to add information on van der Waals radii and atomic partial charges:

.. code-block::bash

    pdb2pqr --ff=amber --chain structure.pdb structure.pqr

Inspect the resulting PQR file in a text editor; the file header will list any error in the protein structure, in particular gaps in the amino acid sequence and missing heavy atoms in the side-chains. Such errors can be fixed in homology modeling tools, such as MODELLER. Warnings about missing parameters for ions are less critical; Epitopsy provides a wrapper function for PDB2PQR that will automatically re-introduce Zn(II), Ca(II) and sulfate ions removed by PDB2PQR. Other ions need to be re-introduced manually.

Caveats: proteins decorated with glycosides, lipids, non-proteinogenic amino acids and non-natural functional groups cannot be processed by PDB2PQR. In that scenario, the PQR file has to be created externally, using charges and radii obtained from specialized biomolecular forcefields, for example the GLYCAM forcefield, using the TLEaP tool instead of PDB2PQR.

Ligand
^^^^^^

For peptides, saccharides and lipids, suitable charges and radii can be extracted from the Amber and GLYCAM forcefields. For drug-like compounds, reasonable charges and geometries can be obtained from quantum chemical calculations in vacuum using the HF/6-31(d) level of theory. The Gasteiger-Marsili method is also an option, however please check the resulting partial charges add up to the expected total charge for the ligand.


Scanning
--------

The calculation is carried out in two steps. First the protein electrostatic potential is computed by running APBS on the protein PQR file. Then Epitopsy scans the APBS grid using the ligand PQR file to generate an Energy Grid. When scanning a library of *N* proteins with *M* ligands, the APBS step will be carried out *N* times and the scanning step *N* :math:`\times` *M* times. Although scans are based on rigid molecules, flexibility can be introduced by scanning multiple conformers of the same protein and/or ligand. The following subsections will provide you with the general procedure to treat each case.


Rigid protein and rigid ligand
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use :func:`epitopsy.EnergyGrid.calculation.scan`:

    >>> from epitopsy import EnergyGrid as EG
    >>> EG.electrostatics('protein.pdb', 'ligand.pqr', mesh_size=3*[0.8],
    ...                   center_pdb=True, cubic_box=False, verbose=True)
    >>> EG.scan('protein.pdb', 'ligand.pqr')


Rigid protein and semi-flexible ligand
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use :func:`epitopsy.EnergyGrid.calculation.scan_conformers`:

    >>> from epitopsy import EnergyGrid as EG
    >>> protein_path = 'protein.pdb'
    >>> ligand_paths = ['ligand_confA.pqr', 'ligand_confB.pqr']
    >>> elec_kwargs = {'mesh_size': 3 * [0.8, ], 'verbose': False}
    >>> scan_kwargs = {'verbose': False}
    >>> EG.scan_conformers(protein_path, ligand_paths,
    ...                    elec_kwargs=elec_kwargs,
    ...                    scan_kwargs=scan_kwargs,
    ...                    merge=True,
    ...                    ligand_weights=[0.23, 0.77])


Semi-flexible protein and semi-flexible ligand
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use :func:`epitopsy.EnergyGrid.calculation.scan_conformers` (note the `centering` argument from **elec_kwargs** transforms ``protein_confA.pdb``, while ``protein_confB.pdb`` get superposed):

    >>> from epitopsy import EnergyGrid as EG
    >>> protein_paths = ['protein_confA.pdb', 'protein_confB.pdb']
    >>> ligand_paths = ['ligand_confA.pqr', 'ligand_confB.pqr']
    >>> elec_kwargs = {'verbose': False, 'mesh_size': 3 * [0.8],
    ...                'centering': 'center+rotate'}
    >>> scan_kwargs = {'verbose': False}
    >>> EG.scan_conformers(protein_paths, ligand_paths,
    ...                    elec_kwargs=elec_kwargs,
    ...                    scan_kwargs=scan_kwargs,
    ...                    merge=True,
    ...                    protein_weights=[0.65, 0.35],
    ...                    ligand_weights=[0.23, 0.77])

Use :func:`epitopsy.EnergyGrid.calculation.merge_conformers` to re-weight a merged EG without re-scanning the entire library of conformers::

    >>> # perform an additional merge with different weights
    >>> EG.merge_conformers(protein_paths, ligand_paths,
    ...                     protein_weights=[0.8, 0.2],
    ...                     ligand_weights=[0.23, 0.77],
    ...                     output_fmt='merge2_{}.dx')

Library of proteins and ligands
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use :func:`epitopsy.EnergyGrid.calculation.scan_multi`:

    >>> import os
    >>> from epitopsy import EnergyGrid as EG
    >>> protein_paths = ['protein1.pdb', 'protein2.pdb']
    >>> ligand_paths = ['ligand1.pqr', 'ligand2.pqr']
    >>> elec_kwargs = {'mesh_size': 3 * [0.8, ], 'verbose': True}
    >>> scan_kwargs = {'verbose': True}
    >>> EnergyGrid.scan_multi(protein_paths, ligand_paths,
    ...                       elec_kwargs=elec_kwargs,
    ...                       scan_kwargs=scan_kwargs)


Visualization
-------------

The following PyMOL script will load the Epitopsy output files, generate isosurfaces and ray-trace the protein::

    >>> reinitialize
    >>> bg_color white
    >>> set ray_shadow, off
    >>> # load files
    >>> load protein.pqr
    >>> load protein_epi.dx
    >>> # generate protein surface and energy isosurfaces
    >>> as surface, protein
    >>> color white, protein
    >>> isosurface epi_pos_1kT, protein_epi, +1
    >>> isosurface epi_pos_2kT, protein_epi, +2
    >>> isosurface epi_neg_1kT, protein_epi, -1
    >>> isosurface epi_neg_2kT, protein_epi, -2
    >>> color skyblue, epi_pos_*kT
    >>> color firebrick, epi_neg_*kT
    >>> set transparency, 0.5, epi_*_1kT
    >>> # ray-trace
    >>> ray 1200


Red isosurfaces represent negative energies (favorable) while blue isosurfaces represent positive energies (defavorable). The translucent isosurface are drawn at 1 |_| |kbT| while the solid isosurfaces are drawn at 2 |_| |kbT|. These thresholds are arbitrary and usually depend on the charge of the ligand. For highly charged ligands such as DNA and RNA fragments, the |pm| |_| 1 |_| |kbT| level may be too small and would encapsulate the entire protein.

It can be helpful to plot the ESP potential on a separate image. The code is very similar, however the red and blue colors have a slightly different hue to distinguish electrostatic potential isosurfaces from energy isosurfaces::

    >>> load protein_esp.dx
    >>> isosurface esp_pos_1kT, protein_epi, +1
    >>> isosurface esp_neg_1kT, protein_epi, -1
    >>> color lightblue, epi_pos_*kT
    >>> color salmon, epi_neg_*kT









