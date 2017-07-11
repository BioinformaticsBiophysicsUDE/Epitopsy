********
Tutorial
********

FFT-accelerated energy scan
===========================

Rigid protein and rigid ligand
------------------------------

Prepare a protein structure as a PDB file and a ligand structure as a PQR file,
with proper partial charges and van der Waals radii. These can be obtained from
a biomolecular forcefield or derived using specialized tools, such as GAFF
or the Prodrug server. Protein charges and radii will be automatically added
by PDB2PQR, with limited support for ions (currently Zn(II), Ca(II) and sulfate
are supported). Select a reasonable grid size, for example 0.8 Angstroms.

To start the electrostatics calculation and the FFT scan::

    >>> from epitopsy import EnergyGrid as EG
    >>> m = 0.80
    >>> EG.electrostatics('protein.pdb', 'ligand.pqr', mesh_size=[m,m,m],
    ...                   center_pdb=True, cubic_box=False, verbose=True)
    >>> EG.scan('protein.pdb', 'ligand.pqr', number_of_rotations=150)


Rigid protein and semi-flexible ligand
--------------------------------------

Prepare a protein structure as a PDB file and several ligand structures as PQR
files. These ligand structure may be obtained from a clustered Molecular
Dynamics trajectory or from a conformational search. Then proceed with::

    >>> from epitopsy import EnergyGrid as EG
    >>> m = 0.8
    >>> protein_paths = ['protein-centered.pdb']
    >>> ligand_paths = ['cluster_1.pqr', 'cluster_2.pqr']
    >>> elec_kwargs = {'mesh_size': [m,m,m], 'cubic_box': False,
    ...               'center_pdb': False, 'rotate_pdb': False,
    ...               'use_pdb2pqr': False, 'verbose': True, 'extend': 30}
    >>> scan_kwargs = {'verbose': True, 'number_of_rotations': 150}
    >>> # compute an electrostatic map for the protein and scan each ligand
    >>> EG.scan_multiconformational(protein_paths, ligand_paths,
    ...                             elec_kwargs=elec_kwargs,
    ...                             scan_kwargs=scan_kwargs)
    >>> # merge the 2 grids with weights derived from a MD simulation
    >>> EG.merge_multiconformational(protein_paths, ligand_paths,
    ...                              ligand_weights=(0.4, 0.6))


Semi-flexible protein and semi-flexible ligand
----------------------------------------------

Prepare several protein structures as PDB files and several ligand structures
as PQR files. These protein structures may be obtained from X-ray structures
or from a clustered Molecular Dynamics trajectory. Then proceed with::

    >>> import os
    >>> from epitopsy import EnergyGrid as EG
    >>> m = 0.6
    >>> ligand_paths  = [os.path.realpath('Ahx-open.pqr'),
    ...                  os.path.realpath('Ahx-folded.pqr')]
    >>> protein_paths = [os.path.realpath('1krn.pdb'),
    ...                  os.path.realpath('2pk4.pdb'),
    ...                  os.path.realpath('4duu.pdb')]
    >>> elec_kwargs = {'mesh_size': [m,m,m], 'verbose': False,
    ...                'center_pdb': False, 'rotate_pdb': False}
    >>> scan_kwargs = {'verbose': True, 'number_of_rotations': 150}
    >>> EnergyGrid.scan_multiconformational(protein_paths, ligand_paths,
    ...                                     elec_kwargs, scan_kwargs,
    ...                                     align_structures=False)
    >>> # merge all 6 grids with equal weights
    >>> EnergyGrid.merge_multiconformational(protein_paths, ligand_paths)





