
:mod:`dca_stuff` --- Direct Coupling Analysis routines
======================================================

.. module:: dca_stuff
   :synopsis: Direct Coupling Analysis routines.
.. moduleauthor:: Christoph Wilms <christoph.wilms@uni-due.de>
.. sectionauthor:: Jean-Noel Grad <jean-noel.grad@uni-due.de>


This module provides additional direct coupling analysis operations.


.. _dca_stuff-syntax:

Module Syntax
-------------

Empty.

.. _contents-of-module-dca_stuff:

Module Contents
---------------

.. function:: get_seq_map(aln, pdb, chain_id, pdb_seq, seq_res, pdb_res)

    Find the matching sequence in the alignment and return a dictionary, which
    maps the alignment positions to the residues of the structure.

    :param aln: alignment
    :param pdb: biopython structure object
    :param chain_id: chain id of interest
    :param pdb_seq: pdb sequence
    :param seq_res: start and end residue ids of the aln sequence
    :param pdb_res: start and end residue ids of the pdb sequence

    :returns: a dictionary, wich maps the alignment positions [1,N] on the structure.
    
.. function:: depreciated_get_seq_map(aln, pdb, chain_id, pdb_seq, seq_res, pdb_res)

    .. deprecated:: v2013
        Use :func:`.get_seq_map` instead.
    
.. function:: repair_pdb(pdb, chain_id)

    :param pdb: biopython structure
    :param chain_id: chain, which contains the correct sequence

    :returns: Repaired pdb object.
    
.. function:: get_sequence(pdb, chain_id)

    :param pdb: biopython structure
    :param chain_id: chain, which contains the correct sequence

    :returns: A list of the amino acids, which build up the sequence.
    
.. function:: get_res_id_aa_dict(pdb, chain_id)

    :param pdb: biopython structure
    :param chain_id: chain, which contains the correct sequence

    :returns: a dictionary of the first chain, which contains the residue id as
      key and the corresponding amino acid as value.
    
.. function:: res_res_distance(res1_coords, res2_coords)

    :param res1_coords: Residue 1 coordinates
    :param res2_coords: Residue 2 coordinates

    :returns: Minimum distance between any atoms.
    
.. function:: get_structure_map(contact_map, dist_map, neighbor_map, pdb, chain_id, seq_map, min_distance=8, seq_distance=5)

    :param contact_map: add a contact
    :param dist_map: add the corresponding distance
    :param neighbor_map: add sequence neighbor information
    :param pdb: biopython structure
    :param chain_id: chain, which contains the correct sequence
    :param min_distance: minimum distance between any two atoms (Morcos2011)

    :returns: A list containing three numpy arrays.
    
    Structure of the three arrays:
        #. contact map :math:`\rightarrow` 0's (no contact) and 1's (contact),
           where a contact is defined as either a c_alpha - c_alpha distance
           less than 0.8 nm or any atom - atom distance less than 0.5 nm.
        #. min distance map :math:`\rightarrow` min distance between all the residues
        #. neighbor map :math:`\rightarrow` is the found minimal distance pair
           a neighbor (1) or not (0).
    
.. function:: add_entropy(aln_path, pdb_path, result_pdb_path, theta)

    This function writes the entropy of the supplied alignment to the b factor
    of each atom/residue.

    :param aln_path: unknown
    :param pdb_path: path to the pdb file (see below)
    :param result_pdb_path: path to the entropy added pdb
    :param theta: reweight sequences with identity threshold

    :returns: None.
    
    Note: **pdb_path** has to be in the form: '<pdb_id>-<chain_id>-<uniprot_start>_<uniprot_end>-<pdb_start>_<pdb_end>-<aln_description>.pdb'
    
.. function:: get_almost_conserved_columns(aln, cutoff=0.95)

    :param aln: biopython alignment

    :returns: A list with the columns, where the amino acid with the highest
      frequency has a frequency larger than the cutoff. The columns
      start with 1!
    
.. function:: get_predictions(dca, seq_map, n_dis, seq_dist, cons_list)

    Return the predicted pairs.

    :param dca: sorted! dca results in the form [i,j,value],  i,j in [1,...,N]
    :param seq_map: dict that maps the aln positions to the sequence
      [1,...,N] :math:`\rightarrow` [1,...,seq_len]
    :param n_dis: number of values for the search
    :param seq_dist: required sequence seperation of the predictions
    :param cons_list: list of conserved columns, which should not be looked at
      [1,...,N]

    :returns: a list [i_pdb, j_pdb, value], where *_seq refers to the coupled
      sequence positions and value is either DI or MI. Sequence positions
      are in the range [1,N].
    
.. function:: read_dca(di_path, method)

    :param di_path: path to the dca file
    :param method: "mi" or "di"

    :returns: an unsorted list [i,j,di].
    
.. function:: read_Meff_from_dca(di_path)

    :param di_path: path to the di result file

    :returns: a float Meff.
    
.. function:: read_M_from_dca(di_path)

    :param di_path: path to the di result file

    :returns: a float Meff.
    
.. function:: get_di_pair_matrix(i,j)

    :param i: alignment position [1,...,N]
    :param j: alignment position [1,...,N]

    :returns: Numpy array :math:`q \times q`, with the size q of the alignment alphabet.
    
.. function:: format_object(obj, n_chars)

    Docstring missing.

.. function:: show_di_pair_matrix(i,j)

    :param i: alignment position [1,...,N]
    :param j: alignment position [1,...,N]

    :returns: ``None``
    
.. function:: get_ROC_data(prot_data_path)

    :param prot_data_path: path to the protein folder, that contains all methods
      (dca_new, mi_new, ...) and the 
      <prot>_DI_update.txt.gz files in each folder.

    :returns: a dictionary with all methods and their corresponding tpr and fpr.
    
