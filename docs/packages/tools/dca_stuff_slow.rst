
:mod:`dca_stuff_slow` --- Direct Coupling Analysis routines
===========================================================

.. module:: dca_stuff_slow
   :synopsis: Direct Coupling Analysis routines.
.. moduleauthor:: Christoph Wilms <christoph.wilms@uni-due.de>
.. sectionauthor:: Jean-Noel Grad <jean-noel.grad@uni-due.de>


This module provides additional direct coupling analysis operations.


.. _dca_stuff_slow-syntax:

Module Syntax
-------------

Empty.

.. _contents-of-module-dca_stuff_slow:

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

    :returns: a dictionary, wich maps the alignment positions on the structure.
    
.. function:: repair_pdb(pdb, chain_id)

    :param pdb: biopython structure
    :param chain_id: chain, which contains the correct sequence

    :returns: repaired pdb object.
    
.. function:: get_sequence(pdb, chain_id)

    :param pdb: biopython structure
    :param chain_id: chain, which contains the correct sequence

    :returns: a list of the amino acids, which build up the sequence.
    
.. function:: get_res_id_aa_dict(pdb, chain_id)

    :returns: a dictionary of the first chain, which contains the residue id as
      key and the corresponding amino acid as value.
    
.. function:: res_res_distance(res1, res2)

    :param res1: residue 1
    :param res2: residue 2 

    :returns: minimum distance between any atoms.
    
.. function:: get_structure_map(contact_map, dist_map, neighbor_map, pdb, chain_id, seq_map, min_distance=8, seq_distance=5)

    :param contact_map: add a contact
    :param dist_map: add the corresponding distance
    :param neighbor_map: add sequence neighbor information
    :param pdb: biopython structure
    :param chain_id: chain, which contains the correct sequence
    :param min_distance: minimum distance between any two atoms (Morcos2011)

    :returns: a list containing three numpy arrays:
	
            #. contact map: 0's (no contact) and 1's (contact), where 
               a contact is defined as either a c_alpha - c_alpha distance
               less than 0.8 nm or any atom - atom distance less than 0.5 nm.
            #. min distance map: min distance between all the residues
            #. neighbor map: is the found minimal distance pair a neighbor (1)
               or not (0).
    
