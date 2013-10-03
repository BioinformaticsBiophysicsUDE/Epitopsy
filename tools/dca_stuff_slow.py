import os
import shutil
import numpy as np

from Bio import AlignIO
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import Residue

def get_seq_map(aln, pdb, chain_id, pdb_seq, seq_res, pdb_res):
    '''
    Find the matching sequence in the alignment and return a dictionary, which
    maps the alignment positions to the residues of the structure.

    Args:
        aln -> Alignment
        pdb -> biopython structure object
        chain_id -> chain id of interest
        pdb_seq -> pdb sequence
        seq_res -> start and end residue ids of the aln sequence
        pdb_res -> start and end residue ids of the pdb sequence 

    Returns:
        A dictionary, wich maps the alignment positions on the structure.
    '''
    seq_id = pdb.id
    res_map = get_res_id_aa_dict(pdb, chain_id)

    ref_seq = None
    ## search ref seq id
    seq_start_id = seq_res[0]
    seq_end_id = seq_res[1]
    pdb_start_id = pdb_res[0]
    pdb_end_id = pdb_res[1]
    # it can happen, that the end residue is not exactly the same as in the
    # pdb instructions
    for record in aln:
        rec_id = record.id.split('/')[0]
        residues_str = record.id.split('/')[1].split('-')
        residues = [int(x) for x in residues_str]
        if rec_id == seq_id:
            if seq_res[0] >= residues[0] and seq_res[1] <= residues[1]:
                ref_seq = record.seq.tostring()
                # sequences do not always start with the pfam databse information!
                real_seq_start = residues[0]
                break

    if ref_seq is None:
        raise ValueError("Could not find sequence of '{0}' and sequence residues '{1}-{2}' with pdb residues '{3}-{4}'".format(seq_id, seq_res[0], seq_res[1], pdb_res[0], pdb_res[1]))

    seq_map = {}
    counter = pdb_start_id
    aa_counter = 0
    # two counters, one counts the pdb residues, the other counts the sequence residues
    for aln_pos, aln_char in enumerate(ref_seq):
        if aln_char != '-':
            # check if the sequence is already at the amino acids of the pdb
            if real_seq_start + aa_counter >= seq_start_id:
                # check if the following position is in the structure
               
                if counter in res_map:
                    # check if the characters match
                    if aln_char == res_map[counter]:
                        seq_map[aln_pos+1] = counter # +1 -> python arrays
                counter += 1
            aa_counter += 1

    return seq_map

three2oneletter = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'ASN': 'N',
                   'GLN': 'Q', 'LYS': 'K', 'ILE': 'I', 'PRO': 'P',
                   'THR': 'T', 'PHE': 'F', 'ALA': 'A', 'GLY': 'G',
                   'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
                   'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M',
                   'MSE':'M'}

def repair_pdb(pdb, chain_id):
    '''
    Args:
        pdb -> biopython structure
        chain_id -> chain, which contains the correct sequence

    Returns:
        Repaired pdb object.
    '''
    for model in pdb:
        damaged_res = []
        # look for non amino acids classified as hetatm
        for i, res in enumerate(pdb[0][chain_id]):
            res_name = res.resname
            if res_name == "MSE" and res.id[0] != " ":
                damaged_res.append([i,res])
        
        # reverse the list, so the last items will be corrected first
        damaged_res.reverse()
        for (i,old_res) in damaged_res:
            new_id = list(old_res.id)
            new_id[0] = " "
            new_id = tuple(new_id)
            new_res = Residue.Residue(new_id, old_res.resname, '')
            for atom in old_res:
                new_res.add(atom)
            pdb[model.id][chain_id].detach_child(old_res.id)
            pdb[model.id][chain_id].add(new_res)

        # sort child list
        pdb[model.id][chain_id].child_list.sort(key=lambda x: x.id[1])
    
    return pdb

def get_sequence(pdb, chain_id):
    '''
    Args:
        pdb -> biopython structure
        chain_id -> chain, which contains the correct sequence

    Returns:
        A list of the amino acids, which build up the sequence.
    '''
    seq_list = []
    for i,res in enumerate(pdb[0][chain_id]):
        res_name = res.resname 
        if res_name in three2oneletter:
            seq_list.append(three2oneletter[res.resname])

    return seq_list

def get_res_id_aa_dict(pdb, chain_id):
    '''
    Returns:
        A dictionary of the first chain, which contains the residue id as
        key and the corresponding amino acid as value.
    '''
    res_map = {}
    non_amino_acid = False
    non_amino_acid_list = []
    for res in pdb[0][chain_id]:
        res_id = res.get_id()[1]
        res_name = res.resname
        # only add the amino acid, if it is one of proteogen amino acids
        # in that way no calciums etc. will be added
        if res_name in three2oneletter:
            aa = three2oneletter[res_name]
            res_map[res_id] = aa
        else:
            if res_name != "HOH" and res_name:
                non_amino_acid = True
                non_amino_acid_list.append(res_name)

#    if non_amino_acid is True:
#        print('Encountered the following non amino acids: {0}'.format(non_amino_acid_list))

    return res_map

def res_res_distance(res1, res2):
    '''
    Args:
        res1 -> Residue 1
        res2 -> Residue 2 

    Returns:
        Minimum distance between any atoms.
    '''
    coord1_list = []
    for atom in res1:
        coord1_list.append(atom.get_coord())

    coord2_list = []
    for atom in res2:
        coord2_list.append(atom.get_coord())

    coord1_list = np.array(coord1_list)
    coord2_list = np.array(coord2_list)

    min_dist = np.inf
    for coord1 in coord1_list:
        dist_list = coord2_list - coord1
        dist_list = np.sqrt(dist_list[:,0]**2 + dist_list[:,1]**2 + dist_list[:,2]**2)
        dist_min = np.min(dist_list)
        if dist_min < min_dist:
            min_dist = dist_min

    return min_dist


def get_structure_map(contact_map, dist_map, neighbor_map, pdb, chain_id,
        seq_map, min_distance=8, seq_distance=5):
    '''
    Args:
        contact_map -> add a contact
        dist_map -> add the corresponding distance
        neighbor_map -> add sequence neighbor information
        pdb -> biopython structure
        chain_id -> chain, which contains the correct sequence
        min_distance -> minimum distance between any two atoms
                        (Morcos2011)


    Returns:
        A list containing three numpy arrays:
            0. contact map -> 0's (no contact) and 1's (contact), where 
                a contact is defined as either a c_alpha - c_alpha distance
                less than 0.8 nm or any atom - atom distance less than 0.5 nm.
            1. min distance map -> min distance between all the residues
            2. neighbor map -> is the found minimal distance pair a neighbor (1)
                or not (0).
    '''
    res_map = seq_map.values()
    aln_map = seq_map.keys()

    n_residues = len(res_map)
    for i in range(n_residues):
        res_i = res_map[i]
        aln_i = aln_map[i] - 1
        for j in range(i+1,n_residues):
            res_j = res_map[j]
            aln_j = aln_map[j] - 1
            # iterate over all models:
            for i_model, model in enumerate(pdb):
                # calculate minimal distance
                any_distance = res_res_distance(pdb[i_model][chain_id][res_i], pdb[i_model][chain_id][res_j])
                ## store neighbor information first, it requires the previous 
                # results
                if neighbor_map[aln_i, aln_j] == -1:
                    # case 1: yet no information available
                    if np.abs(res_i - res_j) >= seq_distance:
                        # no neighbor
                        neighbor_score = 0
                        neighbor_map[aln_i, aln_j] = neighbor_score
                        neighbor_map[aln_j, aln_i] = neighbor_score
                    else:
                        # neighbor
                        neighbor_score = 1
                        neighbor_map[aln_i, aln_j] = neighbor_score
                        neighbor_map[aln_j, aln_i] = neighbor_score
                elif neighbor_map[aln_i, aln_j] == 1: 
                    # case 2: already found a neighbor for these positions
                    # only override if the new information has new information
                    if np.abs(res_i - res_j) >= seq_distance:
                        # no neighbor, but only use this information if this
                        # new information is "better", with better meaning,
                        # that this non neighbor is also a contact or closer
                        # than the previous information
                        if (any_distance < min_distance
                                or any_distance <= dist_map[aln_i,aln_j]):
                            neighbor_score = 0
                            neighbor_map[aln_i, aln_j] = neighbor_score
                            neighbor_map[aln_j, aln_i] = neighbor_score

                elif neighbor_map[aln_i, aln_j] == 0:
                    # case 3: already found a non neighbor for these positions
                    if contact_map[aln_i, aln_j] == 1:
                        # there is already a contact reported for these
                        # positions, so do not do anything
                        pass
                    else:
                        if any_distance < min_distance:
                            # these positions are in contact
                            if np.abs(res_i - res_j) >= seq_distance:
                                # no neighbor
                                neighbor_score = 0
                                neighbor_map[aln_i, aln_j] = neighbor_score
                                neighbor_map[aln_j, aln_i] = neighbor_score
                            else:
                                # neighbor
                                neighbor_score = 1
                                neighbor_map[aln_i, aln_j] = neighbor_score
                                neighbor_map[aln_j, aln_i] = neighbor_score


                ## if the new distance is less than the previous one, or if there
                ## is no distance information available yet
                if dist_map[aln_i, aln_j] == -1 or any_distance < dist_map[aln_i, aln_j]:
                    dist_map[aln_i,aln_j] = any_distance
                    dist_map[aln_j, aln_i] = any_distance

                ## store contact information
                if any_distance < min_distance:
                    contact_score = 1.
                    contact_map[aln_i,aln_j] = contact_score
                    contact_map[aln_j,aln_i] = contact_score
                else:
                    # no contact is also an information, but only, if no contact
                    # has been reported so far
                    if contact_map[aln_i,aln_j] == -1:
                        contact_score = 0
                        contact_map[aln_i,aln_j] = contact_score
                        contact_map[aln_j,aln_i] = contact_score


    return [contact_map, dist_map, neighbor_map]
