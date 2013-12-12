import os
import gzip
import shutil
import collections
import warnings

import numpy as np
import pandas as pd
import numba

from sklearn.metrics import roc_curve, auc

from Bio import AlignIO
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import Residue
from Bio.PDB import PDBIO
from Bio.PDB.PDBExceptions import PDBConstructionWarning

from epitopsy.cython.fast_dca import get_entropies, get_alignment_alphabet

warnings.simplefilter('ignore', PDBConstructionWarning)

## using numba leads to small differences in the accuracy of the results!!!
# probably because it breaks the values down to floats or something like that

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
        A dictionary, which maps the alignment positions [1,N] on the structure.
    '''
    alignment_alphabet = get_alignment_alphabet()
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
        if aln_char != '-' and aln_char != '.':
            # check if the sequence is already at the amino acids of the pdb
            if real_seq_start + aa_counter >= seq_start_id:
                # check if the following position is in the structure

                if counter in res_map:
                    # check if the amino acid is in the alphabet of interest
                    if aln_char in alignment_alphabet:
                        # check if the characters match
                        if aln_char == res_map[counter]:
                            seq_map[aln_pos+1] = counter # +1 -> python arrays
                counter += 1
            aa_counter += 1

    return seq_map

def depreciated_get_seq_map(aln, pdb, chain_id, pdb_seq, seq_res, pdb_res):
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
        A dictionary, wich maps the alignment positions [1,N] on the structure.
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
        for i, res in enumerate(model[chain_id]):
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

@numba.autojit
def res_res_distance(res1_coords, res2_coords):
    '''
    Args:
        res1_coords -> Residue 1 coordinates
        res2_coords -> Residue 2 coordinates

    Returns:
        Minimum distance between any atoms.
    '''
    n_res1_atoms = res1_coords.shape[0]
    n_res2_atoms = res2_coords.shape[0]
    min_dist = 1e9
    for i in range(n_res1_atoms):
        for j in range(n_res2_atoms):
            diff_ij = 0.
            for k in range(3):
                diff_ij = diff_ij + (res1_coords[i,k] - res2_coords[j,k])**2
            diff_ij = np.sqrt(diff_ij)
            if diff_ij < min_dist:
                min_dist = diff_ij

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
                res1_coords = np.array([atom.get_coord() for atom in pdb[i_model][chain_id][res_i]])
                res2_coords = np.array([atom.get_coord() for atom in pdb[i_model][chain_id][res_j]])
                any_distance = res_res_distance(res1_coords, res2_coords)
                ## store neighbor information first, it requires the previous
                # results
                if neighbor_map[aln_i, aln_j] == -1:
                    # case 1: yet no information available
                    if np.abs(res_i - res_j) > seq_distance:
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
                    if np.abs(res_i - res_j) > seq_distance:
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
                            if np.abs(res_i - res_j) > seq_distance:
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


def add_entropy(aln_path, pdb_path, result_pdb_path, theta):
    '''
    This function writes the entropy of the supplied alignment to the b factor
    of each atom/residue.

    Args:
        aln_path ->
        pdb_path -> path to the pdb file. The format has to be like this:
                    '<pdb_id>-<chain_id>-<uniprot_start>_<uniprot_end>-<pdb_start>_<pdb_end>-<aln_description>.pdb'
        result_pdb_path -> path to the entropy added pdb
        theta -> reweight sequences with identity threshold

    Returns:
        None.
    '''
    with gzip.open(aln_path) as f:
        aln = AlignIO.read(f, "fasta")
    pdb_name = pdb_path.replace('.pdb.gz','').split('-')[-1]
    chain_id = pdb_path.replace('.pdb.gz','').split('-')[1]
    uniprot_res = [ int(x) for x in pdb_path.replace('.pdb.gz','').split('-')[2].split('_')]
    pdb_res = [ int(x) for x in pdb_path.replace('.pdb.gz','').split('-')[3].split('_')]
    parser = PDBParser()
    with gzip.open(pdb_path) as f:
        pdb = parser.get_structure(pdb_name, f)
    pdb = repair_pdb(pdb,chain_id)
    pdb_seq = get_sequence(pdb, chain_id)
    # seq_map maps aln positions on the structure
    seq_map = get_seq_map(aln, pdb, chain_id, pdb_seq, uniprot_res, pdb_res)
    entropy = get_entropies(aln_path, theta=theta)
    # set b factor to 0
    delete_chains = []
    delete_waters = []
    for i_model in pdb:
        for chain in i_model:
            if chain.id != chain_id:
                delete_chains.append(chain.id)
            else:
                for res in chain:
                    if res.resname != 'HOH':
                        for atom in res:
                            atom.set_bfactor(0)
                    else:
                        delete_waters.append(res.id)
    # remove chains
    for i_model in pdb:
        for kill_chain in delete_chains:
            i_model.detach_child(kill_chain)

    # remove waters
    for i_model in pdb:
        for chain in i_model:
            for kill_res in delete_waters:
                chain.detach_child(kill_res)

    # set entropy
    pdb_entropy = []
    for i_model in pdb:
        for (aln_pos_prime,res_num) in seq_map.items():
            aln_pos = aln_pos_prime - 1 # python arrays
            aln_entropy = entropy[aln_pos]
            pdb_entropy.append(aln_entropy)
            res = pdb[i_model.id][chain_id][res_num]
            for atom in res:
                atom.set_bfactor(aln_entropy)

    # write structure
    io = PDBIO()
    io.set_structure(pdb)
    io.save(result_pdb_path)


def get_almost_conserved_columns(aln, cutoff=0.95):
    '''
    Args:
        aln -> biopython alignment

    Returns:
        A list with the columns, where the amino acid with the highest
        frequency has a frequency larger than the cutoff. The columns
        start with 1!
    '''
    N = len(aln[0])
    M = len(aln)

    cons_list = []
    for i in range(N):
        aa_list = []
        for j in range(M):
            aa_list.append(aln[j].seq.tostring()[i])

        max_freq = max(collections.Counter(aa_list).values()) / float(len(aa_list))
        if max_freq >= cutoff:
            cons_list.append(i+1)

    return cons_list


def get_predictions(dca, seq_map, n_dis, seq_dist, cons_list):
    '''
    Return the predicted pairs.

    Args:
        dca -> sorted! dca results in the form [i,j,value],  i,j in [1,..,N]
        seq_map -> dict that maps the aln positions to the sequence 
                    [1,..,N] -> [1,..,seq_len]
        n_dis -> number of values for the search
        seq_dist -> required sequence seperation of the predictions
        cons_list -> list of conserved columns, which should not be looked at
                    [1,...,N]
    Returns:
        A list [i_pdb, j_pdb, value], where *_seq refers to the coupled
        sequence positions and value is either DI or MI. Sequence positions
        are in the range [1,N].
    '''
    n_contacts = len(dca)
    # count iterated item to have no inf loop
    item_count = 0
    # count tps
    pred_count = 0
    pred_list = []
    while pred_count < n_dis and item_count < n_contacts:
        dca_item = dca[item_count]
        i = dca_item[0]
        j = dca_item[1]
        value = dca_item[2]
        # check if these columns are almost totaly conserved
        if i not in cons_list and j not in cons_list:
            if i in seq_map and j in seq_map:
                i_pdb = seq_map[i]
                j_pdb = seq_map[j]
                if np.abs(i_pdb-j_pdb) > seq_dist:
                    pred_list.append([i_pdb, j_pdb, value])
                    pred_count += 1
#        pred_count += 1
        item_count += 1

    return pred_list

def read_dca(di_path, method):
    '''
    Args:
        di_path -> path to the dca file
        method -> "mi" or "di"

    Returns:
        An unsorted list [i,j,di].
    '''
    with gzip.open(di_path) as f:
        di_content = f.readlines()

    dca = []
    # filter comments
    di_content = filter(lambda x: not x.startswith('#'), di_content)

    # first line has header information
    header_list = [ x.strip() for x in di_content[0].split(' ')]
    i_pos = header_list.index('i')
    j_pos = header_list.index('j')
    coevo_pos = header_list.index(method)

    for line in di_content[1:]:
        line_content = line.split(' ')
        i = int(line_content[i_pos])
        j = int(line_content[j_pos])
        di = float(line_content[coevo_pos])
        dca.append([i,j,di])

    return dca

def read_Meff_from_dca(di_path):
    '''
    Args:
        di_path -> path to the di result file

    Returns:
        A float Meff.
    '''
    with gzip.open(di_path) as f:
        line = f.next()

    first_split = line.split("Meff=")[-1]
    second_split = first_split.split(' ')[0]
    Meff = float(second_split)
    return Meff

def read_M_from_dca(di_path):
    '''
    Args:
        di_path -> path to the di result file

    Returns:
        A float Meff.
    '''
    with gzip.open(di_path) as f:
        line = f.next()

    first_split = line.split("M=")[-1]
    second_split = first_split.split(' ')[0]
    M = float(second_split)
    return M

def get_di_pair_matrix(i,j):
    '''
    Args:
        i -> alignment position [1,...,N]
        j -> alignment position [1,...,N]

    Returns:
        Numpy array q x q, with the size q of the alignment alphabet.
    '''
    wd = os.getcwd()
    protein_name = os.path.split(wd)[1]
    aln_path = "{0}_clean.fa.gz".format(protein_name)
    with gzip.open(aln_path) as f:
        aln = AlignIO.read(f, "fasta")
    
    aln_alphabet = get_alignment_alphabet()
    # python starts with 0
    for key in aln_alphabet:
        aln_alphabet[key] = aln_alphabet[key] - 1

    # python starts with 0
    i_pos = i - 1
    j_pos = j - 1

    n_aa = len(aln_alphabet)
    matrix_aa = np.zeros([n_aa,n_aa])
    gap = "-"
        
    for record in aln:
        sequence = record.seq.tostring()
        i_aa = sequence[i_pos]
        j_aa = sequence[j_pos]
    
        if i_aa not in aln_alphabet:
            i_index = aln_alphabet[gap]
        else:
            i_index = aln_alphabet[i_aa]
        if j_aa not in aln_alphabet:
            j_index = aln_alphabet[gap]
        else:
            j_index = aln_alphabet[j_aa]

        matrix_aa[i_index, j_index] += 1

    return matrix_aa

def format_object(obj, n_chars):
    obj_str = str(obj)
    n_chars_obj = len(obj_str)
    n_whitespace = n_chars - n_chars_obj
    whitespace = n_whitespace * " "
    return whitespace + obj_str


def show_di_pair_matrix(i,j):
    '''
    Args:
        i -> alignment position [1,...,N]
        j -> alignment position [1,...,N]

    Returns:
        None.
    '''
    matrix_aa = get_di_pair_matrix(i,j)
    x_dim = matrix_aa.shape[0]
    y_dim = matrix_aa.shape[1]
    max_value = int(matrix_aa.max())
    n_chars = len(str(max_value)) + 1
    n_aa = 2
        
    aln_alphabet = get_alignment_alphabet()
    # python starts with 0
    for key in aln_alphabet:
        aln_alphabet[key] = aln_alphabet[key] - 1
    alphabet_mapping = {}
    for k,v in aln_alphabet.items():
        alphabet_mapping[v] = k

    str_list = []
    header_list = [format_object("",n_aa)] + [format_object(value,n_chars) for value in alphabet_mapping.values()]
    header_str = "".join(header_list)
    str_list.append(header_str)
    for x in range(x_dim):
        line_list = [format_object(alphabet_mapping[x], n_aa)]
        for y in range(y_dim):
            line_list.append(format_object(int(matrix_aa[x,y]),n_chars))
        line_str = "".join(line_list)
        str_list.append(line_str)
        
    for line in str_list:
        print(line)


def get_ROC_data(prot_data_path):
    '''
    Args:
        prot_data_path -> path to the protein folder, that contains all methods
                            (dca_new, mi_new, ...) and the 
                            <prot>_DI_update.txt.gz files in each folder.

    Returns:
        A dictionary with all methods and their corresponding tpr and fpr.
    '''
    results = {}
    protein_name = os.path.split(prot_data_path)[1]
    file_list = os.listdir('.')
    folder_list = ["dca_old", "mi_old", "dca_new", "mi_new", "mi"]
#    folder_list = []
#    for file_item in file_list:
#        if os.path.isdir(file_item):
#            if "mi" in file_item or "dca" in file_item:
#                folder_list.append(file_item)
#    folder_list.sort()

    data_template = "{0}_clean_DI_update.txt.gz".format(protein_name)
    for folder in folder_list:
        di_data_path = os.path.join(prot_data_path, folder, data_template)
        with gzip.open(di_data_path) as f:
            data = pd.read_table(f, sep=" ")

        contact_label_nan = data["contact"]
        seq_neighbor = data["seq_neighbor"]
        if folder.startswith("dca"):
            predict_method = "di"
        elif folder.startswith("mi"):
            predict_method = "mi"
        else:
            raise AttributeError("Unkown folder type:\n{0}".format(folder))

        di_probas_nan = data[predict_method]
        contact_label = []
        di_probas = []
        for i, item in enumerate(contact_label_nan):
            if not np.isnan(item) and seq_neighbor[i] == 0:
                contact_label.append(contact_label_nan[i])
                di_probas.append(di_probas_nan[i])
        fpr, tpr, thresholds = roc_curve(contact_label, di_probas)
        roc_auc = auc(fpr, tpr)
        results[folder] = [fpr, tpr,roc_auc]

    return results

def plotDcaMap(filename, method='di', name=None):
    if(method not in ['mi', 'di']):
        raise Exception('The only supported methods are mi and di. Choose one.')            
    file = open(filename)
    lines = file.readlines()
    DataInfo = lines[0].split()
    dim = int(DataInfo[2].split('=')[1])
    MatrixDCA = np.zeros([dim,dim])
    for line in lines[2:]:
        a = line.split()
        if(method=='di'):
            MatrixDCA[int(a[0])-1,int(a[1])-1] = float(a[3])
        if(method=='mi'):
            MatrixDCA[int(a[0]-1),int(a[1])-1] = float(a[2])
    fig = plt.figure(num=None, figsize=(10, 8), dpi=120, facecolor='w', edgecolor='k')
    DCAMap = fig.add_subplot(111)
    #plt.xlabel('Contact A')
    #plt.xlabel('Contact B')
    if name:
        plt.title('DCA pair map of '+name, size=20)
    else:
        plt.title('DCA pair map', size=20)
    DCAMap = imshow(MatrixDCA)
    plt.annotate(DataInfo[2], xy=(dim/4, 2*dim/3), xytext=(dim/4, 2*dim/3), size='20', color='w')
    plt.annotate(DataInfo[3], xy=(dim/4, 4*dim/5), xytext=(dim/4, 4*dim/5), size='20', color='w')
    plt.colorbar()