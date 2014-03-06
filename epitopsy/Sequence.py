import collections

import numpy as np

from Bio import AlignIO, SeqIO, Seq, SeqRecord, pairwise2
from Bio.SubsMat import MatrixInfo as matlist
import urllib, urllib2
import re

alignment_alphabet = {'-':1, 'A':2, 'C':3, 'D':4, 'E':5, 'F':6, 'G':7, 'H':8,
                      'I':9, 'K':10, 'L':11, 'M':12, 'N':13, 'P':14, 'Q':15,
                      'R':16, 'S':17, 'T':18, 'V':19, 'W':20, 'Y':21}

def make_alignment_upper_case(orig_aln):
    '''
    Uses the alignment and transforms all characters to upper case.

    Args:
        orig_aln -> biopython alignment

    Returns:
        An upper case biopython alignment
    '''
    #### copy alignment
    aln = []
    for record in orig_aln:
        aln.append(SeqRecord.SeqRecord(Seq.Seq(record.seq.tostring().upper()), record.id, description=record.description))
    return aln

def make_alignment_clean(orig_aln):
    '''
    Uses the alignment, removes sequences with non standard amino acids and
    reduces the alignment size if unnecessary gaps are there (this may happen,
    if a sequence has been removed).

    Args:
        orig_aln -> biopython alignment

    Returns:
        A clean biopython alignment
    '''
    #### acceptable amino acids
    aa_list = ['-', 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I',
            'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S',
            'T', 'V', 'W', 'Y']

    #### copy alignment
    aln = []
    for record in orig_aln:
        aln.append(SeqRecord.SeqRecord(Seq.Seq(record.seq.tostring()), record.id, description=record.description))

    n_dirty = len(aln[0])
    m_dirty = len(aln)
    print("input size: {0} x {1} (#seq x len(seq))".format(m_dirty, n_dirty))
    #### remove sequences with unkown amino acids
    clean_aln = []
    trash_aln = []
    for record in aln:
        clean_status = True
        for char in record:
            if char not in aa_list:
                clean_status = False
                break

        if clean_status is True:
            clean_aln.append(record)
        else:
            trash_aln.append(record)

    #### look for coloumns with only gaps
    aln_positions = len(clean_aln[0])

    kill_pos = [-1]
    for i in range(aln_positions):
        found_aa = False
        for record in clean_aln:
            if record[i] != '-':
                found_aa = True
                break

        if found_aa is False:
            kill_pos.append(i)

    kill_pos.append(len(record))

    for record in clean_aln:
        rec_seq = record.seq.tostring()
        new_seq = ""
        for i in range(len(kill_pos)-1):
            new_seq = new_seq + rec_seq[kill_pos[i]+1:kill_pos[i+1]]

        record.seq = Seq.Seq(new_seq)

    n_clean = len(clean_aln[0])
    m_clean = len(clean_aln)
    print("output size: {0} x {1}".format(m_clean, n_clean))

    return clean_aln

def format_stockholm_to_fasta(input_path, output_path):
    '''
    Change the format of an stockholm formatted alignment ot an fasta formatted
    alignment.

    Args:
        input_path -> path to the stockholm formatted format
        output_path -> path to the new fasta alignment
    '''
    aln = AlignIO.read(input_path, 'stockholm')
    with open(output_path, "w") as f:
        AlignIO.write(aln, f, 'fasta')


def get_pairwise_alignment(ref_seq, query_seq):
    """
    Args:
        ref_seq -> string of the reference sequence
        query_seq -> string of the query sequence

    Returns:
        A dictionary with keys ("ref_seq", "query_seq") and the corresponding
        aligned sequences. It uses the first calculated alignment, equally
        scored alignments are not considered.
    """
    result = {}
    matrix = matlist.blosum62
    ref = re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]",'',ref_seq.upper())
    if ref != ref_seq:
        print "Warning: removed gaps and/or non-natural amino acids in ref_seq."
    query = re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]",'',query_seq.upper())
    if query != query_seq:
        print "Warning: removed gaps and/or non-natural amino acids in query_seq."
    aln = pairwise2.align.globaldx(ref, query, matrix)
    result["ref_seq"] = aln[0][0]
    result["query_seq"] = aln[0][1]
    return result

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

def uniprot_mapping(fromtype, totype, identifier):
    """Takes an identifier, and types of identifier
    (to and from), and calls the UniProt mapping service.
    To get the right identifiers, visit 
    http://www.uniprot.org/faq/28#conversion.

    Args: 
        fromtype -> type of current identifier
        totype -> needed identifier type
        identifier -> protein identifier

    Returns:
        protein identifier of type specified in totype
    """
    base = 'http://www.uniprot.org'
    tool = 'mapping'
    params = {'from':fromtype,
    'to':totype,
    'format':'tab',
    'query':identifier,
    }
    #urllib turns the dictionary params into an encoded url suffix
    data = urllib.urlencode(params)
    #construct the UniProt URL
    url = base+'/'+tool+'?'+data
    #and grab the mapping
    response = urllib2.urlopen(url)
    #response.read() provides tab-delimited output of the mapping
    ConvertedId = response.read().split("\n")[1]
    return ConvertedId.split("\t")[1]
