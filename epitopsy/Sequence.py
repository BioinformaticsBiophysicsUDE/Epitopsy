import re
import urllib, urllib2
import collections
import numpy as np
from Bio import AlignIO, SeqIO, Seq, SeqRecord, pairwise2
from Bio.SubsMat import MatrixInfo as matlist

alignment_alphabet = {'-':1, 'A':2, 'C':3, 'D':4, 'E':5, 'F':6, 'G':7, 'H':8,
                      'I':9, 'K':10, 'L':11, 'M':12, 'N':13, 'P':14, 'Q':15,
                      'R':16, 'S':17, 'T':18, 'V':19, 'W':20, 'Y':21}

def make_alignment_upper_case(orig_aln):
    '''
    Transform all characters to uppercase.

    :param orig_aln: Biopython alignment
    :type  orig_aln: :class:`Bio.Align.MultipleSeqAlignment`

    :returns: Uppercase Biopython alignment, for
       :class:`Bio.Align.MultipleSeqAlignment`.
    :rtype: list(:class:`Bio.SeqRecord.SeqRecord`)

    Example::

        >>> alignment = AlignIO.read(open('seq.txt', 'rU'), 'fasta')
        >>> print alignment
        SingleLetterAlphabet() alignment with 3 rows and 178 columns
        iigp--gr-gfgkrrhpkkltplaykqfipnvaekt...sgg 3M1N:A|PDBID|CHAIN|SEQUENCE
        iigp--grpgfgkrrhpkkltplaykqfipnvaekt...sgg PRO1:A|NAME1|CHAIN|SEQUENCE
        iigpxxgrcgfgkrrhpkkltplaykqfipnvaekt...sgg PRO2:A|NAME2|CHAIN|SEQUENCE
        >>> new_alignment = Sequence.make_alignment_upper_case(alignment)
        >>> print Align.MultipleSeqAlignment(new_alignment)
        Alphabet() alignment with 3 rows and 178 columns
        IIGP--GR-GFGKRRHPKKLTPLAYKQFIPNVAEKT...SGG 3M1N:A|PDBID|CHAIN|SEQUENCE
        IIGP--GRPGFGKRRHPKKLTPLAYKQFIPNVAEKT...SGG PRO1:A|NAME1|CHAIN|SEQUENCE
        IIGPXXGRCGFGKRRHPKKLTPLAYKQFIPNVAEKT...SGG PRO2:A|NAME2|CHAIN|SEQUENCE

    '''
    aln = []
    for record in orig_aln:
        aln.append(SeqRecord.SeqRecord(Seq.Seq(record.seq.tostring().upper()),
                                   record.id, description=record.description))
    return aln

def make_alignment_clean(orig_aln):
    '''
    Delete sequences with non standard amino acids and reduce alignment size
    if unnecessary gaps are found (may happen when sequences are deleted).

    :param orig_aln: Biopython alignment
    :type  orig_aln: :class:`Bio.Align.MultipleSeqAlignment`

    :returns: Clean Biopython alignment, for
       :class:`Bio.Align.MultipleSeqAlignment`.
    :rtype: list(:class:`Bio.SeqRecord.SeqRecord`)

    Example::

        >>> alignment = AlignIO.read(open('seq.txt', 'rU'), 'fasta')
        >>> print alignment
        SingleLetterAlphabet() alignment with 3 rows and 178 columns
        IIGP--GR-GFGKRRHPKKLTPLAYKQFIPNVAEKT...SGG 3M1N:A|PDBID|CHAIN|SEQUENCE
        IIGP--GRPGFGKRRHPKKLTPLAYKQFIPNVAEKT...SGG PRO1:A|NAME1|CHAIN|SEQUENCE
        IIGPXXGRCGFGKRRHPKKLTPLAYKQFIPNVAEKT...SGG PRO2:A|NAME2|CHAIN|SEQUENCE
        >>> clean_alignment = Sequence.make_alignment_clean(alignment)
        input size: 3 x 178 (#seq x len(seq))
        output size: 2 x 176
        >>> print Align.MultipleSeqAlignment(clean_alignment)
        Alphabet() alignment with 2 rows and 176 columns
        IIGPGR-GFGKRRHPKKLTPLAYKQFIPNVAEKTLG...SGG 3M1N:A|PDBID|CHAIN|SEQUENCE
        IIGPGRPGFGKRRHPKKLTPLAYKQFIPNVAEKTLG...SGG PRO1:A|NAME1|CHAIN|SEQUENCE
        >>> # 3rd sequence with 'X' was deleted; two leading gaps were removed

    ..
        >>> for rec in clean_alignment:
        ...     print "%s...%s %s" % (rec.seq[:44], rec.seq[-3:], rec.id)

    '''
    #### acceptable amino acids
    aa_list = ['-', 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M',
               'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

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
    Open a stockholm formatted alignment file located at **input_path** and
    write the corresponding fasta formatted alignment to **output_path**.

    :param input_path: path to the stockholm formatted alignment
    :type  input_path: str
    :param output_path: path to the new fasta formatted alignment
    :type  output_path: str

    Example::

        >>> Sequence.format_stockholm_to_fasta('stock.txt', 'fasta.txt')

    '''
    aln = AlignIO.read(input_path, 'stockholm')
    with open(output_path, "w") as f:
        AlignIO.write(aln, f, 'fasta')


def get_pairwise_alignment(ref_seq, query_seq):
    '''
    Match the query sequence **query_seq** with the reference sequence
    **ref_seq** and return the corresponding aligned sequences in a
    dictionary with keys ('ref_seq', 'query_seq'). It uses the first
    calculated alignment; equally scored alignments are not considered.

    :param ref_seq: string of the reference sequence
    :type  ref_seq: str
    :param query_seq: string of the query sequence
    :type  query_seq: str

    :returns: Aligned sequences
    :rtype: dict

    Example::

        >>> align = Sequence.get_pairwise_alignment('IIGPTTGRGFGKRR',
        ...                                         'IIGPGGGRGFGKRR')
        >>> for seq in align.values():
        ...     print seq
        IIGPTTG--RGFGKRR
        IIGP--GGGRGFGKRR

    '''
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
    
    Return a list of the positions in the aligned sequence where the amino
    acid with the highest occurence has a frequency larger than the cutoff
    ('-' is also counted as an amino acid). The column index starts with 1!

    :param aln: Biopython alignment
    :type  aln: :class:`Bio.Align.MultipleSeqAlignment`
    :param cutoff: minimal occurence
    :type  cutoff: float

    :returns: list with the columns index
    :rtype: list(int)

    Example::

        >>> alignment = AlignIO.read(open('fasta.txt', 'rU'), 'fasta')
        >>> c = Sequence.get_almost_conserved_columns(alignment)
        >>> print c
        [1, 3, 4, 7, 8, 10, 11, 12, 14, 16]
        >>> # Result visualization
        >>> print alignment; print ''.join('+' if x in c else '-' for x
        ...  in range(1,alignment.get_alignment_length()+1)) + ' (at 95+%)'
        SingleLetterAlphabet() alignment with 3 rows and 16 columns
        IKGP--GR-GFGPRRH prot1
        IIGPRGGRPGFGKRVH prot2
        IIGPGGGRCGFGKRVH prot3
        +-++--++-+++-+-+ (at 95+%)

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
    '''
    Take an identifier and types of identifier (to and from) and call the
    UniProt mapping service. To get the right identifiers, visit
    http://www.uniprot.org/faq/28#conversion

    :param fromtype: type of current identifier
    :type  fromtype: str
    :param totype: needed identifier type
    :type  totype: str
    :param identifier: protein identifier
    :type  identifier: str

    :returns: Protein identifier of type **totype**.
    :rtype: str
    '''
    base = 'http://www.uniprot.org'
    tool = 'mapping'
    params = {'from': fromtype, 'to': totype,
              'format': 'tab', 'query': identifier}
    #urllib turns the dictionary params into an encoded url suffix
    data = urllib.urlencode(params)
    #construct the UniProt URL
    url = base+'/'+tool+'?'+data
    #and grab the mapping
    response = urllib2.urlopen(url)
    #response.read() provides tab-delimited output of the mapping
    ConvertedId = response.read().split("\n")[1]
    return ConvertedId.split("\t")[1]
