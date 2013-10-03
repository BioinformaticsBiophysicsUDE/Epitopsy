"""
@author: Christoph Wilms

This code is based on a matlab code provided by the authors of:
%     F Morcos, A Pagnani, B Lunt, A Bertolino, DS Marks, C Sander, 
%     R Zecchina, JN Onuchic, T Hwa, M Weigt (2011), Direct-coupling
%     analysis of residue co-evolution captures native contacts across 
%     many protein families, Proc. Natl. Acad. Sci. 108:E1293-1301.
%
% Permission is granted for anyone to copy, use, or modify this
% software and accompanying documents for any uncommercial
% purposes, provided this copyright notice is retained, and note is
% made of any changes that have been made. This software and
% documents are distributed without any warranty, express or
% implied. All use is entirely at the user's own risk.
%
% Any publication resulting from applications of DCA should cite:
%
%     F Morcos, A Pagnani, B Lunt, A Bertolino, DS Marks, C Sander, 
%     R Zecchina, JN Onuchic, T Hwa, M Weigt (2011), Direct-coupling
%     analysis of residue co-evolution captures native contacts across 
%     many protein families, Proc. Natl. Acad. Sci. 108:E1293-1301.

The code was ported to cython/python and the re-weighting algorithm was
optimized to use a quasi pairwise alignment for the determination of the 
pairwise sequence identities.
"""
import numpy as np

# "cimport" is used to import special compile-time information
# about the numpy module (this is stored in a file numpy.pxd which is
# currently part of the Cython distribution).
cimport numpy as np
cimport cython

from cython.parallel import prange, parallel

from libc.math cimport exp, log, round

from numpy cimport NPY_INT8 as NPY_int8
from numpy cimport NPY_INT16 as NPY_int16
from numpy cimport NPY_INT32 as NPY_int32
from numpy cimport NPY_INT64 as NPY_int64
from numpy cimport NPY_FLOAT16 as NPY_float16
from numpy cimport NPY_FLOAT32 as NPY_float32
from numpy cimport NPY_FLOAT64 as NPY_float64

int8 = np.dtype(np.int8)
int16 = np.dtype(np.int16)
int32 = np.dtype(np.int32)
int64 = np.dtype(np.int64)
float16 = np.dtype(np.float16)
float32 = np.dtype(np.float32)
float64 = np.dtype(np.float64)

## special imports
import time
import sys
import gzip
import multiprocessing

from Bio import AlignIO


cdef dict _get_alignment_alphabet():
    alignment_alphabet = {'-':1, 'A':2, 'C':3, 'D':4, 'E':5, 'F':6, 'G':7,
                          'H':8, 'I':9, 'K':10, 'L':11, 'M':12, 'N':13,
                          'P':14, 'Q':15, 'R':16, 'S':17, 'T':18, 'V':19,
                          'W':20, 'Y':21}
    return alignment_alphabet


def get_alignment_alphabet():
    return _get_alignment_alphabet()


def dca_new(inputfile, outputfile, identity_threshold= 0.8,
            pseudocount = 1.,
            num_threads=None):
    '''
    Direct Coupling Analysis (DCA)

    Args:
        inputfile -> file with fasta formatted sequences
        outputfile -> file containing N(N-1)/2 rows and 4 columns:
            residue i (column 1)
            residue j (column 2)
            MI(i,j) (Mutual Information between i and j)
            DI(i,j) (Direct Information between i and j)
        pseudocount -> lambda = pseudocount * Meff
        identity_treshold -> threshold for sequence id in reweighting

    F Morcos, A Pagnani, B Lunt, A Bertolino, DS Marks, C Sander,
    R Zecchina, JN Onuchic, T Hwa, M Weigt (2011), Direct-coupling
    analysis of residue co-evolution captures native contacts across
    many protein families, Proc. Natl. Acad. Sci. 108:E1293-1301.
    '''
    start_time = time.time()
    # settings
    dca_version = 'new' # -> dca identifier = 1
    dca_identifier = 1
    if num_threads is None:
        num_threads = multiprocessing.cpu_count()
    elif num_threads > multiprocessing.cpu_count():
        num_threads = multiprocessing.cpu_count()

    # run dca
    run_dca(inputfile,
            outputfile,
            identity_threshold,
            pseudocount,
            dca_version,
            dca_identifier,
            num_threads)

    end_time = time.time()
    print('finished after {0}s ...'.format(end_time-start_time))


def dca_old(inputfile, outputfile, identity_threshold= 0.8,
            pseudocount = 1.,
            num_threads=None):
    '''
    Direct Coupling Analysis (DCA)

    Args:
        inputfile -> file with fasta formatted sequences
        outputfile -> file containing N(N-1)/2 rows and 4 columns:
            residue i (column 1)
            residue j (column 2)
            MI(i,j) (Mutual Information between i and j)
            DI(i,j) (Direct Information between i and j)
        pseudocount_weight -> relative weight of pseudo count
            pseudocount = lambda / (Meff + lambda)
        theta -> threshold for sequence id in reweighting

    F Morcos, A Pagnani, B Lunt, A Bertolino, DS Marks, C Sander,
    R Zecchina, JN Onuchic, T Hwa, M Weigt (2011), Direct-coupling
    analysis of residue co-evolution captures native contacts across
    many protein families, Proc. Natl. Acad. Sci. 108:E1293-1301.
    '''
    start_time = time.time()
    # settings
    print('read alignment ...')
    [N, M, q, Z] = return_alignment(inputfile)
    print('compute frequencies ...')
    dca_version = 'old' # -> dca identifier = 0
    dca_identifier = 0
    if num_threads is None:
        num_threads = multiprocessing.cpu_count()
    elif num_threads > multiprocessing.cpu_count():
        num_threads = multiprocessing.cpu_count()

    # run dca
    run_dca(inputfile,
            outputfile,
            identity_threshold,
            pseudocount,
            dca_version,
            dca_identifier,
            num_threads)

    end_time = time.time()
    print('finished after {0}s ...'.format(end_time-start_time))


cdef void run_dca(inputfile, outputfile, identity_threshold,
            pseudocount,
            dca_version,
            dca_identifier,
            num_threads):
    '''
    This functions implements the calculation of DCA.
    '''
    # read alignment
    print('read alignment ...')
    [N, M, q, Z] = return_alignment(inputfile)

    ## settings
    cdef double n_pairs = 0.5 * ((N * N) - N)
    cdef double[:,:] results = np.zeros([n_pairs, 2], dtype=np.float64) # mi di

    # compute frequencies
    print('compute frequencies ...')
    # compute the weight for each count based on the identity to all
    # other sequences
    W = Compute_Meff(Z, M, N, identity_threshold, num_threads, dca_identifier)
    # use these individual counts to calculate single and pair site
    # frequencies
    [counts_ij_true, counts_i_true, Meff] = Compute_True_Frequencies(Z,M,N,q, W)

    # release memory
    del(Z)
    del(W)

    header_string = "#### {0} DCA N={1} M={2} Meff={3} q={4} pseudocount={5} x={6}".format(dca_version,N,M,Meff,q, pseudocount, identity_threshold)
    print(header_string)

    # calculate MI
    print('compute MI ...')
    results = calculate_mi(results, counts_i_true, counts_ij_true, Meff, N, q)

    # add a pseudo count
    print("add pseudocount ...")
    [Pij, Pi] = with_pc(counts_ij_true,counts_i_true,pseudocount,N,q,Meff)

    # realease memory
    del(counts_i_true)
    del(counts_ij_true)

    # calculate the correlation matrix
    print('compute C ...')
    C = Compute_C(Pij,Pi,N,q)

    # invert the matrix
    print('invert C ...')
    invC = np.linalg.inv(C)

    # realease memory
    del(C)

    # start the DI calculation
    print('compute DI ...')
    results = Compute_Results(Pij, Pi, invC, N, q, results)


    cdef int i,j
    cdef double MI,DI
    cdef long index_counter = 0
    output_list = []
    for i in range(N-1):
        for j in range(i+1,N):
            MI = results[index_counter,0]
            DI = results[index_counter,1]
            index_counter += 1
            output_list.append("{0} {1} {2} {3}\n".format(i+1,j+1, MI, DI))

    # save data
    if outputfile.endswith('.gz'):
        with gzip.open(outputfile, 'w') as f:
            f.write(header_string + "\n")
            f.write("i j mi di\n")
            for item in output_list:
                f.write(item)

    else:
        with open(outputfile, 'w') as f:
            f.write(header_string + "\n")
            f.write("i j mi di\n")
            for item in output_list:
                f.write(item)


def mi(inputfile, outputfile):
    '''
    Mutual Information (mi), calls dca_new with theta = 0

    Args:
        inputfile -> file with fasta formatted sequences
        outputfile -> file containing N(N-1)/2 rows and 3 columns:
            residue i (column 1)
            residue j (column 2)
            MI(i,j) (Mutual Information between i and j)
            DI(i,j) (Direct Information between i and j)
    '''
    dca_new(inputfile, outputfile, identity_threshold=1.)


def read_alignment_and_calc_Meff(inputfile, identity_threshold, num_threads=None):
    '''
    '''
    start_time = time.time()
    print('read alignment ...')
    [N, M, q, Z] = return_alignment(inputfile)
    print('compute frequencies ...')
    dca_version = 'new' # -> dca identifier = 1
    dca_identifier = 1
    if num_threads is None:
        num_threads = multiprocessing.cpu_count()
    elif num_threads > multiprocessing.cpu_count():
        num_threads = multiprocessing.cpu_count()
    W = Compute_Meff(Z, M, N, identity_threshold, num_threads, dca_identifier)
    end_time = time.time()
    print('finished after {0}s ...'.format(end_time-start_time))
    return [N,M,q,Z,W]

def get_entropies(inputfile, identity_threshold = 1., num_threads=None):
    '''
    Args:
        inputfile -> path to the alignment
        identity_threshold -> default is 1
        num_threads -> default is None, i.e. all available processors
    Returns:
        A list with the entropies.
    '''
    cdef int M,N
    dca_version = 'new' # -> dca identifier = 1
    dca_identifier = 1
    if num_threads is None:
        num_threads = multiprocessing.cpu_count()
    elif num_threads > multiprocessing.cpu_count():
        num_threads = multiprocessing.cpu_count()
    start_time = time.time()
    #print('read alignment ...')
    [N, M, q, Z] = return_alignment(inputfile)
    #print("Compute Meff")
    W = Compute_Meff(Z, M, N, identity_threshold, num_threads, dca_identifier)
    #print('compute entropies ...')
    entropy = Compute_Entropies(Z, M, N, q, W)
    return entropy


cdef long map_upper_matrix(int i, int j, int M) nogil:
    '''
    Map indices from upper triuangular matrix (without diagonal to
    one dimensional vector.
    i starts with 0 ... M-1
    j starts with 0 ... M-1
    '''
    cdef double ip = <double> i + 1
    cdef double jp = <double> j + 1
    cdef double Mp = <double> M
    cdef double score = (ip - 1.) * Mp - ((ip-1.)**2 + (ip-1.))*0.5+jp-ip
    return <long> (score - 1)

def get_upper_seq_identity_matrix(inputfile):
    '''
    Args:
        inputfile -> path to the alignment

    Returns:
        An 1d array with 100*seq_identity (i.e. 0.34) = 34.
    '''
    start_time = time.time()
    print('read alignment ...')
    cdef long M,N
    cdef double q
    [N, M, q, Z] = return_alignment(inputfile)
    print('compute frequencies ...')
    dca_version = 'new' # -> dca identifier = 1
    dca_identifier = 1
    cdef int i,j,k
    cdef double vec_len = 0.5 * ((M * M) - M)
    cdef np.ndarray[np.int8_t,ndim=1] identity_vec_np = np.zeros(vec_len,dtype=np.int8)
    cdef signed char[:] identity_vec = identity_vec_np
    cdef double identity_score
    cdef double accuracy = 100
    cdef signed char[:,:] Z_view = Z
    cdef long counter = 0
    cdef long vec_ind
    print("vec len: {0}".format(vec_len))
    try:
        with nogil:
            with parallel():
                for i in prange(M - 1):
                    for j in range(i + 1, M):
                        identity_score = get_sequence_identity(Z_view, i, j, N)
                        vec_ind = map_upper_matrix(i,j,M)
                        identity_vec[vec_ind] = <signed char> round(identity_score*accuracy)
    except:
        print("error")
        print(i,j,vec_ind, vec_len)

    end_time = time.time()
    print('finished after {0}s ...'.format(end_time-start_time))
    return identity_vec_np


def get_upper_aln_identity_matrix(inputfile):
    '''
    Args:
        inputfile -> path to the alignment

    Returns:
        An 1d array with 100*seq_identity (i.e. 0.34) = 34.
    '''
    start_time = time.time()
    print('read alignment ...')
    cdef long M,N
    cdef double q
    [N, M, q, Z] = return_alignment(inputfile)
    print('compute frequencies ...')
    dca_version = 'new' # -> dca identifier = 1
    dca_identifier = 1
    cdef int i,j,k
    cdef double vec_len = 0.5 * ((M * M) - M)
    cdef np.ndarray[np.int8_t,ndim=1] identity_vec_np = np.zeros(vec_len,dtype=np.int8)
    cdef signed char[:] identity_vec = identity_vec_np
    cdef double identity_score
    cdef double accuracy = 100
    cdef signed char[:,:] Z_view = Z
    cdef long counter = 0
    cdef long vec_ind
    print("vec len: {0}".format(vec_len))
    try:
        with nogil:
            with parallel():
                for i in prange(M - 1):
                    for j in range(i + 1, M):
                        identity_score = get_alignment_identity(Z_view, i, j, N)
                        vec_ind = map_upper_matrix(i,j,M)
                        identity_vec[vec_ind] = <signed char> round(identity_score*accuracy)
    except:
        print("error")
        print(i,j,vec_ind, vec_len)

    end_time = time.time()
    print('finished after {0}s ...'.format(end_time-start_time))
    return identity_vec_np


def return_alignment(inputfile):
    if inputfile.endswith('.fa'):
        alignment = AlignIO.read(inputfile, "fasta")
    elif inputfile.endswith('.fa.gz'):
        f = gzip.open(inputfile)
        alignment = AlignIO.read(f, "fasta")
        f.close()

    # num of sequences
    cdef int M = len(alignment)
    # num of alignment positions
    cdef int N = len(alignment[0])
    cdef int q = 21
    aa_dict = _get_alignment_alphabet()
    cdef np.ndarray[np.int8_t, ndim=2] Z_np = np.zeros([M,N],dtype=np.int8)
    cdef signed char[:,:] Z = Z_np

    cdef int i,j
#    cdef int score
    cdef signed char score

    for i in range(M):
        record = alignment[i]
        sequence = record.seq.tostring()
        for j in range(N):
            aa_char = sequence[j]
            if aa_char in aa_dict:
                score = aa_dict[aa_char]
            else:
                score = 1
            Z[i,j] = score

    return [N,M,q,Z_np]

# This version should work fine, if one changes something, please comment these
# options
@cython.boundscheck(False)
@cython.cdivision(True)
cdef np.ndarray Compute_Entropies(np.ndarray[np.int8_t, ndim=2] Z, int M,
        int N, int q, np.ndarray[np.float64_t,ndim=1] W):
    '''
    Compute the real unmodified frequencies for the MI calculation.
    '''
    cdef int i,j,l,k,alpha,beta
    cdef double Meff = np.sum(W)
    cdef np.ndarray[np.float64_t, ndim=2] Pi_true_np = np.zeros([N,q],dtype=np.float64)
    cdef np.ndarray[np.float64_t, ndim=1] entropy_np = np.zeros(N,dtype=np.float64)
    cdef double [:,:] Pi_true = Pi_true_np
    cdef double [:] entropy = entropy_np
    cdef int z_ind

    for j in range(M):
        for i in range(N):
            z_ind = (Z[j,i] - 1)
            Pi_true[i, z_ind] = Pi_true[i, z_ind] + W[j]

    # make frequencies or probabilities
    for i in range(N):
        for alpha in range(q):
            Pi_true[i,alpha] = Pi_true[i,alpha] / Meff

    # make entropies
    cdef double s = 0
    cdef double freq = 0
    cdef double q_float = <double> q
    for i in range(N):
        s = 0
        for alpha in range(q):
            freq = Pi_true[i, alpha]
            if freq != 0:
                s = s - freq * log(freq)/log(q_float)
        entropy[i] = s

    return entropy_np

@cython.cdivision(True)
@cython.boundscheck(False)
# cdivision leads to small deviations in the calculations
cdef double get_alignment_identity(signed char[:,:] Z, int i, int j, int N) nogil:
    '''
    Does not check the length of both objects.
    Args:
        ref -> iterable object which will be compared with query
        query -> iterable object which will be compared with ref

    Returns:
        Alignment similarity [0,1].
    '''
    cdef double counter = 0.
    cdef double n_items = <double> N
    cdef int k
    cdef double score = 0.
    cdef double return_value = 0.

    for k in range(N):
        if Z[i,k] == Z[j,k]:
            counter += 1.

    score = counter / n_items

    return score

@cython.cdivision(True)
@cython.boundscheck(False)
# cdivision leads to small deviations in the calculations
cdef double get_sequence_identity(signed char[:,:] Z, int i, int j, int N) nogil:
    '''
    Does not check the length of both objects.
    Args:
        ref -> iterable object which will be compared with query
        query -> iterable object which will be compared with ref

    Returns:
        Alignment similarity [0,1].
    '''
    cdef double counter = 0.
    cdef double n_items = <double> N
    cdef int k
    cdef double score = 0.
    cdef double return_value = 0.

    for k in range(N):
        if Z[i,k] == Z[j,k] == 1:
            n_items -= 1.
            continue
        else:
            if Z[i,k] == Z[j,k]:
                counter += 1.

    score = counter / n_items

    return score

@cython.cdivision(True)
@cython.boundscheck(False)
cdef double get_seq_weight(signed char[:,:] Z_view, int i, int k, int N,
                          double identity_threshold, int dca_identifier) nogil:
    cdef double identity_score
    cdef double return_value
    if dca_identifier == 1:
        identity_score = get_sequence_identity(Z_view, i, k, N)
    elif dca_identifier == 0:
        identity_score = get_alignment_identity(Z_view, i, k, N)

    if identity_score < identity_threshold:
        return_value = 0.
    else:
        return_value = 1.

    return return_value

# This version should work fine, if one changes something, please comment these
# options
@cython.cdivision(True)
@cython.boundscheck(False)
cdef np.ndarray Compute_Meff(np.ndarray[np.int8_t, ndim=2] Z, int M, int N,
        double identity_threshold, int num_threads, int dca_identifier):
    cdef np.ndarray[np.float64_t, ndim = 1] W
    cdef signed char[:, :] Z_view = Z
    cdef np.ndarray[np.float64_t, ndim = 1] hamming_sum_np = np.zeros(M,dtype=np.float64)
    cdef np.ndarray[np.float64_t, ndim = 1] util_sum_np = np.zeros(M,dtype=np.float64)
    cdef double[:] util_sum = util_sum_np
    cdef double[:] hamming_sum = hamming_sum_np
    cdef int i, j, l, k
    cdef double identity_score

    if (identity_threshold < 1.):
        # W = (1./(1.+np.sum(distance.squareform((distance.pdist(Z,'hamming') < theta).astype(double)),0)))
#        for i in range(M):
        with nogil:
            for i in range(M - 1):
                # set util sum to zeros
                for l in range(M):
                    util_sum[l] = 0.

                with parallel(num_threads=num_threads):
                    for k in prange(i + 1, M):
                        util_sum[k] = get_seq_weight(Z_view, i, k, N, identity_threshold, dca_identifier)

                for j in range(i + 1, M):
                    hamming_sum[j] += util_sum[j]
                    hamming_sum[i] += util_sum[j]

        W = (1. / (1. + hamming_sum_np))
    else:
        W = np.ones(M,dtype=np.float64)

    return W

# This version should work fine, if one changes something, please comment these
# options
@cython.boundscheck(False)
cdef list Compute_True_Frequencies(np.ndarray[np.int8_t, ndim=2] Z, int M,
        int N, int q,
        np.ndarray[np.float64_t, ndim=1] W):
    cdef int i,j,l,k,alpha,beta
    cdef double Meff = np.sum(W)
    cdef np.ndarray[np.float64_t, ndim=4] counts_ij_true_np = np.zeros([N,N,q,q],dtype=np.float64)
    cdef double [:,:,:,:] counts_ij_true = counts_ij_true_np
    cdef np.ndarray[np.float64_t, ndim=2] counts_i_true_np = np.zeros([N,q],dtype=np.float64)
    cdef double [:,:] counts_i_true = counts_i_true_np
    cdef int z_ind

    for j in range(M):
        for i in range(N):
            z_ind = (Z[j,i] - 1)
            counts_i_true[i, z_ind] = counts_i_true[i, z_ind] + W[j]


    cdef int z_ind_i, z_ind_j

    for l in range(M):
        for i in range(N-1):
            z_ind_i = Z[l,i] - 1 # python index
            for j in range(i+1,N):
                z_ind_j = Z[l,j] - 1 # python index
                counts_ij_true[i,j,z_ind_i, z_ind_j] += W[l]
                counts_ij_true[j,i,z_ind_j, z_ind_i] = counts_ij_true[i,j,z_ind_i, z_ind_j]

    cdef np.ndarray[np.float64_t, ndim=2] scra_np = np.eye(q,dtype=np.float64)
    cdef double[:,:] scra = scra_np
    for i in range(N):
        for alpha in range(q):
            for beta in range(q):
                counts_ij_true[i,i,alpha,beta] = counts_i_true[i,alpha] * scra[alpha,beta]

    return [counts_ij_true_np, counts_i_true_np, Meff]

# This version should work fine, if one changes something, please comment these
# options
@cython.boundscheck(False)
@cython.cdivision(True)
cdef list with_pc(np.ndarray[np.float64_t, ndim=4] counts_ij_true_np,
                  np.ndarray[np.float64_t, ndim=2] counts_i_true_np,
                  double pseudocount, int N, int q, double Meff):
    cdef int i, j, alpha, beta
    cdef double q_float = <double> q
    cdef double lambda_value = pseudocount * Meff
    cdef double norm_factor = (1. / (lambda_value + Meff))

    cdef double [:,:] counts_i_true = counts_i_true_np
    cdef double [:,:,:,:] counts_ij_true = counts_ij_true_np
    cdef np.ndarray[np.float64_t, ndim=4] Pij_np = np.zeros([N,N,q,q],dtype=np.float64)
    cdef double [:,:,:,:] Pij = Pij_np
    cdef np.ndarray[np.float64_t, ndim=2] Pi_np = np.zeros([N,q],dtype=np.float64)
    cdef double [:,:] Pi = Pi_np
    cdef np.ndarray[np.float64_t, ndim=2] scra_np = np.eye(q,dtype=np.float64)
    cdef double [:,:] scra = scra_np

    for i in range(N):
        for alpha in range(q):
            Pi[i,alpha] = norm_factor * ((lambda_value/q_float) + counts_i_true[i,alpha])
            for j in range(N):
                for beta in range(q):
                    Pij[i,j,alpha,beta] = norm_factor * ((lambda_value/(q_float*q_float)) + counts_ij_true[i,j,alpha,beta])

    for i in range(N):
        for alpha in range(q):
            for beta in range(q):
                Pij[i,i,alpha,beta] = norm_factor * ((lambda_value/q_float) * scra[alpha,beta] + counts_ij_true[i,i,alpha,beta])

    return [Pij_np, Pi_np]


# This version should work fine, if one changes something, please comment these
# options
@cython.boundscheck(False)
cdef int mapkey(int i, int alpha, int q):
    return ((q-1)*(i-1) + alpha)

# This version should work fine, if one changes something, please comment these
# options
@cython.boundscheck(False)
cdef np.ndarray Compute_C(np.ndarray[np.float64_t, ndim=4] Pij_np,
        np.ndarray[np.float64_t, ndim=2] Pi_np, int N, int q):
    cdef int qm = q-1
    cdef np.ndarray[np.float64_t, ndim=2] C_np = np.zeros([N*(q-1),N*qm],dtype=np.float64)

    # memory views for efficient indexing
    cdef double[:,:] C = C_np
    cdef double [:,:,:,:] Pij = Pij_np
    cdef double [:,:] Pi = Pi_np

    cdef int i,j,alpha,beta

    for i in range(N):
        for j in range(N):
            for alpha in range(qm):
                for beta in range(qm):
                    C[mapkey(i,alpha,q),mapkey(j,beta,q)] = (Pij[i,j,alpha,beta]
                            - Pi[i,alpha] * Pi[j, beta])

    return C_np


cdef double[:,:] Compute_Results(np.ndarray[np.float64_t, ndim=4] Pij,
        np.ndarray[np.float64_t, ndim=2] Pi,
        np.ndarray[np.float64_t, ndim=2] invC,
        int N,
        int q,
        double[:,:] results):

    cdef int i,j
    cdef np.ndarray[np.float64_t, ndim=2] W_mf
    cdef double DI_mf_pc
    cdef long index_counter = 0

    for i in range(N-1):
        for j in range(i+1, N):
            ## direct information
            W_mf = ReturnW(invC, i, j, q)
            DI_mf_pc = bp_link(i,j,W_mf,Pi,q)
            results[index_counter,1] = DI_mf_pc
            index_counter += 1

    return results

cdef np.ndarray ReturnW(np.ndarray[np.float64_t, ndim=2] invC_np, int i,
        int j, int q):
    cdef np.ndarray[np.float64_t, ndim=2] W_np = np.ones([q,q],dtype=np.float64)
    cdef double[:,:] W = W_np
    cdef double[:,:] invC = invC_np
    cdef int qm = q-1
    cdef int qmi = qm * (i-1)
    cdef int qmj = qm * (j-1)
    cdef int qmia

    cdef int alpha,beta

    for alpha in range(qm):
        qmia = qmi + alpha
        for beta in range(qm):
            W[alpha, beta] = exp(-invC[qmia, qmj+beta])
    return W_np

cdef double bp_link(int i,int j, np.ndarray[np.float64_t, ndim=2] W_mf,
        np.ndarray[np.float64_t, ndim=2] Pi, int q):
    cdef np.ndarray[np.float64_t, ndim=1] mu1
    cdef np.ndarray[np.float64_t, ndim=1] mu2
    cdef double DI
    [mu1, mu2] = compute_mu(i,j,W_mf,Pi,q)
    DI = compute_di(i,j,W_mf, mu1,mu2,Pi)
    return DI

cdef list compute_mu(int i, int j, np.ndarray[np.float64_t, ndim=2] W,
        np.ndarray[np.float64_t, ndim=2] Pi, int q):
    cdef double epsilon = 1e-4
    cdef double diff = 1.0
    cdef double q_float = <double> q
    cdef np.ndarray[np.float64_t, ndim=1] mu1 = np.ones(q,dtype=np.float64)/q_float
    cdef np.ndarray[np.float64_t, ndim=1] mu2 = np.ones(q,dtype=np.float64)/q_float
    cdef np.ndarray[np.float64_t, ndim=1] pi = Pi[i,:]
    cdef np.ndarray[np.float64_t, ndim=1] pj = Pi[j,:]
    cdef np.ndarray[np.float64_t, ndim=1] new1
    cdef np.ndarray[np.float64_t, ndim=1] new2
    cdef np.ndarray[np.float64_t, ndim=1] scra1
    cdef np.ndarray[np.float64_t, ndim=1] scra2

    while diff > epsilon:
        scra1 = np.dot(mu2, np.transpose(W))
        scra2 = np.dot(mu1, W)

        new1 = pi / scra1
        new1 = new1 / np.sum(new1)

        new2 = pj / scra2
        new2 = new2 / np.sum(new2)

        diff = np.max([np.abs(new1-mu1), np.abs(new2-mu2)])

        mu1 = new1
        mu2 = new2

    return [mu1, mu2]

cdef double compute_di(int i, int j, np.ndarray[np.float64_t, ndim=2] W,
        np.ndarray[np.float64_t, ndim=1] mu1,
        np.ndarray[np.float64_t, ndim=1] mu2,
        np.ndarray[np.float64_t, ndim=2] Pi):
    cdef double DI
    tiny = 1.0e-100
    Pdir = W * np.dot(mu1.reshape(-1,1),mu2.reshape(1,-1))
    Pdir = Pdir / np.sum(Pdir)
    Pfac = np.dot(Pi[i,:].reshape(-1,1), Pi[j,:].reshape(1,-1))
    DI = np.trace( np.dot(Pdir.transpose(), np.log( (Pdir+tiny) / (Pfac+tiny))))
    return DI

@cython.boundscheck(False)
@cython.cdivision(True)
cdef double[:,:] calculate_mi(double[:,:] results,
        np.ndarray[np.float64_t,ndim=2] counts_i_true_np,
        np.ndarray[np.float64_t,ndim=4] counts_ij_true_np,
        double Meff,
        int N,
        int q):
    cdef int i,j,alpha,beta
    cdef double MI = 0.
    cdef double q_float = <double> q
    cdef long index_counter = 0

    # memory views
    cdef double [:,:,:,:] counts_ij_true = counts_ij_true_np
    cdef double [:,:] counts_i_true = counts_i_true_np

    # calculate MI
    cdef double f_i,f_j,f_ij
    for i in range(N-1):
        for j in range(i+1,N):
            MI = 0.
            for alpha in range(q):
                for beta in range(q):
                    if counts_ij_true[i,j,alpha,beta] > 0:
                        # normalize counts -> frequencies/probabilities
                        f_i = counts_i_true[i,alpha] / Meff
                        f_j = counts_i_true[j,beta] / Meff
                        f_ij = counts_ij_true[i,j,alpha,beta] / Meff
                        MI = MI + f_ij * log(f_ij / f_i / f_j) / log(q_float)

            results[index_counter,0] = MI
            index_counter += 1

    return results
