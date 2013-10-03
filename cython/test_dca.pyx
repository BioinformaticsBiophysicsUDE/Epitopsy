"""
@author: Christoph Wilms
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
from scipy.linalg import sqrtm


cdef dict _get_alignment_alphabet():
    alignment_alphabet = {'-':1, 'A':2, 'C':3, 'D':4, 'E':5, 'F':6, 'G':7,
                          'H':8, 'I':9, 'K':10, 'L':11, 'M':12, 'N':13,
                          'P':14, 'Q':15, 'R':16, 'S':17, 'T':18, 'V':19,
                          'W':20, 'Y':21}
    return alignment_alphabet


def get_alignment_alphabet():
    return _get_alignment_alphabet()


def hopf_new(inputfile, outputfile, identity_threshold=0.8,
            pseudocount=0.5,
            patterns=128,
            num_threads=None):
    '''
    Args:
        inputfile -> file with fasta formatted sequences
        outputfile -> file containing N(N-1)/2 rows and 4 columns:
            residue i (column 1)
            residue j (column 2)
            MI(i,j) (Mutual Information between i and j)
            DI(i,j) (Direct Information between i and j)
        pseudocount -> lambda = pseudocount * Meff
        identity_treshold -> threshold for sequence id in reweighting
        patterns -> ?

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
    run_hopf(inputfile,
            outputfile,
            identity_threshold,
            pseudocount,
            patterns,
            dca_version,
            dca_identifier,
            num_threads)

    end_time = time.time()
    print('finished after {0}s ...'.format(end_time-start_time))



cdef void run_hopf(char* inputfile,char* outputfile, float identity_threshold,
            pseudocount,
            int patterns,
            char* dca_version,
            int dca_identifier,
            int num_threads):
    '''
    This functions implements the calculation of DCA.
    '''
    # read alignment
    print('read alignment ...')
    cdef int N,M,i,j
    cdef float hp
    [N, M, q, Z] = return_alignment(inputfile)

    ## settings
    cdef double n_pairs = 0.5 * ((N * N) - N)
    cdef float[:,:] results = np.zeros([n_pairs, 2], dtype=np.float32) # mi di

    # compute frequencies
    print('compute frequencies ...')
    # compute the weight for each count based on the identity to all
    # other sequences
    W = Compute_Meff(Z, M, N, identity_threshold, num_threads, dca_identifier)
    # use these individual counts to calculate single and pair site
    # frequencies
    print('compute true frequencies ...')
    [Pij_true, Pi_true, Meff] = Compute_True_Frequencies(Z,M,N,q, W)

    # release memory
    del(Z)
    del(W)

    header_string = "#### {0} DCA float32 N={1} M={2} Meff={3} q={4} pseudocount={5} x={6}".format(dca_version,N,M,Meff,q, pseudocount, identity_threshold)
    print(header_string)

    # add a pseudo count
    print("add pseudocount ...")
    [Pij, Pi] = with_pc(Pij_true,Pi_true,pseudocount,N,q,Meff)

    # realease memory
    del(Pi_true)
    del(Pij_true)

    # calculate the correlation matrix
    print('compute C ...')
    [C,C_self] = Compute_C(Pij,Pi,N,q)
    print('compute D ...')
    D = sqrtm(C_self)
    print("###debug D")
    np.save("d.npy", D)
#    print("###debug1")
    print('compute gamma ...')
    # matlab / -> matlab_matrix_division
    # matlab \ -> np.linalg.solve()
    Gamma = matlab_matrix_division(np.linalg.solve(D,C),D)
    print('compute eigenvectors ...')
    [V, Lambda] = np.linalg.eig(Gamma)
    Vtilde = np.linalg.solve(D,V)

    print("potts part 2 ...")
    ## potts 2
    Lambda1 = np.eye(N*(q-1))
    ll = np.diag(Lambda1) - np.ones([N*(q-1),1]) - np.log(np.diag(Lambda))
    b_ascending = np.argsort(ll)
    b = []
    for item in reversed(b_ascending):
        b.append(item)
    b = np.array(b)

    for i in range(patterns):
        Lambda1[b[i],b[i]] = Lambda[b[i],b[i]]

#    print("###debug2")
    Lambda2 = matlab_matrix_division((Lambda1 - np.eye(N*(q-1))) , Lambda1)

    invC = -np.matrix(Vtilde) * np.matrix(Lambda2) * np.matrix(np.transpose(Vtilde))
    invC = np.array(invC)

    cdef np.ndarray[np.float32_t, ndim=2] F_apc_np
    F_apc_np = calc_norm_apc(invC, N, q)

    # save data
    if outputfile.endswith('.gz'):
        with gzip.open(outputfile, 'w') as f:
            f.write(header_string + "\n")
            f.write("i j mi\n")
            for i in range(N-1):
                for j in range(i+1,N):
                    hp = F_apc_np[i,j]
                    f.write("{0} {1} {2}\n".format(i,j,hp))

    else:
        with open(outputfile, 'w') as f:
            f.write(header_string + "\n")
            f.write("i j mi\n")
            for i in range(N-1):
                for j in range(i+1,N):
                    hp = F_apc_np[i,j]
                    f.write("{0} {1} {2}\n".format(i,j,hp))


def matlab_matrix_division(B, A):
    '''
    Solve:
        C = B/A

    In python:
       C = B/A          | *A
       C * A = B        | transpose
       A' * C' = B'
    this can be interpreted as a standard problem Ax=B
    '''
    return np.linalg.lstsq(A.T,B.T)[0].T

    
cdef np.ndarray calc_norm_apc(np.ndarray[np.float32_t, ndim=2] invC_np,
                              int N, int q):
    cdef np.ndarray[np.float32_t, ndim=2] F_np = np.zeros([N,N])
    cdef np.ndarray[np.float32_t, ndim=2] F_apc_np = np.zeros([N,N])
    cdef float[:,:] F = F_np
    cdef float[:,:] F_apc = F_apc_np
    cdef int i,j,a,b
    for i in range(N-1):
        for j in range(i+1,N):
            J_mf = ReturnW(invC_np,i,j,q)
            J_j = np.mean(J_mf,0)
            J_i = np.mean(J_mf,1)
            J = np.mean(J_i)
            for a in range(q):
                for b in range(q):
                    J_mf[a,b] = J_mf[a,b] - J_i[a] - J_j[b] + J

            F[i,j] = np.linalg.norm(J_mf,"fro")
            F[j,i] = F[i,j]

    F_i = np.mean(F,0)
    F_av = np.mean(F_i,0)
    for i in range(N-1):
        for j in range(i+1,N):
            F_apc[i,j] = F[i,j] - ((F_i[i]*F_i[j]) / F_av)
            F_apc[j,i] = F_apc[i,j]

    return F_apc_np

cdef list return_alignment(char* inputfile):
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


@cython.cdivision(True)
@cython.boundscheck(False)
# cdivision leads to small deviations in the calculations
cdef float get_alignment_identity(signed char[:,:] Z, int i, int j, int N) nogil:
    '''
    Does not check the length of both objects.
    Args:
        ref -> iterable object which will be compared with query
        query -> iterable object which will be compared with ref

    Returns:
        Alignment similarity [0,1].
    '''
    cdef float counter = 0.
    cdef float n_items = <float> N
    cdef int k
    cdef float score = 0.
    cdef float return_value = 0.

    for k in range(N):
        if Z[i,k] == Z[j,k]:
            counter += 1.

    score = counter / n_items

    return score

@cython.cdivision(True)
@cython.boundscheck(False)
# cdivision leads to small deviations in the calculations
cdef float get_sequence_identity(signed char[:,:] Z, int i, int j, int N) nogil:
    '''
    Does not check the length of both objects.
    Args:
        ref -> iterable object which will be compared with query
        query -> iterable object which will be compared with ref

    Returns:
        Alignment similarity [0,1].
    '''
    cdef float counter = 0.
    cdef float n_items = <float> N
    cdef int k
    cdef float score = 0.
    cdef float return_value = 0.

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
cdef float get_seq_weight(signed char[:,:] Z_view, int i, int k, int N,
                          float identity_threshold, int dca_identifier) nogil:
    cdef float identity_score
    cdef float return_value
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
        float identity_threshold, int num_threads, int dca_identifier):
    cdef np.ndarray[np.float32_t, ndim = 1] W
    cdef signed char[:, :] Z_view = Z
    cdef np.ndarray[np.float32_t, ndim = 1] hamming_sum_np = np.zeros(M,dtype=np.float32)
    cdef np.ndarray[np.float32_t, ndim = 1] util_sum_np = np.zeros(M,dtype=np.float32)
    cdef float[:] util_sum = util_sum_np
    cdef float[:] hamming_sum = hamming_sum_np
    cdef int i, j, l, k
    cdef float identity_score

    if (identity_threshold < 1.):
        # W = (1./(1.+np.sum(distance.squareform((distance.pdist(Z,'hamming') < theta).astype(float)),0)))
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
        W = np.ones(M,dtype=np.float32)

    return W

# This version should work fine, if one changes something, please comment these
# options
#@cython.boundscheck(False)
cdef list Compute_True_Frequencies(np.ndarray[np.int8_t, ndim=2] Z, int M,
        int N, int q,
        np.ndarray[np.float32_t, ndim=1] W):
    cdef int i,j,alpha,beta,l,z_ind, z_ind_i, z_ind_j
    cdef np.ndarray[np.float32_t,ndim=2] Pi_true = np.zeros([N,q],dtype=np.float32)
    cdef np.ndarray[np.float32_t,ndim=4] Pij_true = np.zeros([N,N,q,q],dtype=np.float32)
    cdef float Meff = np.sum(W)
    
    for j in range(M):
        for i in range(N):
            z_ind = Z[j,i] - 1
            Pi_true[i,z_ind] = Pi_true[i,z_ind] + W[j]

    Pi_true = Pi_true / Meff

    for l in range(M):
        for i in range(N):
            for j in range(i+1,N):
                z_ind_i = Z[l,i] - 1
                z_ind_j = Z[l,j] - 1
                Pij_true[i,j, z_ind_i, z_ind_j] = Pij_true[i,j,z_ind_i,z_ind_j] + W[j]
                Pij_true[j,i,z_ind_j,z_ind_i] = Pij_true[i,j,z_ind_i,z_ind_j]
    
    Pij_true = Pij_true / Meff

    cdef np.ndarray[np.float32_t,ndim=2] scra = np.eye(q,dtype=np.float32)
    for i in range(N):
        for alpha in range(q):
            for beta in range(q):
                Pij_true[i,i,alpha,beta] = Pi_true[i,alpha] * scra[alpha,beta]

    return [Pij_true, Pi_true, Meff]

# This version should work fine, if one changes something, please comment these
# options
@cython.boundscheck(False)
@cython.cdivision(True)
cdef list with_pc(np.ndarray[np.float32_t, ndim=4] Pij_true,
                  np.ndarray[np.float32_t, ndim=2] Pi_true,
                  float pseudocount, int N, int q, float Meff):
    cdef int i, alpha, beta
    cdef float q_float = <float> q

    cdef np.ndarray[np.float32_t, ndim=4] Pij = (1.-pseudocount) * Pij_true + pseudocount/q_float/q_float * np.ones([N,N,q,q],dtype=np.float32)
    cdef np.ndarray[np.float32_t, ndim=2] Pi = (1.-pseudocount) * Pi_true + pseudocount/q_float * np.ones([N,q],dtype=np.float32)
    cdef np.ndarray[np.float32_t, ndim=2] scra = np.eye(q,dtype=np.float32)

    for i in range(N):
        for alpha in range(q):
            for beta in range(q):
                Pij[i,i,alpha,beta] = (1.-pseudocount) * Pij_true[i,i,alpha,beta] + pseudocount/q_float * scra[alpha,beta]
    
    return [Pij, Pi]


# This version should work fine, if one changes something, please comment these
# options
@cython.boundscheck(False)
cdef int mapkey(int i, int alpha, int q):
    return ((q-1)*(i-1) + alpha)

# This version should work fine, if one changes something, please comment these
# options
@cython.boundscheck(False)
cdef list Compute_C(np.ndarray[np.float32_t, ndim=4] Pij_np,
        np.ndarray[np.float32_t, ndim=2] Pi_np, int N, int q):
    cdef int qm = q-1
    cdef np.ndarray[np.float32_t, ndim=2] C_np = np.zeros([N*qm,N*qm],dtype=np.float32)
    cdef np.ndarray[np.float32_t, ndim=2] C_self_np = np.zeros([N*qm,N*qm],dtype=np.float32)

    # memory views for efficient indexing
    cdef float[:,:] C = C_np
    cdef float[:,:] C_self = C_self_np
    cdef float [:,:,:,:] Pij = Pij_np
    cdef float [:,:] Pi = Pi_np

    cdef int i,j,alpha,beta,i_mapkey,j_mapkey

    for i in range(N):
        for j in range(N):
            for alpha in range(qm):
                for beta in range(qm):
                    i_mapkey = mapkey(i,alpha,q)
                    j_mapkey = mapkey(j,beta,q)
                    C[i_mapkey,j_mapkey] = (Pij[i,j,alpha,beta]
                            - Pi[i,alpha] * Pi[j, beta])
                    if i==j:
                        C_self[i_mapkey,j_mapkey] = C[i_mapkey,j_mapkey]

#    print('###debug Pij')
#    for item in C_self_np[0,:]:
#        print(item)
#    print("###debug")
#    np.save("cself.npy", C_self_np)
    return [C_np,C_self_np]




cdef np.ndarray ReturnW(np.ndarray[np.float32_t, ndim=2] C_np, int i,
        int j, int q):
    cdef np.ndarray[np.float32_t, ndim=2] W_np = np.zeros([q,q],dtype=np.float32)
    cdef float[:,:] W = W_np
    cdef float[:,:] C = C_np
    cdef int qm = q-1
    cdef int alpha,beta

    for alpha in range(qm):
        for beta in range(qm):
            W[alpha, beta] = C[mapkey(i,alpha,q),mapkey(i,beta,q)]
    return W_np

