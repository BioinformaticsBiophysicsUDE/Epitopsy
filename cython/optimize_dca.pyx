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

from libc.math cimport exp, log

# We now need to fix a datatype for our arrays. I've used the variable
# DTYPE for this, which is assigned to the usual NumPy runtime
# type info object.
DTYPE_int = np.int
DTYPE_float = np.float
DTYPE_cmplx = np.complex

# "ctypedef" assigns a corresponding compile-time type to DTYPE_t. For
# every type in the numpy module there's a corresponding compile-time
# type with a _t-suffix.
ctypedef np.int_t DTYPE_int_t
ctypedef np.float_t DTYPE_float_t
ctypedef np.complex_t DTYPE_cmplx_t


## special imports
import time
import sys
import gzip
import multiprocessing

from Bio import AlignIO

def dca(Z, W, outputfile, pseudocount_weight = 0.5, 
        num_threads=None, q=21.):
    '''
    Direct Coupling Analysis (DCA)
    
    Args:
        Z -> numpy array containing the alignment translated to [0,...,q]
        W -> 
        outputfile -> file containing N(N-1)/2 rows and 4 columns:
            residue i (column 1)
            residue j (column 2)
            MI(i,j) (Mutual Information between i and j)
            DI(i,j) (Direct Information between i and j)
        pseudocount_weight -> relative weight of pseudo count
            i.e. lambda = pseudocount_weight * Meff
        theta -> threshold for sequence id in reweighting

    F Morcos, A Pagnani, B Lunt, A Bertolino, DS Marks, C Sander, 
    R Zecchina, JN Onuchic, T Hwa, M Weigt (2011), Direct-coupling
    analysis of residue co-evolution captures native contacts across 
    many protein families, Proc. Natl. Acad. Sci. 108:E1293-1301.
    '''
    start_time = time.time()
    N = Z.shape[1]
    M = Z.shape[0]
    if num_threads is None:
        num_threads = multiprocessing.cpu_count()
    elif num_threads > multiprocessing.cpu_count():
        num_threads = multiprocessing.cpu_count()
    [Pij_true, Pi_true, Meff] = Compute_True_Frequencies(Z,W,M,N,q,num_threads)
    header_string = "#### optimize DCA N={0} M={1} Meff={2} q={3} pseudo_weight={4}".format(N,M,Meff,q, pseudocount_weight)
    [Pij, Pi] = with_pc(Pij_true,Pi_true,pseudocount_weight,N,q)
    C = Compute_C(Pij,Pi,N,q)
    invC = np.linalg.inv(C)
    Compute_Results(Pij, Pi,Pij_true, Pi_true, invC, N, q, outputfile, header_string)
    end_time = time.time()


#def get_entropies(inputfile, pseudocount_weight = 0.5, theta = 0,
#        num_threads=None):
#    '''
#
#    '''
#    cdef int M,N
#    dca_version = 'new' # -> dca identifier = 1
#    dca_identifier = 1
#    if num_threads is None:
#        num_threads = multiprocessing.cpu_count()
#    elif num_threads > multiprocessing.cpu_count():
#        num_threads = multiprocessing.cpu_count()
#    start_time = time.time()
#    print('read alignment ...')
#    [N, M, q, Z] = return_alignment(inputfile)
#    print('compute entropies ...')
#    entropy = Compute_Entropies(Z, M, N, q, theta, dca_identifier, num_threads)
#    return entropy



# This version should work fine, if one changes something, please comment these
# options
@cython.boundscheck(False)
cdef np.ndarray Compute_Entropies(np.ndarray[DTYPE_int_t, ndim=2] Z,
        np.ndarray[DTYPE_float_t, ndim=1] W, int M,
        int N, int q, int num_threads):
    '''
    Compute the real unmodified frequencies for the MI calculation.
    '''
    cdef int i,j,l,k,alpha,beta
    cdef long[:,:] Z_view = Z
    
    cdef float Meff = np.sum(W)
    cdef np.ndarray[DTYPE_float_t, ndim=2] Pi_true_np = np.zeros([N,q])
    cdef np.ndarray[DTYPE_float_t, ndim=1] entropy_np = np.zeros(N)
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
    cdef float s = 0
    cdef float freq = 0
    cdef float q_float = <float> q
    for i in range(N):
        s = 0
        for alpha in range(q):
            freq = Pi_true[i, alpha]
            if freq != 0:
                s = s - freq * log(freq)/log(q_float)
        entropy[i] = s
    
    return entropy_np


# This version should work fine, if one changes something, please comment these
# options
@cython.boundscheck(False)
cdef list Compute_True_Frequencies(np.ndarray[DTYPE_int_t, ndim=2] Z, 
        np.ndarray[DTYPE_float_t, ndim=1] W, int M,
        int N, int q, int num_threads):
    cdef int i,j,l,k,alpha,beta
    cdef long[:,:] Z_view = Z
    
    cdef float Meff = np.sum(W)
    cdef np.ndarray[DTYPE_float_t, ndim=4] Pij_true_np = np.zeros([N,N,q,q])
    cdef double [:,:,:,:] Pij_true = Pij_true_np
    cdef np.ndarray[DTYPE_float_t, ndim=2] Pi_true_np = np.zeros([N,q])
    cdef double [:,:] Pi_true = Pi_true_np
    cdef int z_ind
    
    for j in range(M):
        for i in range(N):
            z_ind = (Z[j,i] - 1)
            Pi_true[i, z_ind] = Pi_true[i, z_ind] + W[j]


    cdef int z_ind_i, z_ind_j 

    for l in range(M):
        for i in range(N-1):
            z_ind_i = Z[l,i] - 1 # python index
            for j in range(i+1,N):
                z_ind_j = Z[l,j] - 1 # python index
                Pij_true[i,j,z_ind_i, z_ind_j] += W[l]
                Pij_true[j,i,z_ind_j,z_ind_i] = Pij_true[i,j,z_ind_i,z_ind_j]

    for i in range(N):
        for alpha in range(q):
            Pi_true[i,alpha] = Pi_true[i,alpha] / Meff
            for j in range(N):
                for beta in range(q):
                    Pij_true[i,j,alpha,beta] = Pij_true[i,j,alpha,beta] / Meff

    cdef np.ndarray[DTYPE_float_t, ndim=2] scra_np = np.eye(q)
    cdef double[:,:] scra = scra_np
    for i in range(N):
        for alpha in range(q):
            for beta in range(q):
                Pij_true[i,i,alpha,beta] = Pi_true[i,alpha] * scra[alpha,beta]

    return [Pij_true_np, Pi_true_np, Meff]

# This version should work fine, if one changes something, please comment these
# options
@cython.boundscheck(False)
cdef list with_pc(np.ndarray[DTYPE_float_t, ndim=4] Pij_true_np, np.ndarray[DTYPE_float_t, ndim=2] Pi_true_np, float pseudocount_weight, int N, int q):
    cdef float q_float = <float> q
    cdef np.ndarray[DTYPE_float_t, ndim=4] Pij_np = (1. - pseudocount_weight)*Pij_true_np + pseudocount_weight/q_float/q_float * np.ones([N,N,q,q])
    cdef np.ndarray[DTYPE_float_t, ndim=2] Pi_np = (1. - pseudocount_weight)*Pi_true_np + pseudocount_weight/q_float * np.ones([N,q])

    cdef np.ndarray[DTYPE_float_t, ndim=2] scra_np = np.eye(q)
    cdef int i, alpha, beta
    
    # memory views for efficient indexing
    cdef double [:,:,:,:] Pij_true = Pij_true_np
    cdef double [:,:,:,:] Pij = Pij_np
    cdef double [:,:] Pi_true = Pi_true_np
    cdef double [:,:] Pi = Pi_np
    cdef double [:,:] scra = scra_np
    for i in range(N):
        for alpha in range(q):
            for beta in range(q):
                Pij[i,i,alpha,beta] = ((1.-pseudocount_weight)*Pij_true[i,i,alpha,beta]
                        + pseudocount_weight/q_float*scra[alpha,beta])
    
    return [Pij_np, Pi_np]

# This version should work fine, if one changes something, please comment these
# options
@cython.boundscheck(False)
cdef int mapkey(int i, int alpha, int q):
    return ((q-1)*(i-1) + alpha)

# This version should work fine, if one changes something, please comment these
# options
@cython.boundscheck(False)
cdef np.ndarray Compute_C(np.ndarray[DTYPE_float_t, ndim=4] Pij_np,
        np.ndarray[DTYPE_float_t, ndim=2] Pi_np, int N, int q):
    cdef int qm = q-1 
    cdef np.ndarray[DTYPE_float_t, ndim=2] C_np = np.zeros([N*(q-1),N*qm])
    
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




cdef Compute_Results(np.ndarray[DTYPE_float_t, ndim=4] Pij,
        np.ndarray[DTYPE_float_t, ndim=2] Pi,
        np.ndarray[DTYPE_float_t, ndim=4] Pij_true,
        np.ndarray[DTYPE_float_t, ndim=2] Pi_true,
        np.ndarray[DTYPE_float_t, ndim=2] invC, int N, int q,
        char* outputfile,
        char* header_string):
    output_list = []

    cdef int i,j

    for i in range(N-1):
        for j in range(i+1, N):
            ## mutual information
            MI_true, si_true, sj_true = calculate_mi(i,j,Pij_true, Pi_true, q)

            ## direct information
            W_mf = ReturnW(invC, i, j, q)
            DI_mf_pc = bp_link(i,j,W_mf,Pi,q)
            output_list.append("{0} {1} {2} {3}\n".format(i+1,j+1, MI_true, DI_mf_pc)) # python arrays ...
    
    if outputfile.endswith('.gz'):
        f = gzip.open(outputfile, 'w')
        f.write(header_string + "\n")
        f.write("i j mi di\n")
        for item in output_list:
            f.write(item)
        f.close()
    else:
        f = open(outputfile, 'w')
        f.write(header_string + "\n")
        f.write("i j mi di\n")
        for item in output_list:
            f.write(item)
        f.close()


cdef np.ndarray ReturnW(np.ndarray[DTYPE_float_t, ndim=2] invC_np, int i,
        int j, int q):
    cdef np.ndarray[DTYPE_float_t, ndim=2] W_np = np.ones([q,q])
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

cdef float bp_link(int i,int j, np.ndarray[DTYPE_float_t, ndim=2] W_mf,
        np.ndarray[DTYPE_float_t, ndim=2] Pi, int q):
    [mu1, mu2] = compute_mu(i,j,W_mf,Pi,q)
    DI = compute_di(i,j,W_mf, mu1,mu2,Pi)
    return DI

cdef list compute_mu(int i, int j, np.ndarray[DTYPE_float_t, ndim=2] W,
        np.ndarray[DTYPE_float_t, ndim=2] Pi, int q):
    cdef float epsilon = 1e-4
    cdef float diff = 1.0
    cdef float q_float = <float> q
    cdef np.ndarray[DTYPE_float_t, ndim=1] mu1 = np.ones(q)/q_float
    cdef np.ndarray[DTYPE_float_t, ndim=1] mu2 = np.ones(q)/q_float
    cdef np.ndarray[DTYPE_float_t, ndim=1] pi = Pi[i,:]
    cdef np.ndarray[DTYPE_float_t, ndim=1] pj = Pi[j,:]
    cdef np.ndarray[DTYPE_float_t, ndim=1] new1
    cdef np.ndarray[DTYPE_float_t, ndim=1] new2
    cdef np.ndarray[DTYPE_float_t, ndim=1] scra1
    cdef np.ndarray[DTYPE_float_t, ndim=1] scra2
    
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

cdef float compute_di(int i, int j, np.ndarray[DTYPE_float_t, ndim=2] W,
        np.ndarray[DTYPE_float_t, ndim=1] mu1,
        np.ndarray[DTYPE_float_t, ndim=1] mu2,
        np.ndarray[DTYPE_float_t, ndim=2] Pi):
    tiny = 1.0e-100
    Pdir = W * np.dot(mu1.reshape(-1,1),mu2.reshape(1,-1))
    Pdir = Pdir / np.sum(Pdir)
    Pfac = np.dot(Pi[i,:].reshape(-1,1), Pi[j,:].reshape(1,-1))
    DI = np.trace( np.dot(Pdir.transpose(), np.log( (Pdir+tiny) / (Pfac+tiny))))
    return DI

cdef list calculate_mi(int i, int j, np.ndarray[DTYPE_float_t, ndim=4] P2_np,
        np.ndarray[DTYPE_float_t, ndim=2] P1_np, int q):
    cdef float M = 0.
    cdef float q_float = <float> q
    cdef int alpha, beta

    # memory views
    cdef double [:,:,:,:] P2 = P2_np
    cdef double [:,:] P1 = P1_np

    for alpha in range(q):
        for beta in range(q):
            if P2[i,j,alpha, beta] > 0:
                M = M + P2[i,j,alpha,beta] * log(P2[i,j,alpha,beta] / P1[i,alpha] / P1[j,beta]) / log(q_float)

    cdef float s1 = 0.
    cdef float s2 = 0.

    for alpha in range(q):
        if P1[i,alpha] > 0:
            s1 = s1 - P1[i,alpha] * log(P1[i,alpha])/log(q_float)

        if P1[j,alpha] > 0:
            s2 = s2 - P1[j,alpha] * log(P1[j,alpha])/log(q_float)

    return [M, s1, s2]
