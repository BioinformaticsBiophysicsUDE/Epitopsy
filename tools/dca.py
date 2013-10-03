import time
import gzip
import sys

import numpy as np
from scipy.spatial import distance

from Bio import AlignIO


def dca(inputfile, outputfile):
    '''
    Direct Coupling Analysis (DCA)

    Args:
        inputfile -> file with fasta formatted sequences
        outputfile -> file containing N(N-1)/2 rows and 4 columns:
            residue i (column 1)
            residue j (column 2)
            MI(i,j) (Mutual Information between i and j)
            DI(i,j) (Direct Information between i and j)

    F Morcos, A Pagnani, B Lunt, A Bertolino, DS Marks, C Sander,
    R Zecchina, JN Onuchic, T Hwa, M Weigt (2011), Direct-coupling
    analysis of residue co-evolution captures native contacts across
    many protein families, Proc. Natl. Acad. Sci. 108:E1293-1301.
    '''
    # relative weight of pseudo count -> lambda = pseudocount_weight * Meff
    pseudocount_weight = 0.5
    # threshold for sequence id in reweighting
    theta = 0.2
    print('read alignment ...')
    [N, M, q, Z] = return_alignment(inputfile)
    print('compute frequencies ...')
    [Pij_true, Pi_true, Meff] = Compute_True_Frequencies(Z,M,N,q, theta)
    # free memory -> Z can be large for "good" alignments i.e. many sequences
    del(Z)
    print("#### DCA N={0} M={1} Meff={2} q={3}".format(N,M,Meff,q))
    with open(outputfile, 'w') as f:
        f.write("#### DCA N={0} M={1} Meff={2} q={3}\n".format(N,M,Meff,q))
        # r header
        f.write("i j mi di\n")

    [Pij, Pi] = with_pc(Pij_true,Pi_true,pseudocount_weight,N,q)
    print('compute C ...')
    C = Compute_C(Pij,Pi,N,q)
    print('invert C ...')
    invC = np.linalg.inv(C)
    print('compute DI+MI ...')
    Compute_Results(Pij, Pi,Pij_true, Pi_true, invC, N, q, outputfile)

def return_alignment(inputfile):
    if inputfile.endswith('.fa'):
        alignment = AlignIO.read(inputfile, "fasta")
    elif inputfile.endswith('.fa.gz'):
        with gzip.open(inputfile) as f:
            alignment = AlignIO.read(f, "fasta")

    # num of sequences
    M = len(alignment)
    # num of alignment positions
    N = len(alignment[0])
    q = 21
    aa_dict = {'-':1, 'A':2, 'C':3, 'D':4, 'E':5, 'F':6, 'G':7, 'H':8, 'I':9,
            'K':10, 'L':11, 'M':12, 'N':13, 'P':14, 'Q':15, 'R':16, 'S':17,
            'T':18, 'V':19, 'W':20, 'Y':21}
    Z = np.zeros([M,N])
    for i, record in enumerate(alignment):
        sequence = record.seq.tostring()
        for j,char in enumerate(sequence):
            if char in aa_dict:
                score = aa_dict[char]
            else:
                score = 1
            Z[i,j] = score

    return [N,M,q,Z]


def get_similarity(ref, query, theta):
    '''
    Does not check the length of both objects.
    Args:
        ref -> iterable object which will be compared with query
        query -> iterable object which will be compared with ref

    Returns:
        0 if greater than theta, 1 if less than theta.
    '''
    counter = 0.
    n_items = float(len(ref))
    for r,q in zip(ref, query):
        if r == q:
            counter += 1.

    score = 1. - counter / n_items

    if score > theta:
        return 0.
    else:
        return 1.


def Compute_True_Frequencies(Z,M,N,q, theta):
    W = np.ones(M)
    if (theta > 0.0):
        W = (1./(1.+np.sum(distance.squareform((distance.pdist(Z,'hamming') < theta).astype(float)),0)))

#        ## modification, because for many sequences there will be an memory error
#        hamming_sum = np.zeros(M)
#        for i in range(M):
#            for j in range(i+1,M):
#                #hamming_dist = (distance.pdist([Z[i],Z[j]],'hamming') < theta).astype(float)
#                hamming_dist = get_similarity(Z[i], Z[j], theta)
#                hamming_sum[i] += hamming_dist
#                hamming_sum[j] += hamming_dist
#        W = (1./(1.+hamming_sum))


    Meff = np.sum(W)
    Pij_true = np.zeros([N,N,q,q])
    Pi_true = np.zeros([N,q])

    for j in range(M):
        for i in range(N):
            z_ind = Z[j,i] - 1
            Pi_true[i, z_ind] = Pi_true[i, z_ind] + W[j]

    Pi_true = Pi_true / Meff

    for l in range(M):
        for i in range(N-1):
            for j in range(i+1,N):
                z_ind_i = Z[l,i] - 1 # python index
                z_ind_j = Z[l,j] - 1 # python index
                Pij_true[i,j,z_ind_i, z_ind_j] += W[l]
                Pij_true[j,i,z_ind_j,z_ind_i] = Pij_true[i,j,z_ind_i,z_ind_j]

    Pij_true = Pij_true / Meff

    scra = np.eye(q)
    for i in range(N):
        for alpha in range(int(q)):
            for beta in range(int(q)):
                Pij_true[i,i,alpha,beta] = Pi_true[i,alpha] * scra[alpha,beta]

    return [Pij_true, Pi_true, Meff]

def with_pc(Pij_true,Pi_true,pseudocount_weight,N,q):
    Pij = (1. - pseudocount_weight)*Pij_true + pseudocount_weight/q/q * np.ones([N,N,q,q])
    Pi = (1. - pseudocount_weight)*Pi_true + pseudocount_weight/q * np.ones([N,q])

    scra = np.eye(q)

    for i in range(N):
        for alpha in range(int(q)):
            for beta in range(int(q)):
                Pij[i,i,alpha,beta] = ((1.-pseudocount_weight)*Pij_true[i,i,alpha,beta]
                        + pseudocount_weight/q*scra[alpha,beta])

    return [Pij, Pi]

def Compute_C(Pij,Pi,N,q):
    def mapkey(i, alpha, q):
        return ((q-1)*(i-1) + alpha)

    qm = int(q)-1
    C = np.zeros([N*(q-1),N*qm])
    for i in range(N):
        for j in range(N):
            for alpha in range(qm):
                for beta in range(qm):
                    C[mapkey(i,alpha,q),mapkey(j,beta,q)] = (Pij[i,j,alpha,beta]
                            - Pi[i,alpha] * Pi[j, beta])

    return C




def Compute_Results(Pij, Pi,Pij_true, Pi_true, invC, N, q, outputfile):
    output_list = []
    for i in range(N-1):
        for j in range(i+1, N):
            ## mutual information
            MI_true, si_true, sj_true = calculate_mi(i,j,Pij_true, Pi_true, q)

            ## direct information
            W_mf = ReturnW(invC, i, j, q)
            DI_mf_pc = bp_link(i,j,W_mf,Pi,q)
            output_list.append("{0} {1} {2} {3}\n".format(i+1,j+1, MI_true, DI_mf_pc)) # python arrays ...

    with open(outputfile, 'a') as f:
        for item in output_list:
            f.write(item)


def calculate_mi(i,j,P2,P1,q):
    M = 0.
    for alpha in range(q):
        for beta in range(q):
            if P2[i,j,alpha, beta] > 0:
                M = M + P2[i,j,alpha,beta] * np.log(P2[i,j,alpha,beta] / P1[i,alpha] / P1[j,beta])/np.log(q)

    s1 = 0.
    s2 = 0.
    for alpha in range(q):
        if P1[i,alpha] > 0:
            s1 = s1 - P1[i,alpha] * np.log(P1[i,alpha])/np.log(q)

        if P1[j,alpha] > 0:
            s2 = s2 - P1[j,alpha] * np.log(P1[j,alpha])/np.log(21)

    return [M, s1, s2]

def ReturnW(invC, i, j, q):
    W = np.ones([q,q])
    qm = int(q)-1
    qmi = qm * (i-1)
    qmj = qm * (j-1)
    for alpha in range(qm):
        qmia = qmi + alpha
        for beta in range(qm):
            W[alpha, beta] = np.exp(-invC[qmia, qmj+beta])
    return W

def bp_link(i,j,W_mf,Pi,q):
    [mu1, mu2] = compute_mu(i,j,W_mf,Pi,q)
    DI = compute_di(i,j,W_mf, mu1,mu2,Pi)
    return DI

def compute_mu(i,j,W,Pi,q):
    epsilon = 1e-4
    diff = 1.0
    mu1 = np.ones(q)/q
    mu2 = np.ones(q)/q
    pi = Pi[i,:]
    pj = Pi[j,:]

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

def compute_di(i,j,W, mu1,mu2,Pi):
    tiny = 1.0e-100
    Pdir = W * np.dot(mu1.reshape(-1,1),mu2.reshape(1,-1))
    Pdir = Pdir / np.sum(Pdir)
    Pfac = np.dot(Pi[i,:].reshape(-1,1), Pi[j,:].reshape(1,-1))
    DI = np.trace( np.dot(Pdir.transpose(), np.log( (Pdir+tiny) / (Pfac+tiny))))
    return DI


if __name__ == "__main__":
    start_time = time.time()
    error_string = "Usage:\npython dca.py -in=<alignment_path> -out=<di_output_path>"
    if len(sys.argv) != 3:
        print(error_string)

    inputfile = None
    outputfile = None
    for item in sys.argv[1:]:
        content = item.split('=')
        if content[0] == '-in':
            inputfile = content[1]
        elif content[0] == '-out':
            outputfile = content[1]

    if inputfile is None or outputfile is None:
        raise AttributeError(error_string)

    dca(inputfile, outputfile)
    end_time = time.time()
    print('finished after {0}s'.format(end_time - start_time))


