
:mod:`dca` --- Direct Coupling Analysis operations
==================================================

.. module:: dca
   :synopsis: Direct Coupling Analysis operations.
.. moduleauthor:: Christoph Wilms <christoph.wilms@uni-due.de>
.. sectionauthor:: Jean-Noel Grad <jean-noel.grad@uni-due.de>


This module provides direct coupling analysis operations.


.. seealso::
        
    Direct-coupling analysis of residue co-evolution captures native 
    contacts across many protein families
        F Morcos, A Pagnani, B Lunt, A Bertolino, DS Marks, C Sander,
        R Zecchina, JN Onuchic, T Hwa, M Weigt (2011), Proc. Natl.
        Acad. Sci. 108:E1293-1301.


.. _dca-syntax:

Module Syntax
-------------

Empty.

.. _contents-of-module-dca:

Module Contents
---------------

.. function:: dca(inputfile, outputfile)

    Direct Coupling Analysis (DCA)

    :param file inputfile: file with fasta formatted sequences
    :param file outputfile: file containing N(N-1)/2 rows and 4 columns
        
    Structure of **outputfile**:

    * residue i (column 1)
    * residue j (column 2)
    * MI(i,j) (Mutual Information between i and j)
    * DI(i,j) (Direct Information between i and j)
    
    
.. function:: get_similarity(ref, query, theta)
    
    Does not check the length of both objects.
    
    :param ref: iterable object which will be compared with query
    :type ref: obj
    :param query: iterable object which will be compared with ref
    :type query: obj
    :param theta: reweight sequences with identity threshold
    :type theta: float
    
    :returns: 0 if greater than theta, 1 if less than theta.
    
    
.. function:: return_alignment(inputfile)
    
    Docstring missing.
    
    
.. function:: Compute_True_Frequencies(Z,M,N,q, theta)
    
    Docstring missing.
    
    
.. function:: with_pc(Pij_true,Pi_true,pseudocount_weight,N,q)
    
    Docstring missing.
    
    
.. function:: Compute_C(Pij,Pi,N,q)
    
    Docstring missing.
    
    
.. function:: Compute_Results(Pij, Pi,Pij_true, Pi_true, invC, N, q, outputfile)
    
    Docstring missing.
    
    
.. function:: calculate_mi(i,j,P2,P1,q)
    
    Docstring missing.
    
    
.. function:: ReturnW(invC, i, j, q)
    
    Docstring missing.
    
    
.. function:: bp_link(i,j,W_mf,Pi,q)
    
    Docstring missing.
    
    
.. function:: compute_mu(i,j,W,Pi,q)
    
    Docstring missing.
    
    
.. function:: compute_di(i,j,W, mu1,mu2,Pi)
    
    Docstring missing.
    
    
