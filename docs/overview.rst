********
Overview
********

Goals
=====

Epitopsy is a Python package designed to extend the reach of the `Biopython
<http://biopython.org/wiki/Biopython>`_ environment. Epitopsy supports
additional file formats (PQR, lattice protein representation) and additional
third party softwares (APBS, Clustal Omega, PDB2PQR, *etc.*, :doc:`see list
<install/third-party>`), and implements various algorithms relevant to
the field of bioinformatics (see for example :doc:`Direct Coupling Analysis
<packages/tools/dca>` or :doc:`FFT-driven protein-ligand energy scanning
<epitopsy.EnergyGrid>`).

Unless otherwise stated, the whole project is open source (see
:doc:`license <license>`) and available on
`Google Code <https://code.google.com/p/epitopsy/>`_.

History
=======

The Epitopsy project started in 2006 with the diploma thesis of J. Nikolaj
Dybowski at the Center of Advanced European Studies and Research (CAESAR).
The original programming language was Java, but the project was ported to
Python in 2011 thanks to the efforts of Christoph Wilms, during his time as
a *Ph.D.* student at the Centre for Medical Biotechnology (Zentrum für
Medizinische Biotechnologie, ZMB), Faculty of Biology -- University of
Duisburg-Essen.

The code was extended in 2012 by integrating features from the `Biopython
<http://biopython.org/wiki/Biopython>`_ package and implementing new
techniques found in the scientific literature. The project is currently
maintained by Ludwig Ohl and Jean-Noël Grad, *Ph.D.* students at the ZMB.

Technical aspects
=================

This package shares common functionalities with `Biopython
<http://biopython.org/wiki/Biopython>`_, yet provides much more new features.
It is further enhanced using `Cython <http://cython.org/>`_ for faster
execution times and `NumPy <http://www.numpy.org/>`_ for most operations on
arrays. `Sphinx <http://sphinx-doc.org/>`_ was used to generate this website.

Development team:

    * Christoph Wilms <christoph.wilms@uni-due.de>
    * Ludwig Ohl <ludwig.ohl@uni-due.de>
    * Jean-Noël Grad <jean-noel.grad@uni-due.de>

Further references:

    * Dybowski, J. N. *Development of a method for optimal superposition of
      pairs of similar macromolecules*, diploma thesis, Fachhochschule Bingen,
      2006.
    * Wilms, C. *Methods for the prediction of complex biomolecular
      structures*, PhD thesis, Universität Duisburg-Essen, 2013.
    * Cock, P. J. *et al.* *Biopython: freely available Python tools for
      computational molecular biology and bioinformatics*, *Bioinformatics*
      **2009**, *25*, 1422-1423.
    * Behnel, S. *et al.* *Cython: The Best of Both Worlds*, *Comput Sci Eng.*
      **2011**, *13*, 31--39.
    * Oliphant, T. E. *Python for Scientific Computing*, *Comput Sci Eng.*
      **2007**, *9*, 10--20.


