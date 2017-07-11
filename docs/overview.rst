********
Overview
********

Description
===========

Epitopsy is a Python package designed to carry out :doc:`FFT-driven
protein-ligand energy scanning <epitopsy.EnergyGrid>` and :doc:`Direct
Coupling Analysis <packages/tools/dca>`.

The whole project is open source (see :doc:`license <license>`) and available
on `GitHub <https://github.com/BioinformaticsBiophysicsUDE/Epitopsy>`_.

History
=======

The Epitopsy project started in 2006 with the diploma thesis of Jan Nikolaj
Dybowski at the Center of Advanced European Studies and Research (CAESAR).
The original programming language was Java, but the project was ported to
Python in 2011 by Christoph Wilms, *Ph.D.* student at the Centre for
Medical Biotechnology, University of Duisburg-Essen.

The code was extended in 2012 by integrating features from the `Biopython
<http://biopython.org/wiki/Biopython>`_ package and implementing new
techniques found in the scientific literature. The project is currently
maintained by Ludwig Ohl and Jean-Noël Grad, *Ph.D.* students at the 
Bioinformatics and Computational Biophysics department, University of
Duisburg-Essen.

Technical aspects
=================

This package shares common functionalities with `Biopython
<http://biopython.org/wiki/Biopython>`_. It is further enhanced using
`Cython <http://cython.org/>`_ for faster execution and
`NumPy <http://www.numpy.org/>`_ for most operations on arrays. It also
provides an interface to common bioinformatics tools (APBS, Clustal Omega,
PDB2PQR, *etc.*, :doc:`see list <install/third-party>`).
`Sphinx <http://sphinx-doc.org/>`_ was used for documentation.

Current development team:

    * Ludwig Ohl <ludwig.ohl@uni-due.de>
    * Jean-Noël Grad <jean-noel.grad@uni-due.de>

References on Epitopsy:

    * Jean-Noël Grad, Alba Gigante, Christoph Wilms, Jan Nikolaj Dybowski,
      Ludwig Ohl, Christian Ottmann, Carsten Schmuck and Daniel Hoffmann,
      *Locating large flexible ligands on proteins*, **2017**,
      `arXiv:1707.02614 <https://arxiv.org/abs/1707.02614>`_
    * Christoph Wilms, *Methods for the Prediction of Complex Biomolecular
      Structures* `[online]
      <https://duepublico.uni-duisburg-essen.de/servlets/DocumentServlet?id=33166>`_,
      PhD thesis, University of Duisburg-Essen, **2013**.
    * Jan Nikolaj Dybowski, *Development of a method for optimal superposition
      of pairs of similar macromolecules*, diploma thesis, Fachhochschule
      Bingen, **2006**.

Further references:

    * Cock, P. J. *et al.* *Biopython: freely available Python tools for
      computational molecular biology and bioinformatics*, *Bioinformatics*
      **2009**, *25*, 1422-1423.
    * Behnel, S. *et al.* *Cython: The Best of Both Worlds*, *Comput Sci Eng.*
      **2011**, *13*, 31--39.
    * Oliphant, T. E. *Python for Scientific Computing*, *Comput Sci Eng.*
      **2007**, *9*, 10--20.

