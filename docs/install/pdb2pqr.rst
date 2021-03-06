***********************************
PDB2PQR --- Generation of PQR files
***********************************

Epitopsy provides several functions to create PQR files from PDB files using
the command line version of PDB2PQR.

.. note::

    In order to run install PDB2PQR, following
    dependency should be present on the system:

    * APBS (see :doc:`../../install/apbs`)

Installation procedure
======================

.. highlight:: bash

Before downloading the source code from SourceForge, you are encouraged to
go first on the official website of the PDB2PQR project to register your name
and organization, and tell if your work is supported by a US federal funding
(`registration form <http://www.poissonboltzmann.org/pdb2pqr/d/downloads>`_).

Version 1.9::

    cd $HOME/Downloads/
    mkdir pdb2pqr
    wget -O - http://downloads.sourceforge.net/project/pdb2pqr/pdb2pqr/pdb2pqr-1.9.0/pdb2pqr-src-1.9.0.tar.gz | tar xfz - -C pdb2pqr --strip-components=1
    cd pdb2pqr/
    python scons/scons.py PREFIX="$HOME/bin/pdb2pqr-1.9" APBS="$HOME/bin/apbs" BUILD_PDB2PKA=True
    python scons/scons.py install
    cd ..; rm -rf pdb2pqr/
    ln -s $HOME/bin/pdb2pqr-1.9/pdb2pqr.py $HOME/bin/pdb2pqr
    echo 'export APBS_PDB2PQR_DIR="$HOME/bin/pdb2pqr-1.9/"' >> ~/.bashrc # for PyMOL

Version 1.8::

    cd $HOME/Downloads/
    mkdir pdb2pqr
    wget -O - wget http://sourceforge.net/projects/pdb2pqr/files/pdb2pqr/pdb2pqr-1.8/pdb2pqr-1.8.tar.gz | tar xfz - -C pdb2pqr --strip-components=1
    cd pdb2pqr/
    ./configure --prefix $HOME/bin/pdb2pqr-1.8/ --with-apbs=$HOME/bin/apbs NUMPY=/usr/local/lib/python2.7/dist-packages/numpy # change this according to your system
    make
    #make test # outputs "simple test passed for PQR files"
    make install
    make clean; cd ..; rm -rf pdb2pqr/
    ln -s $HOME/bin/pdb2pqr-1.8/pdb2pqr.py $HOME/bin/pdb2pqr

The shell command ``pdb2pqr --version`` should give ``pdb2pqr (Version 1.8)``
or ``pdb2pqr (Version 1.9)``.

.. highlight:: python


