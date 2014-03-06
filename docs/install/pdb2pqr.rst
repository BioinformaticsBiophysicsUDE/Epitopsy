***********************************
PDB2PQR --- Generation of PQR files
***********************************

..
    gedit /etc/environment
    PATH="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games"
    echo 'pdb2pqr="/home/grad/pdb2pqr/pdb2pqr.py"' >> /etc/environment

Epitopsy provides several functions to create PQR files from PDB files using
the command line version of PDB2PQR.

Installation procedure
======================

.. highlight:: bash

For the impatients
------------------

Open a shell and type::

    cd $HOME/Downloads/
    wget http://sourceforge.net/projects/pdb2pqr/files/pdb2pqr/pdb2pqr-1.8/pdb2pqr-1.8.tar.gz
    tar xfz pdb2pqr-1.8.tar.gz
    cd pdb2pqr-1.8/
    ./configure --prefix $HOME/bin/pdb2pqr/ --with-apbs=/usr/bin/apbs NUMPY=$HOME/.anaconda/lib/python2.7/ # or any directory containing site-packages/numpy/
    make
    make install
    make clean
    cd ..
    rm -rf pdb2pqr-1.8/ pdb2pqr-1.8.tar.gz
    ln -s $HOME/bin/pdb2pqr/pdb2pqr.py $HOME/bin/pdb2pqr/pdb2pqr
    echo -e "\n# added by $(whoami) on $(date) to source PDB2PQR executables" >> $HOME/.bashrc
    echo 'export PATH="$HOME/bin/pdb2pqr:$PATH"' >> $HOME/.bashrc

Step-by-step installation
-------------------------

You may download the source code from SourceForge, but you are encouraged to
go first on the official website of the PDB2PQR project to register your name
and organization, and tell if your work is supported by a US federal funding
(`registration form <http://www.poissonboltzmann.org/pdb2pqr/d/downloads>`_).

Then proceed to the download and installation::

    cd $HOME/Downloads/
    wget http://sourceforge.net/projects/pdb2pqr/files/pdb2pqr/pdb2pqr-1.8/pdb2pqr-1.8.tar.gz
    tar xvfz pdb2pqr-1.8.tar.gz
    cd pdb2pqr-1.8/
    ./configure --prefix $HOME/bin/pdb2pqr/ --with-apbs=/usr/bin/apbs NUMPY=$HOME/.anaconda/lib/python2.7/ # or any directory containing site-packages/numpy/
    make
    make test # outputs "simple test passed for PQR files"
    make install

Check that the shell command ``pdb2pqr --version`` outputs ``pdb2pqr
(Version 1.8)``, and remove the installation files and the archive::

    make clean
    cd ..
    rm -rf pdb2pqr-1.8/ pdb2pqr-1.8.tar.gz

You still have to update your :file:`.bashrc` to be able to call pdb2pqr in
the shell::

    ln -s $HOME/bin/pdb2pqr/pdb2pqr.py $HOME/bin/pdb2pqr/pdb2pqr
    echo -e "\n# added by $(whoami) on $(date) to source PDB2PQR executables" >> $HOME/.bashrc
    echo 'export PATH="$HOME/bin/pdb2pqr:$PATH"' >> $HOME/.bashrc

.. highlight:: python


