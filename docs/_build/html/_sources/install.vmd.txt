***************************
VMD --- 3D molecular viewer
***************************

Epitopsy cannot call VMD yet.

Installation procedure
======================

.. highlight:: bash

For the impatients
------------------

You should obtain a copy of the VMD source code on the official website of the
Theoretical and Computational Biophysics Group after (`registration form
<http://www.ks.uiuc.edu/Development/Download/download.cgi>`_). Store it in
:file:`$HOME/Downloads`. Open a shell and type::

    cd $HOME/Downloads/
    tar xfz pdb2pqr-1.8.tar.gz
    cd pdb2pqr-1.8/
    ./configure --prefix $HOME/pdb2pqr/ --disable-propbka --disable-pdb2pka \
    --with-apbs=/usr/bin/apbs NUMPY=$HOME/.anaconda/lib/python2.7/site-packages/numpy/core/include/numpy/
    make
    make install
    make clean
    cd ..
    rm -rf pdb2pqr-1.8/ pdb2pqr-1.8.tar.gz
    sudo ln -s $HOME/pdb2pqr/pdb2pqr.py /usr/local/bin/pdb2pqr

Step-by-step installation
-------------------------

You should obtain a copy of the VMD source code on the official website of the
Theoretical and Computational Biophysics Group after (`registration form
<http://www.ks.uiuc.edu/Development/Download/download.cgi>`_). Store it in
:file:`$HOME/Downloads` and proceed to the compilation::

    cd $HOME/Downloads/
    tar xvfz pdb2pqr-1.8.tar.gz
    cd pdb2pqr-1.8/
    ./configure --prefix $HOME/pdb2pqr/ --disable-propbka --disable-pdb2pka \
    --with-apbs=/usr/bin/apbs NUMPY=$HOME/.anaconda/lib/python2.7/site-packages/numpy/core/include/numpy/
    make
    make test # outputs "simple test passed for PQR files"
    make install

Check that the shell command ``pdb2pqr --version`` outputs ``pdb2pqr
(Version 1.8)``, and remove the installation files and the archive::

    make clean
    cd ..
    rm -rf pdb2pqr-1.8/ pdb2pqr-1.8.tar.gz

You still have to update your environment variables to be able to call pdb2pqr
in the shell. Alternatively, you may add a symbolic link to the main.py file
in the /usr/local/bin directory::

    sudo ln -s $HOME/pdb2pqr/pdb2pqr.py /usr/local/bin/pdb2pqr

.. highlight:: python


