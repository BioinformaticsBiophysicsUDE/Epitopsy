
..
    rsync -auv sewczyk@132.252.170.160:Amber12.tar.bz2 .

************************
Amber --- MD simulations
************************

Epitopsy cannot call Amber yet.

Installation procedure
======================

.. highlight:: bash

For the impatients
------------------

Get your licensed version of Amber12 as a .tar.bz2 file and download the
AmberTools13 package from the official website
(`registration form <http://ambermd.org/AmberTools-get.html>`_).
Store them in $HOME/Downloads for example. Open a shell and type::

    sudo apt-get install csh flex gfortran g++ xorg-dev \
                         zlib1g-dev libbz2-dev
    cd $HOME
    tar xfj $HOME/Downloads/AmberTools13.tar.bz2
    tar xfj $HOME/Downloads/Amber12.tar.bz2
    export AMBERHOME=$HOME/amber12
    cd $AMBERHOME
    ./configure gnu
    make install
    make test


Step-by-step installation
-------------------------

Get your licensed version of Amber12 as a .tar.bz2 file and download the
AmberTools13 package from the official website
(`registration form <http://ambermd.org/AmberTools-get.html>`_).
Store them in $HOME/Downloads for example.

Download all dependencies and extract the two Amber archives in your
installation directory::

    sudo apt-get install csh flex gfortran g++ xorg-dev \
                         zlib1g-dev libbz2-dev
    INSTALLDIR=$HOME # installation directory, may be changed
    cd $INSTALLDIR
    tar xvfj $HOME/Downloads/AmberTools13.tar.bz2
    tar xvfj $HOME/Downloads/Amber12.tar.bz2

Compile Amber using these variables (you may change $INSTALLDIR)::

    export AMBERHOME=$INSTALLDIR/amber12
    cd $AMBERHOME
    ./configure gnu
    make install
    make test # takes some time...

The ``make test`` instruction may take some time depending on your machine.
The results are stored in $AMBERHOME/logs/

The shell command ``$AMBERHOME/update_amber --version`` should output this::

    Version is reported as <version>.<patches applied>
        AmberTools version 13.22
        Amber version 12.21



