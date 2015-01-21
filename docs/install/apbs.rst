***********************
APBS --- Electrostatics
***********************

Epitopsy can call APBS.

Installation procedure
======================

.. highlight:: bash

You may download the source code from SourceForge, but you are encouraged to
go first on the official website of the APBS project to register your name
and organization, and tell if your work is supported by a US federal funding
(`registration form <http://www.poissonboltzmann.org/apbs/downloads>`_). In
the shell::

    cd $HOME/Downloads/
    mkdir apbs
    wget -O - http://downloads.sourceforge.net/project/apbs/apbs/apbs-1.4.0/APBS-1.4-source.tar.gz | tar xfz - -C apbs --strip-components=1
    cd apbs/build/
    cmake .. -DCMAKE_INSTALL_PREFIX:PATH=$HOME/bin/apbs-1.4
    make
    make install
    make clean; cd ../..; rm -rf apbs/
    ln -s $HOME/bin/apbs-1.4/bin/apbs $HOME/bin/apbs
    echo LD_LIBRARY_PATH="$HOME/bin/apbs-1.4/lib:$LD_LIBRARY_PATH" >> ~/.bashrc

The shell command ``apbs --version 2>&1 | grep -e "Version [0-9.]*"`` shoud give ``1.4.1``.

.. highlight:: python


