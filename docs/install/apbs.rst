***********************
APBS --- Electrostatics
***********************

Epitopsy can call APBS.

Installation procedure
======================

.. highlight:: bash

For the impatients
------------------

Open a shell and type::

    wget http://downloads.sourceforge.net/project/apbs/apbs/apbs-1.4.0/APBS-1.4-source.tar.gz
    tar xfz APBS-1.4-source.tar.gz 
    cd apbs/build/
    cmake ..
    make
    sudo make install
    make clean; ./cleanup.sh; cd ../..; rm -rf apbs/ APBS-1.4-source.tar.gz

Step-by-step installation
-------------------------

You may download the source code from SourceForge, but you are encouraged to
go first on the official website of the APBS project to register your name
and organization, and tell if your work is supported by a US federal funding
(`registration form <http://www.poissonboltzmann.org/apbs/downloads>`_).

Then proceed to the download and installation::

    cd Downloads/
    wget http://downloads.sourceforge.net/project/apbs/apbs/apbs-1.4.0/APBS-1.4-source.tar.gz
    tar xvfz APBS-1.4-source.tar.gz
    cd apbs/build/
    cmake ..
    make
    sudo make install

Check that the shell command
``apbs --version 2>&1 | grep -e "Version [0-9.]*"``
outputs
``APBS 1.4.1``
, and remove the installation files and the archive::

    make clean
    ./cleanup.sh
    cd ../..
    rm -rf apbs/ APBS-1.4-source.tar.gz

.. highlight:: python


