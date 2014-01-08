**************************
Gromacs --- MD simulations
**************************

Epitopsy can call Gromacs.

Installation procedure
======================

.. highlight:: bash

For the impatients
------------------

.. note::

    Before upgrading to another version of Gromacs, always remove the old
    version with ``sudo rm -rf /usr/local/gromacs`` first.

Open a shell and type::

    cd $HOME/Downloads/
    wget ftp://ftp.gromacs.org/pub/gromacs/gromacs-4.6.5.tar.gz
    tar xfz gromacs-4.6.5.tar.gz
    cd gromacs-4.6.5
    mkdir build; cd build
    cmake .. -DGMX_DOUBLE=ON -DGMX_BUILD_OWN_FFTW=ON
    make -j $(nproc)
    sudo make install
    echo -e '\nsource /usr/local/gromacs/bin/GMXRC' >> $HOME/.bashrc
    make clean; cd ../..; rm -rf gromacs-4.6.5/ gromacs-4.6.5.tar.gz

Step-by-step installation
-------------------------

.. note::

    Before upgrading to another version of Gromacs, always remove the old
    version with ``sudo rm -rf /usr/local/gromacs`` first.

You may follow the `official Gromacs installation procedure
<http://www.gromacs.org/Documentation/Installation_Instructions>`_ to get
started. It should answer all your questions. Here is the
procedure, put in a nutshell::

    cd $HOME/Downloads/
    wget ftp://ftp.gromacs.org/pub/gromacs/gromacs-4.6.5.tar.gz
    tar xvfz gromacs-4.6.5.tar.gz
    cd gromacs-4.6.5
    mkdir build
    cd build
    cmake .. -DGMX_DOUBLE=ON -DGMX_BUILD_OWN_FFTW=ON # to disable quotes: -DGMX_COOL_QUOTES=OFF
    make -j $(nproc) # or "make -j N" with N the number of CPU cores to use
    echo "Finished"
    sudo make install

You should now source the location of Gromacs in your .bashrc using the GMXRC
bash file::

    echo -e "\n# added by $(whoami) on $(date) to source gromacs executables" >> $HOME/.bashrc
    echo -e '\nsource /usr/local/gromacs/bin/GMXRC' >> $HOME/.bashrc

..    echo -e '\nsource /usr/local/gromacs/bin/GMXRC' | sudo tee -a /etc/bash.bashrc # does not change anything, the python interpreter does not read it!!

Check if the shell command
``mdrun_d -version 2>/dev/null | grep -e "VERSION [0-9\.]*"``
outputs
``VERSION 4.6.5``
, and remove the installation files and the archive::

    make clean
    cd ../..
    rm -rf gromacs-4.6.5/ gromacs-4.6.5.tar.gz

To completely erase Gromacs from your machine, you may follow the `official
Gromacs uninstallation procedure
<http://www.gromacs.org/Documentation/Removing_Installation>`_.

.. highlight:: python


