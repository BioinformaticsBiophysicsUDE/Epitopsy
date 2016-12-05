*********************************************
OpenBabel --- Molecular file format converter
*********************************************

To get the package and its Python module::

    sudo apt-get install -y openbabel python-openbabel


Installation procedure
======================

.. highlight:: bash

From the sources
----------------

Eigen is a dependency. To locally install version 3.2.2::

    mkdir ~/Downloads/eigen; cd $_
    wget -O - http://bitbucket.org/eigen/eigen/get/3.2.2.tar.bz2 | tar xjf - -C . --strip-components=1
    mkdir build_dir; cd $_
    cmake .. -DCMAKE_INSTALL_PREFIX=$HOME/bin/eigen-3.2.2
    make install
    make clean; cd ../..; rm -rf eigen

OpenBabel::

    # dependencies
    sudo apt-get install -y libcairo2-dev libwxgtk2.8-dev libwxgtk2.8-dbg libxml2-dev
    mkdir ~/Downloads/obabel; cd $_
    wget -O - http://downloads.sourceforge.net/project/openbabel/openbabel/2.3.2/openbabel-2.3.2.tar.gz | tar xfz - -C . --strip-components=1
    mkdir build; cd $_
    cmake .. -DEIGEN3_INCLUDE_DIR=$HOME/bin/eigen-3.2.2/include/eigen3 -DPYTHON_BINDINGS=ON -DCMAKE_INSTALL_PREFIX=$HOME/bin/openbabel-2.3.2
    make -j$(nproc)
    make install
    make clean; cd ../..; rm -rf openbabel
    echo -e "\n# added by $(whoami) on $(date) to source openbabel" >> $HOME/.bashrc
    echo 'export PATH="$PATH:$HOME/bin/openbabel-2.3.2/bin"' >> ~/.bashrc

The shell command ``obabel -V`` should output ``Open Babel 2.3.2``.


.. highlight:: python

