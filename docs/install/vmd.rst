***************************
VMD --- 3D molecular viewer
***************************

Epitopsy cannot call VMD yet.

.. highlight:: bash

Installation procedure from the binaries
========================================

For the impatients
------------------

You should obtain a copy of the VMD pre-compiled sources on the official
website of the Theoretical and Computational Biophysics Group (`registration
form <http://www.ks.uiuc.edu/Development/Download/download.cgi>`_). Store it
in :file:`$HOME/Downloads`. Open a shell and type::

    cd $HOME/Downloads
    tar xfz vmd-1.9.1.bin.LINUXAMD64.opengl.tar.gz
    cd vmd-1.9.1/
    ./configure LINUXAMD64 # change here according to your architecture
    cd src/
    sudo make install
    cd ../..
    rm -rf vmd-1.9.1/ vmd-1.9.1.bin.LINUXAMD64.opengl.tar.gz

Step-by-step installation
-------------------------

You should obtain a copy of the VMD pre-compiled sources on the official
website of the Theoretical and Computational Biophysics Group (`registration
form <http://www.ks.uiuc.edu/Development/Download/download.cgi>`_). Store it
in :file:`$HOME/Downloads`.

To extract the archive::

    cd $HOME/Downloads
    tar xfz vmd-1.9.1.bin.LINUXAMD64.opengl.tar.gz
    cd vmd-1.9.1/

When necessary (e.g. if you don't have root privileges), change the
installation location::

    export VMDINSTALLNAME="vmd"
    export VMDINSTALLBINDIR="$HOME/bin"
    export VMDINSTALLLIBRARYDIR="$HOME/vmd"

Compile and install::

    ./configure LINUXAMD64 # change here according to your architecture
    cd src/
    sudo make install
    cd ../..
    rm -rf vmd-1.9.1/ vmd-1.9.1.bin.LINUXAMD64.opengl.tar.gz

Installation procedure from the sources
=======================================

.. warning::

    JN: Could not make it work.

For the impatients
------------------

.. seealso::

    The `TCBG official compilation documentation <http://www.ks.uiuc.edu/Research/vmd/doxygen/>`_.

You first need `Tcl <http://www.tcl.tk/>`_ installed::

    wget http://prdownloads.sourceforge.net/tcl/tcl8.6.1-src.tar.gz
    tar xfz tcl8.6.1-src.tar.gz
    cd tcl8.6.1/unix/
    ./configure --enable-threads --enable-shared --enable-64bit \
        --prefix=$HOME/test/tlc --exec_prefix=$HOME/test/tlc-exec
    make
    make install
    make clean

You should obtain a copy of the VMD source code on the official website of the
Theoretical and Computational Biophysics Group after (`registration form
<http://www.ks.uiuc.edu/Development/Download/download.cgi>`_). Store it in
:file:`$HOME/Downloads`. Open a shell and type::

    cd $HOME/Downloads/
    tar xfz vmd-1.9.1.src.tar.gz
    cd plugins/
    make LINUXAMD64 # change here according to your architecture; print the list with "make"
    # got an error message, I'm stuck here
    export TCLLIB=/home/grad/test/tlc-exec/lib/
    export TCLINC=/home/grad/test/tlc/include/
    cd ../vmd-1.9.1/
    ./configure
    make
    [...]
    rm -rf vmd-1.9.1/ plugins/ vmd-1.9.1.src.tar.gz

.. highlight:: python


