********************
tcl/tk --- GUI tools
********************

The tcl/tk package is necessary for VMD and for AutoDockTools.

More elaborated instructions are available on `Beyond Linux From Scratch
<http://www.linuxfromscratch.org/blfs/>`_ for 
`tcl <http://www.linuxfromscratch.org/blfs/view/stable/general/tcl.html>`_
and `tk <http://www.linuxfromscratch.org/blfs/view/stable/general/tk.html>`_.

.. highlight:: bash

Installation procedure from the binaries
========================================

For the impatients
------------------
Open a shell and type::

    cd $HOME/Downloads
    wget http://prdownloads.sourceforge.net/tcl/tcl8.6.1-src.tar.gz
    tar xfz tcl8.6.1-src.tar.gz
    cd tcl8.6.1/unix/
    ./configure --enable-threads --enable-shared \
        $([ $(uname -m) = x86_64 ] && echo --enable-64bit) \
        --prefix=$HOME/bin/tlc --exec_prefix=$HOME/bin/tlc-exec
    make; make install
    make clean; cd ../..; rm -rf tcl8.6.1/ tcl8.6.1-src.tar.gz 
    wget http://downloads.sourceforge.net/project/tcl/Tcl/8.6.1/tk8.6.1-src.tar.gz
    tar xfz tk8.6.1-src.tar.gz
    cd tk8.6.1/unix/
    ./configure --enable-threads --enable-shared \
        $([ $(uname -m) = x86_64 ] && echo --enable-64bit) \
        --prefix=$HOME/bin/tk --exec_prefix=$HOME/bin/tk-exec \
        --with-tcl=$HOME/bin/tlc-exec/lib
    make; make install
    make clean; cd ../..; rm -rf tk8.6.1/ tk8.6.1-src.tar.gz

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


.. highlight:: python


