***************************
VMD --- 3D molecular viewer
***************************

Epitopsy cannot call VMD yet.

.. highlight:: bash

Installation procedure from the binaries
========================================

You should obtain a copy of the VMD pre-compiled sources on the official
website of the Theoretical and Computational Biophysics Group (`registration
form <http://www.ks.uiuc.edu/Development/Download/download.cgi>`_). Store it
in :file:`$HOME/Downloads`. Open a shell and type::

    export VMDINSTALLNAME="vmd"
    export VMDINSTALLBINDIR="$HOME/bin/vmd-1.9.1"
    export VMDINSTALLLIBRARYDIR="$HOME/bin/vmd-1.9.1"
    cd $HOME/Downloads
    mkdir vmd
    tar xfz vmd-1.9.1.bin.LINUXAMD64.opengl.tar.gz -C vmd --strip-components=1
    cd vmd
    ./configure LINUXAMD64 # change here according to your architecture
    cd src/
    make install
    cd ../..; rm -rf vmd
    ln -s $HOME/bin/vmd-1.9.1/vmd $HOME/bin/vmd


Installation procedure from the sources
=======================================

.. note::

    In order to compile VMD, following dependency must be present on the system:

    * tcl/tk (see :doc:`../../install/tcl_tk`)

.. seealso::

    The `TCBG official compilation documentation <http://www.ks.uiuc.edu/Research/vmd/doxygen/>`_.

You should obtain a copy of the VMD source code on the official website of the
Theoretical and Computational Biophysics Group after (`registration form
<http://www.ks.uiuc.edu/Development/Download/download.cgi>`_). Store it in
:file:`$HOME/Downloads`. Open a shell and type::

    export VMDINSTALLNAME="vmd"
    export VMDINSTALLBINDIR="$HOME/bin/vmd-1.9.1"
    export VMDINSTALLLIBRARYDIR="$HOME/bin/vmd-1.9.1"
    cd $HOME/Downloads/
    mkdir vmd
    tar xfz vmd-1.9.1.src.tar.gz -C vmd
    cd vmd/plugins/
    for file in $(ls $HOME/bin/tcl/include/*.h); do cp $file ./include; done # only useful if tcl.h and other headers from $TCLINC aren't found by make
    make LINUXAMD64 TCLINC="$HOME/bin/tcl/include/" TCLLIB="$HOME/bin/tcl-exec/lib/" NETCDFLIB="" NETCDFINC="" NETCDFLDFLAGS="" # change this line according to your architecture; print the list with "make"
    # got an error message with webpdbplugin.o, I'm stuck here
    export TCLLIB=/home/grad/test/tcl-exec/lib/
    export TCLINC=/home/grad/test/tcl/include/
    cd ../vmd-1.9.1/
    ./configure --prefix=$HOME/bin/vmd-1.9.1
    make
    make install
    cd ../..; rm -rf vmd
    ln -s $HOME/bin/vmd-1.9.1/vmd $HOME/bin/vmd

.. highlight:: python


