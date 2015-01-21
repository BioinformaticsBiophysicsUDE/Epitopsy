********************
tcl/tk --- GUI tools
********************

The tcl/tk package is necessary for VMD and for AutoDockTools.

More elaborated instructions are available on `Beyond Linux From Scratch
<http://www.linuxfromscratch.org/blfs/>`_ for 
`tcl <http://www.linuxfromscratch.org/blfs/view/stable/general/tcl.html>`_
and `tk <http://www.linuxfromscratch.org/blfs/view/stable/general/tk.html>`_.

.. highlight:: bash

Installation procedure from the sources
=======================================

Open a shell and type::

    prefix_tcl="$HOME/bin/tcl"
    prefix_tk="$HOME/bin/tk"
    cd $HOME/Downloads
    wget http://prdownloads.sourceforge.net/tcl/tcl8.6.1-src.tar.gz
    wget http://downloads.sourceforge.net/project/tcl/Tcl/8.6.1/tk8.6.1-src.tar.gz
    tar xfz tcl8.6.1-src.tar.gz
    tar xfz tk8.6.1-src.tar.gz
    cd tcl8.6.1/unix/
    ./configure --enable-threads --enable-shared \
        $([ $(uname -m) = x86_64 ] && echo --enable-64bit) \
        --prefix=$prefix_tcl --exec_prefix=${prefix_tcl}-exec
    make; make install
    make clean; cd ../..
    cd tk8.6.1/unix/
    ./configure --enable-threads --enable-shared \
        $([ $(uname -m) = x86_64 ] && echo --enable-64bit) \
        --prefix=$prefix_tk --exec_prefix=${prefix_tk}-exec \
        --with-tcl=${prefix_tcl}-exec/lib
    make; make install
    make clean; cd ../..
    rm -rf tcl8.6.1/ tcl8.6.1-src.tar.gz tk8.6.1/ tk8.6.1-src.tar.gz # remove tcl8.6.1/ only at the very end (files required for the compilation of tk)


.. highlight:: python


