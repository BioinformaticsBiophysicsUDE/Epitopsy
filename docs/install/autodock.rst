***********************************
AutoDock Vina --- Molecular Docking
***********************************

Epitopsy does not support AutoDock yet.

.. highlight:: bash

Installation of AutoDock
========================

From the binaries
-----------------

You may download the binaries from the official website, but you are
encouraged to register first your name and organization (`registration
form <http://autodock.scripps.edu/downloads/autodock-registration>`_).

Open a shell and type::

    cd $HOME/Downloads/
    wget http://autodock.scripps.edu/downloads/autodock-registration/tars/dist4251/autodocksuite-4.2.5.1-x86_64Linux2.tar.gz
    tar xfz autodocksuite-4.2.5.1-x86_64Linux2.tar.gz
    mv i86_64Linux2/ $HOME/bin/AutoDock
    rm autodocksuite-4.2.5.1-x86_64Linux2.tar.gz

From the sources
----------------

For the impatients
^^^^^^^^^^^^^^^^^^

You may download the source code from the official website, but you are
encouraged to register first your name and organization (`registration
form <http://autodock.scripps.edu/downloads/autodock-registration>`_).

Open a shell and type::

    cd $HOME/Downloads/
    wget http://autodock.scripps.edu/downloads/autodock-registration/tars/dist4251/autodocksuite-4.2.5.1-src.tar.gz
    tar xfz autodocksuite-4.2.5.1-src.tar.gz 
    for file in {autodock,autogrid}
    do
        cd $HOME/Downloads/src/$file
        mkdir x86_64Linux2; cd $_
        ../configure # --prefix=$HOME/
        make
        sudo make install
        make clean
    done
    cd $HOME/Downloads/
    mv src/USERGUIDES/AutoDock4.2_UserGuide.pdf $HOME/
    rm -rf src/ autodocksuite-4.2.5.1-src.tar.gz

Step-by-step installation
^^^^^^^^^^^^^^^^^^^^^^^^^

You may download the source code from the official website, but you are
encouraged to register first your name and organization (`registration
form <http://autodock.scripps.edu/downloads/autodock-registration>`_).

Then proceed to the download and installation::

    cd $HOME/Downloads/
    wget http://autodock.scripps.edu/downloads/autodock-registration/tars/dist4251/autodocksuite-4.2.5.1-src.tar.gz
    tar xvfz autodocksuite-4.2.5.1-src.tar.gz 
    cd src/autodock
    mkdir x86_64Linux2; cd $_
    ../configure
    make
    make check
    sudo make install
    make clean
    cd ../../autogrid/
    mkdir x86_64Linux2; cd $_
    ../configure
    make
    make check
    sudo make install
    make clean
    cd ../../..
    mv src/USERGUIDES/AutoDock4.2_UserGuide.pdf $HOME/
    rm -rf src/ autodocksuite-4.2.5.1-src.tar.gz

The command ``make check`` should return in both cases ``OK Nothing to be done
for 'check-am'``. The installation directory is :file:`/usr/local/bin` by
default, but it can be changed using the commutator ``../configure
--prefix=$HOME/`` to install in :file:`$HOME/bin`.

Installation of AutoDockTools
=============================

From the binaries
-----------------

Open a shell and type::

    cd $HOME/Downloads/
    wget http://mgltools.scripps.edu/downloads/downloads/tars/releases/REL1.5.6/mgltools_x86_64Linux2_1.5.6.tar.gz
    tar xfz mgltools_x86_64Linux2_1.5.6.tar.gz
    cd mgltools_x86_64Linux2_1.5.6
    ./install.sh -d $HOME/bin/MGLTools
    cd ..; rm -rf mgltools_x86_64Linux2_1.5.6/ mgltools_x86_64Linux2_1.5.6.tar.gz

Read the instructions dispayed at the end of the installation to add
AutoDockTools to your $PATH.

From the sources
----------------

For the impatients
^^^^^^^^^^^^^^^^^^

The SWIG package (`homepage <http://www.swig.org/>`_) is a dependency of
AutoDockTools. You may download the source code from the official website,
but you are encouraged to fill in a survey first (`document
<http://swig.org/survey.html>`_).

Open a shell and type::

    cd $HOME/Downloads
    wget http://downloads.sourceforge.net/project/swig/swig/swig-2.0.11/swig-2.0.11.tar.gz
    tar xfz swig-2.0.11.tar.gz
    cd swig-2.0.11
    ./configure # --prefix=$HOME/bin/swig --exec_prefix=$HOME/bin/swig-exec
    make
    sudo make install
    make clean; make distclean; cd ..; rm -rf swig-2.0.11/ swig-2.0.11.tar.gz

And for AutoDockTools::

    if [[ -z "$LD_LIBRARY_PATH" ]]; then
        export LD_LIBRARY_PATH="$HOME/bin/tlc-exec/lib:$HOME/bin/tk-exec/lib/"
      else
        export LD_LIBRARY_PATH="$HOME/bin/tlc-exec/lib:$HOME/bin/tk-exec/lib/:$LD_LIBRARY_PATH"
    fi
    cd $HOME/Downloads
    wget http://mgltools.scripps.edu/downloads/downloads/tars/releases/REL1.5.6/mgltools_source_1.5.6.tar.gz
    tar xfz mgltools_source_1.5.6.tar.gz
    cd mgltools_source_1.5.6/
    sed -i '1s/python2.5/python2.7/' InstTools/install
    ./startInstallation --instDir=$HOME/bin/MGLTools <<EOF
    n
    n
    y
    EOF
    cd ..; rm -rf mgltools_source_1.5.6/ mgltools_source_1.5.6.tar.gz

You should now source the location of MGLTools executables in your
:file:`.bashrc` file::

    echo -e "\n# added by $(whoami) on $(date) to source MGLTools executables (AutoDockTools)" >> $HOME/.bashrc
    echo 'export PATH="$HOME/bin/MGLTools/bin:$PATH"' >> $HOME/.bashrc

Installation of Vina
====================

From the binaries
-----------------

Open a shell and type::

    cd $HOME/Downloads
    wget http://vina.scripps.edu/download/autodock_vina_1_1_2_linux_x86.tgz
    tar xfz autodock_vina_1_1_2_linux_x86.tgz
    mv autodock_vina_1_1_2_linux_x86 $HOME/bin/vina
    rm autodock_vina_1_1_2_linux_x86.tgz

You should now source the location of Vina in your :file:`.bashrc` file::

    echo -e "\n# added by $(whoami) on $(date) to source vina executables" >> $HOME/.bashrc
    echo 'export PATH="$HOME/bin/vina/bin:$PATH"' >> $HOME/.bashrc

From the sources
----------------

.. warning:: Not tested.

Vina depends on makedepend (`compilation
<http://www.linuxfromscratch.org/blfs/view/stable/x/makedepend.html>`_),
which itself depends on pkg-config (`compilation
<http://www.linuxfromscratch.org/blfs/view/stable/general/pkgconfig.html>`_)
and Xorg Protocol Headers (`compilation
<http://www.linuxfromscratch.org/blfs/view/stable/x/x7proto.html>`_),
which itself depends on util-macros (`compilation
<http://www.linuxfromscratch.org/blfs/view/stable/x/util-macros.html>`_),
which depends on Xorg (`compilation
<http://www.linuxfromscratch.org/blfs/view/svn/x/xorg7.html#xorg-env>`_).

Once all dependencies are installed, install the Boost libraries v1.55
(`homepage <http://www.boost.org/>`_,
`doc <http://www.boost.org/doc/libs/1_55_0/?view=categorized>`_,
`download <http://www.boost.org/users/history/version_1_55_0.html>`_) with::

    cd $HOME/Downloads
    tar xfz boost_1_55_0.tar.gz
    cd boost_1_55_0
    ./bootstrap.sh
    ./b2
    ./b2 install --prefix=$HOME/bin/boost/ --prefix-exec=$HOME/bin/boost-exec/
    cd ..; rm -rf boost_1_55_0/ boost_1_55_0.tar.gz

To install Vina, open a shell and type::

    cd $HOME/Downloads
    wget http://vina.scripps.edu/download/autodock_vina_1_1_2.tgz
    tar xfz autodock_vina_1_1_2.tgz
    cd autodock_vina_1_1_2/build/linux/release/
    sed -i 's/BOOST_VERSION=1_41/BOOST_VERSION=1_55/' Makefile
    sed -i 's/BOOST_INCLUDE = \$(BASE)\/include/BOOST_INCLUDE = \$HOME\/bin\/boost\/include/' Makefile
    sed -i 's/GPP=\/usr\/local\/bin\/g++/GPP=\/usr\/bin\/g++/' Makefile
    make depend
    make
    [to be continued...]

.. highlight:: python


