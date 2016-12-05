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
form <http://autodock.scripps.edu/downloads/autodock-registration>`__).

Open a shell and type::

    cd $HOME/Downloads/
    mkdir autodock
    wget -O - http://autodock.scripps.edu/downloads/autodock-registration/tars/dist426/autodocksuite-4.2.6-x86_64Linux2.tar | tar xf - -C autodock --strip-components=1
    mv autodock/ $HOME/bin/AutoDock

From the sources
----------------

For the impatients
^^^^^^^^^^^^^^^^^^

You may download the source code from the official website, but you are
encouraged to register first your name and organization (`registration
form <http://autodock.scripps.edu/downloads/autodock-registration>`__).

Open a shell and type::

    mkdir ~/Downloads/autodock; cd $_
    wget -O - http://autodock.scripps.edu/downloads/autodock-registration/tars/dist426/autodocksuite-4.2.6-src.tar.gz | tar xfz - -C autodock --strip-components=1
    for dir in autodock autogrid
    do
        mkdir $dir/build; cd $_
        ../configure --prefix=$HOME/bin/autodock/ --bindir=$HOME/bin/autodock/
        make
        # make check
        make install
        make clean
        cd ../..
    done
    cd ~/Downloads/
    mv autodock/USERGUIDES/AutoDock4.2_UserGuide.pdf ~/
    rm -rf autodock

The command ``make check`` should return in both cases ``OK Nothing to be done
for 'check-am'``.

Installation of AutoDockTools
=============================

From the binaries
-----------------

Open a shell and type::

    cd $HOME/Downloads/
    mkdir adt
    wget -O - http://mgltools.scripps.edu/downloads/downloads/tars/releases/REL1.5.6/mgltools_x86_64Linux2_1.5.6.tar.gz | tar xfz - -C adt --strip-components=1
    cd adt
    ./install.sh -d $HOME/bin/MGLTools
    cd ..; rm -rf adt

Read the instructions dispayed at the end of the installation to add
AutoDockTools to your $PATH.

From the sources
----------------

The SWIG package (`homepage <http://www.swig.org/>`__) is a dependency of
AutoDockTools. You may download the source code from the official website,
but you are encouraged to fill in a survey first (`document
<http://swig.org/survey.html>`__).

Open a shell and type::

    sudo apt-get install -y libpcre++-dev
    cd $HOME/Downloads
    mkdir swig
    wget -O - http://downloads.sourceforge.net/project/swig/swig/swig-2.0.11/swig-2.0.11.tar.gz | tar xfz - -C swig --strip-components=1
    cd swig
    ./configure # --prefix=$HOME/bin/swig --exec_prefix=$HOME/bin/swig-exec
    make
    sudo make install
    make clean; make distclean; cd ..; rm -rf swig

Alternatively:

    sudo apt-get install -y swig

And for AutoDockTools::

    if [[ -z "$LD_LIBRARY_PATH" ]]; then
        export LD_LIBRARY_PATH="$HOME/bin/tcl"
      else
        export LD_LIBRARY_PATH="$HOME/bin/tcl:$LD_LIBRARY_PATH"
    fi
    cd $HOME/Downloads
    mkdir adt
    wget -O - http://mgltools.scripps.edu/downloads/downloads/tars/releases/REL1.5.6/mgltools_source_1.5.6.tar.gz | tar xfz - -C adt --strip-components=1
    cd adt
    sed -i '1s/python2.5/python2.7/' InstTools/install
    ./startInstallation --instDir=$HOME/bin/MGLTools <<EOF
    n
    n
    y
    EOF
    sudo sh -c 'echo "" >> /usr/lib/python2.7/sitecustomize.py; cat sitecustomize.py >> /usr/lib/python2.7/sitecustomize.py'
    cd ..; rm -rf adt
    echo -e "\n# added by $(whoami) on $(date) to source MGLTools binaries (AutoDockTools)" >> $HOME/.bashrc
    echo 'export PATH="$HOME/bin/MGLTools/bin:$PATH"' >> $HOME/.bashrc


Troubleshooting
^^^^^^^^^^^^^^^

Libraries which failed to compile are displayed at the end of the installation.
Here are some solutions to common issues:

    # mslib-1.5.6
    # x86_64-linux-gnu-gcc: error: Togl2.1/togl.c: No such file or directory
    sed -i "298,299s/2.1/2.0/g" $HOME/Downloads/adt/MGLPACKS/opengltk-1.5.6/setup.py
    # opengltk-1.5.6
    # Togl2.0/togl.c:19:20: fatal error: tclInt.h: No such file or directory
    sed -i "124,124s/#//g"      $HOME/Downloads/adt/MGLPACKS/opengltk-1.5.6/setup.py

For issues upon running ADT::

    # ImportError: cannot import name ImageTk
    # global name 'ImageTk' is not defined
    sudo apt-get install python-imaging-tk
    # NameError: global name 'ViewerFrameworkGUI' is not defined
    # ImportError: No module named idlelib
    sudo apt-get install idle-python2.7
    # RuntimeError: opengltk.OpenGL Package not found
    sed -i "72,72s|try:|try:\n  sys.path.append(\"$HOME/bin/MGLTools/MGLToolsPckgs/\")|g" $HOME/bin/MGLTools/MGLToolsPckgs/DejaVu/__init__.py
    # SyntaxError: cannot assign to __debug__
    sed 's/\(self\.\)\?__debug__/#\1__debug__/'
    sed -i '288,288s/\(self\.\)\?__debug__ = 0/\#\1__debug__ = 0/g' $HOME/bin/MGLTools/MGLToolsPckgs/ViewerFramework/VF.py

Compilation of opengltk may fail due to the way libraries :file:`Tcl_InitStubs`
and :file:`Tk_InitStubs` of tcl/tk 8.6.1 are loaded. In version 8.4, they were
named :file:`tclstub8.4` resp. :file:`tkstub8.4`. If the error message says::

    # opengltk-1.5.6
    /usr/bin/ld: cannot find -ltclstub8.4
    /usr/bin/ld: cannot find -ltkstub8.4

You'll have to manually compile file :file:`togl.so` by copy-pasting the
failing compilation line and changing the include property ``-ltclstub8.4
-ltkstub8.4`` to ``-L/usr/local/share/man/man3/Tcl_InitStubs.3
-L/usr/local/share/man/man3/Tk_InitStubs.3``::

    cd $HOME/Downloads/adt/MGLPACKS/opengltk-1.5.6
    x86_64-linux-gnu-gcc -pthread -shared -Wl,-O1 -Wl,-Bsymbolic-functions -Wl,-Bsymbolic-functions \
    -Wl,-z,relro -fno-strict-aliasing -DNDEBUG -g -fwrapv -O2 -Wall -Wstrict-prototypes -D_FORTIFY_SOURCE=2 \
    -g -fstack-protector --param=ssp-buffer-size=4 -Wformat -Werror=format-security \
    build/temp.linux-x86_64-2.7/Togl2.0/togl.o \
    build/temp.linux-x86_64-2.7/Togl2.0/toglStubInit.o \
    build/temp.linux-x86_64-2.7/Togl2.0/toglProcAddr.o \
    -L/usr/lib/tcl8.4 -L/usr/local/lib -L/usr/X11R6/lib -L/usr/X11/lib -L/usr/X11R6/lib64 -L/usr/X11/lib64 \
    -lGLU -lGL -lX11 -lXmu -lXext -lXt -lm -ldl -L/usr/local/share/man/man3/Tcl_InitStubs.3 \
    -L/usr/local/share/man/man3/Tk_InitStubs.3 -o build/lib.linux-x86_64-2.7/opengltk/OpenGL/Tk/Togl/togl.so
    cd ../..
    ./startInstallation --instDir=$HOME/bin/MGLTools <<EOF
    n
    n
    y
    EOF

If you have any other error and look forward to manually review the compilation
commands, you'll have to set the verbosity level of distutil to 1 by adding
this code to your installation script (where "opengltk" should be replaced with
a list of the failing packages)::

    sed -i '239,239s/)$/)\n    if pkgname in ["opengltk"]:\n        cmd = "DISTUTILS_DEBUG=1 " + cmd/g' $HOME/Downloads/adt/InstTools/install
 

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
<http://www.linuxfromscratch.org/blfs/view/stable/x/makedepend.html>`__),
which itself depends on pkg-config (`compilation 
<http://www.linuxfromscratch.org/blfs/view/stable/general/pkgconfig.html>`__)
and Xorg Protocol Headers (`compilation 
<http://www.linuxfromscratch.org/blfs/view/stable/x/x7proto.html>`__),
which itself depends on util-macros (`compilation 
<http://www.linuxfromscratch.org/blfs/view/stable/x/util-macros.html>`__),
which depends on Xorg (`compilation 
<http://www.linuxfromscratch.org/blfs/view/svn/x/xorg7.html#xorg-env>`__).

Once all dependencies are installed, install the Boost libraries v1.55
(`homepage <http://www.boost.org/>`__,
`doc <http://www.boost.org/doc/libs/1_55_0/?view=categorized>`__,
`download <http://www.boost.org/users/history/version_1_55_0.html>`__) with::

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


