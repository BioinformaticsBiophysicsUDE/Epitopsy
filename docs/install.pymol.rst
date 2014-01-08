*****************************
PyMOL --- 3D molecular viewer
*****************************

Epitopsy can call PyMOL.

Installation procedure
======================

.. highlight:: bash

For the impatients
------------------

.. note::

    See the :ref:`PyMOL-troubleshooting` section below if you have Anaconda
    installed on your system.

Open a shell and type::

    sudo apt-get install build-essential python-dev python-pmw \
                         libglew-dev freeglut3-dev libpng-dev libfreetype6-dev
    cd $HOME/Downloads/
    wget http://downloads.sourceforge.net/project/pymol/pymol/1.6/pymol-v1.6.0.0.tar.bz2
    tar xfj pymol-v1.6.0.0.tar.bz2; cd pymol/
    prefix=$HOME/pymol
    python setup.py build install --{home,install-scripts}=$prefix  --install-lib=$prefix/modules
    sudo ln -s $prefix/pymol /usr/bin/pymol
    cd ..; rm -rf pymol/ pymol-v1.6.0.0.tar.bz2

Step-by-step installation
-------------------------

.. note::

    See the :ref:`PyMOL-troubleshooting` section below if you have Anaconda
    installed on your system.

You may follow the `PyMOL Wiki installation procedure
<http://www.pymolwiki.org/index.php/Linux_Install>`_ to get the compilation
done from the svn repository. The procedure described below uses the
SourceForge repository.

First retrieve the dependencies and the archive::

    sudo apt-get install build-essential python-dev python-pmw \
                         libglew-dev freeglut3-dev libpng-dev libfreetype6-dev
    # sudo apt-get autoremove
    cd $HOME/Downloads/
    wget http://downloads.sourceforge.net/project/pymol/pymol/1.6/pymol-v1.6.0.0.tar.bz2
    tar xvfj pymol-v1.6.0.0.tar.bz2
    cd pymol/

Compile PyMOL using these variables (you may set ``prefix=/usr/share/pymol``
to make PyMOL available to everyone, in which case the build should be done
with root privileges)::

    prefix=$HOME/pymol
    python setup.py build install --home=$prefix --install-lib=$prefix/modules --install-scripts=$prefix

Create a soft link to the executable in your ``/usr/bin/``::

    sudo ln -s $prefix/pymol /usr/bin/pymol

Check if the shell command ``pymol -c | grep -e "Version [0-9\.]*"`` outputs
``Version 1.6.0.0``, and remove the installation files and the archive::

    cd ..
    rm -rf pymol/ pymol-v1.6.0.0.tar.bz2

Optionally, remove the default splash screen (or replace it by one of your own
flavor)::

    mv $prefix/modules/pymol/pymol_path/data/pymol/splash{,_old}.png

.. _PyMOL-troubleshooting:

Troubleshooting
---------------

If there is a copy of Anaconda installed on your system (which we strongly
recommend), you may encounter an issue while compiling PyMOL. Typically,
the compilation will run smoothly, but the PyMOL build obtained cannot be
correctly executed, as it will search for a pymol module inside
:file:`.anaconda/lib/python2.7/`. When running this build of PyMOL from the
shell, several errors will print out, usually involving shared libraries which
cannot be accessed (libpng15.so.15, libc.so.6, or _cmd.so). The pymol module
is actually located inside the ``$prefix`` directory, since it is where all
PyMOL files were copied to.

To resolve this, you must use another installation of python, which you should
already have if you did ``sudo apt-get python-dev`` as explained. Then,
temporarily deactivate shell sourcing for Anaconda by manually editing the
following line in your :file:`.bashrc` file::

    export PATH="/home/<your_username>/.anaconda/bin:$PATH"

Add a `#` character in front of it, save the file and start a new shell
interpreter to build PyMOL. The correct pymol module directory will be used
by your build of PyMOL. Once finished, simply remove the `#` character in
the :file:`.bashrc` file.

You might also find yourself unable to start PyMOL directly from your $HOME
folder (error message: ``*** buffer overflow detected ***:
/usr/bin/python terminated``), but any other directory within your $HOME will
do it. This is independent from the location of the PyMOL directory as well
as from the sourcing of your Anaconda distribution.
The only workaround is not to start PyMOL from your $HOME. You
can achieve this by opening the launcher (``gedit $prefix/pymol``) and typing
the following ``if`` statement anywhere before the line executing PyMOL::

    if [ $PWD = $HOME ]; then
      cd Documents/
      echo "Moved to $HOME/Documents/ since this build of PyMOL cannot start directly from your home folder."
    fi

.. highlight:: python



