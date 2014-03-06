*****************************
PyMOL --- 3D molecular viewer
*****************************

You may open PyMOL in a minimal GUI using this code::

    import pymol
    pymol.finish_launching()
    # further code waiting for user to quit the GUI

But this minimalist GUI lacks the command prompt, and the script execution is
frozen until the user exits the window. Alternatively, the PyMOL core can be
started in the command line mode::

    import __main__
    __main__.pymol_argv = ['pymol','-qc'] # -c: command line mode, -q: quiet mode (low verbosity)
    import pymol
    
    if __name__ == "__main__":
        # start PyMOL if not already running in the current interpreter
        try:
            a = pymol.glutThread # see if PyMOL's variables are loaded in memory
            del a
        except:
            pymol.finish_launching()
        
        # PyMOL instructions...

.. note::

    In order to use the latest version of PyMOL in Anaconda, see how to update
    your Anaconda library in the :ref:`updating-anaconda` section below.

Installation procedure
======================

.. highlight:: bash

For the impatients
------------------

.. note::

    See the :ref:`troubleshooting` section below if you have Anaconda
    installed on your system.

Open a shell and type::

    sudo apt-get install build-essential python-dev python-pmw \
                         libglew-dev freeglut3-dev libpng-dev libfreetype6-dev -y
    cd $HOME/Downloads/
    wget http://downloads.sourceforge.net/project/pymol/pymol/1.6/pymol-v1.6.0.0.tar.bz2
    tar xfj pymol-v1.6.0.0.tar.bz2; cd pymol/
    prefix=/usr/share/pymol
    sudo python setup.py build install --{home,install-scripts}=$prefix  --install-lib=$prefix/modules
    sudo ln -s $prefix/pymol /usr/bin/pymol
    cd ..; rm -rf pymol/ pymol-v1.6.0.0.tar.bz2

Step-by-step installation
-------------------------

.. note::

    See the :ref:`troubleshooting` section below if you have Anaconda
    installed on your system.

You may follow the `PyMOL Wiki installation procedure
<http://www.pymolwiki.org/index.php/Linux_Install>`_ to get the compilation
done from the svn repository. The procedure described below uses the
SourceForge repository.

First retrieve the dependencies and the archive::

    sudo apt-get install build-essential python-dev python-pmw \
                         libglew-dev freeglut3-dev libpng-dev libfreetype6-dev -y
    # sudo apt-get autoremove
    cd $HOME/Downloads/
    wget http://downloads.sourceforge.net/project/pymol/pymol/1.6/pymol-v1.6.0.0.tar.bz2
    tar xvfj pymol-v1.6.0.0.tar.bz2
    cd pymol/

Build PyMOL and add a soft link in your ``/usr/bin/``::

    prefix=/usr/share/pymol
    sudo python setup.py build install --home=$prefix --install-lib=$prefix/modules --install-scripts=$prefix
    sudo ln -s $prefix/pymol /usr/bin/pymol

Alternatively, if you prefer a local copy or if you lack root privileges::

    prefix=$HOME/bin/pymol
    python setup.py build install --home=$prefix --install-lib=$prefix/modules --install-scripts=$prefix
    echo -e "\n# added by $(whoami) on $(date) to source PyMOL binaries" >> $HOME/.bashrc
    echo 'export PATH="$prefix/:$PATH"' >> $HOME/.bashrc

Check if the shell command ``pymol -c | grep -e "Version [0-9\.]*"`` outputs
``Version 1.6.0.0``, and remove the installation files and the archive::

    cd ..
    rm -rf pymol/ pymol-v1.6.0.0.tar.bz2

Optionally, remove the default splash screen (or replace it by one of your own
flavor)::

    mv $prefix/modules/pymol/pymol_path/data/pymol/splash{,_old}.png

.. _updating-anaconda:

Updating Anaconda
=================

Anaconda is sometimes shipped with an old version of the PyMOL module. Check
if your Anaconda library is up to date by typing this command in the shell:
``python -c "import __main__;__main__.pymol_argv = ['pymol','-c'];import
pymol;pymol.finish_launching()" | grep -e "Version [0-9\.]*"`` --- if the
output is ``PyMOL(TM) Molecular Graphics System, Version 1.5.0.1.``, feel free
to create a symbolic link to the newer module::

    oldversion=$(python -c "import pymol;print pymol.__file__")  # /home/user/.anaconda/lib/python2.7/pymol/__init__.pyc
    oldversion=${oldversion%/*}                                  # /home/user/.anaconda/lib/python2.7/pymol
    mv $oldversion ${oldversion}_v1.5                            # backup
    ln -s /usr/share/pymol/modules/pymol/ $oldversion            # for a hard copy: cp -r

.. _troubleshooting:

Troubleshooting
===============

If there is a copy of Anaconda installed on your system, you may encounter an
issue while compiling PyMOL. Typically, the compilation will run smoothly, but
the PyMOL build obtained cannot be correctly executed, as it will search for a
pymol module inside :file:`.anaconda/lib/python2.7/`.
When running this build of PyMOL from the shell, several errors will print out,
usually involving shared libraries which cannot be accessed (libpng15.so.15,
libc.so.6, or _cmd.so). The pymol module is actually located inside the
``$prefix`` directory, since it is where all PyMOL files were copied to.

In order to solve this issue, you must use another installation of python,
which you should already have if you did ``sudo apt-get python-dev`` as
explained. Then, temporarily deactivate shell sourcing for Anaconda by
manually editing the following line in your :file:`.bashrc` file::

    export PATH="$HOME/.anaconda/bin:$PATH"

Add a `#` character in front of it, save the file and start a new shell
interpreter to build PyMOL. The correct pymol module directory will be used
by your build of PyMOL. Once finished, simply remove the `#` character in
the :file:`.bashrc` file.

You might also find yourself unable to start PyMOL directly from your $HOME
folder (error message: ``*** buffer overflow detected ***:
/usr/bin/python terminated``), but any other directory within your $HOME will
do it. The only workaround is not to start PyMOL from your $HOME. You
can achieve this by opening the launcher (gedit $prefix/pymol) and typing
the following if statement anywhere before the line executing PyMOL::

    if [ $PWD = $HOME ]; then
      cd Documents/
      echo "Moved to $HOME/Documents/ since this build of PyMOL cannot start directly from your home folder."
    fi

Anaconda v1.70 may be unable to import the PyMOL module (see this `bug report
<https://groups.google.com/a/continuum.io/forum/#!topic/anaconda/-DLG2ZdTkw0>`_).
In this case, the following error message is displayed:

.. code-block:: none

    Traceback (most recent call last):
    File "<stdin>", line 1, in <module>
    File "/home/user/.anaconda/lib/python2.7/pymol/__init__.py", line 491, in <module>
    from pymol import _cmd
    ImportError: /home/user/.anaconda/bin/../lib/libm.so.6: version `GLIBC_2.15' not found (required by /home/user/.anaconda/lib/python2.7/pymol/_cmd.so)

A workaround provided by the Anaconda development team consists in deleting
(or renaming) the :file:`libm.so.6` symbolic link from the Anaconda directory::

    mv $HOME/.anaconda/lib/libm.so.6 $HOME/.anaconda/lib/libm.so.6-old

Anaconda might still be unable to import PyMOL if the chempy package is
missing, with the following error message:

.. code-block:: none

    Error: unable to initalize the pymol.cmd module
    Traceback (most recent call last):
    File "/home/user/.anaconda/lib/python2.7/pymol/cmd.py", line 117, in <module>
    from chempy import io
    ImportError: No module named chempy

Correct this by copying the chempy package from your Python installation::

    cp -r /usr/share/pymol/modules/chempy $HOME/.anaconda/lib/python2.7/site-packages/chempy


.. highlight:: python


