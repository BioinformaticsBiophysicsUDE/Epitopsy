*****************************
PyMOL --- 3D molecular viewer
*****************************

To import the PyMOL module::

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

Python module ``chempy`` is a dependency of PyMOL, but it is shipped with it.

Installation procedure
======================

.. highlight:: bash

.. note::

    See the :ref:`troubleshooting` section if you experience any difficulty.

You may follow the `PyMOL Wiki installation procedure
<http://www.pymolwiki.org/index.php/Linux_Install>`_ to get the compilation
done from the svn repository. The procedure described below uses the
SourceForge repository.

Version 1.7.X::

    sudo apt-get install -y build-essential python-dev python-pmw \
                            libglew-dev freeglut3-dev libpng-dev libfreetype6-dev
    version=1.7.4.0
    prefix=$HOME/bin/pymol-${version}
    mkdir $HOME/Downloads/pymol; cd $_
    wget -O - http://downloads.sourceforge.net/project/pymol/pymol/${version:0:3}/pymol-v${version}.tar.bz2 | tar xfj - -C . --strip-components=1
    python setup.py build install --home=$prefix --install-lib=$prefix/modules --install-scripts=$prefix
    cd ..; rm -rf pymol
    ln -s $prefix/pymol $HOME/bin/pymol
    mv $prefix/modules/pymol/pymol_path/data/pymol/splash{,_old}.png # remove splash screen (optional)
    ln -s ${prefix}/modules/pymol $HOME/.local/lib/python2.7/site-packages/pymol
    ln -s ${prefix}/modules/chempy $HOME/.local/lib/python2.7/site-packages/chempy # if not already installed

The command ``pymol -c | grep -e "Version [0-9\.]*"`` should output
``$version``.

.. _troubleshooting:

Troubleshooting
===============

Shared libraries cannot be accessed (libpng15.so.15, libc.so.6, or _cmd.so)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Description: several shared libraries cannot be imported when starting the
PyMOL GUI (libpng15.so.15, libc.so.6, _cmd.so).

Explanation: PyMOL binaries compiled using python2.7 from the Anaconda package
will try to import the PyMOL python module from
:file:`.anaconda/lib/python2.7/`, but it is actually stored in ``$prefix``.

Solution: build with python-dev. First deactivating shell sourcing for
Anaconda by manually editing the following line in your :file:`.bashrc` file::

    #export PATH="$HOME/.anaconda/bin:$PATH"

Add a `#` character to comment out the line, save and start a new shell
interpreter to build PyMOL. The correct PyMOL module directory will be used
by your build of PyMOL. Once finished, simply remove the `#` character in
the :file:`.bashrc` file.

Cannot start PyMOL GUI from $HOME (buffer overflow detected)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Description: the GUI refuses to launch in the $HOME folder (error message:
``*** buffer overflow detected ***: /usr/bin/python terminated``), but any
other directory will do it (PyMOL < v1.7.0.0).

Solution: the only workaround is not to start PyMOL from your $HOME. You
can achieve this by opening the launcher (gedit $prefix/pymol) and typing
the following if statement anywhere before the line executing PyMOL::

    if [ $PWD = $HOME ]; then
      cd Documents/
      echo "Moving to $HOME/Documents/ since this build of PyMOL cannot start directly from your home folder."
    fi

Cannot import the PyMOL module in Anaconda due to GLIBC_2.15
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Description: python2.7 from Anaconda v1.70 refuses to import the PyMOL module
and displays this error message:

.. code-block:: none

    Traceback (most recent call last):
    File "<stdin>", line 1, in <module>
    File "/home/user/.anaconda/lib/python2.7/pymol/__init__.py", line 491, in <module>
    from pymol import _cmd
    ImportError: /home/user/.anaconda/bin/../lib/libm.so.6: version `GLIBC_2.15' not found (required by /home/user/.anaconda/lib/python2.7/pymol/_cmd.so)

Solution: A workaround provided by the Anaconda development team (see this
`bug report
<https://groups.google.com/a/continuum.io/forum/#!topic/anaconda/-DLG2ZdTkw0>`_)
consists in deleting (or renaming) the :file:`libm.so.6` symbolic link from
the Anaconda directory::

    mv $HOME/.anaconda/lib/libm.so.6 $HOME/.anaconda/lib/libm.so.6-old

Cannot import the PyMOL module in Anaconda due to chempy
""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Description: python2.7 from Anaconda v1.70 refuses to import the chempy module
and displays this error message:

.. code-block:: none

    Error: unable to initalize the pymol.cmd module
    Traceback (most recent call last):
    File "/home/user/.anaconda/lib/python2.7/pymol/cmd.py", line 117, in <module>
    from chempy import io
    ImportError: No module named chempy

Solution: the chempy module is missing and needs to be copied from the normal
python modules::

    cp -r /usr/share/pymol/modules/chempy $HOME/.anaconda/lib/python2.7/site-packages/chempy

Python imports an outdated version of PyMOL
"""""""""""""""""""""""""""""""""""""""""""

Error: Anaconda is sometimes shipped with an old version of the PyMOL module.

Solution: check if the Anaconda library is up to date by typing this command::

    ``python -c "import __main__;__main__.pymol_argv = ['pymol','-c'];import
    pymol;pymol.finish_launching()" | grep -e "Version [0-9\.]*"``

If the output is ``PyMOL(TM) Molecular Graphics System, Version 1.5.0.1.``, 
create a symbolic link to the newer module::

    oldversion=$(python -c "import pymol;print pymol.__file__")  # stored: /home/user/.anaconda/lib/python2.7/pymol/__init__.pyc
    oldversion=${oldversion%/*}                                  # stored: /home/user/.anaconda/lib/python2.7/pymol
    mv $oldversion ${oldversion}_v1.5                            # make a backup
    ln -s /usr/share/pymol/modules/pymol/ $oldversion

Cannot open the APBS plugin window
""""""""""""""""""""""""""""""""""

Description: the APBS plugin crashes due to a conflict in Tkinter and displays
this error message (PyMOL < v1.7.4.0):

.. code-block:: none

    Traceback (innermost last):
      File "/usr/lib/python2.7/dist-packages/Pmw/Pmw_1_3/lib/PmwBase.py", line 1747, in __call__
        return apply(self.func, args)
      File "/home/user/bin/pymol-1.7.2.1/modules/pmg_tk/startup/apbs_tools.py", line 680, in __init__
        group.pack(fill='both',expand=1, padx=4, pady=5)

Solution: methods :meth:`Pmw.Group.pack` and :meth:`Pmw.Group.grid` shouldn't
be used together with Tk >= 8.6 (see the first warning in
`<http://effbot.org/tkinterbook/grid.htm>`_). Comment out all calls to this
method in :file:`apbs_tools.py`::

    sed -i "s/group\.grid/#group.grid/" $HOME/bin/pymol-1.7.2.1/modules/pmg_tk/startup/apbs_tools.py

.. highlight:: python



