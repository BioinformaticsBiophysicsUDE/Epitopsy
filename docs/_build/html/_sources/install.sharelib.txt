*****************
Sharing libraries
*****************

.. highlight:: bash

Numerous third-party softwares require the installation of shared libraries
(\*.so) on your system. Depending on your configuration, these libraries may
not be accessible from your shell and from Python. If you find yourself in
this situation, you will end up with error messages such as "<software_name>:
error while loading shared libraries: <library_name.so>: cannot open shared
object file: No such file or directory" everytime you try to run your software
(i.e. APBS, Clustal Omega, etc.). Several options are available to you.

First, if you have root privileges, you should modify the $PATH environment
variable of your system::

    echo -e '\n# libc default configuration\n/usr/local/lib' | sudo tee -a /etc/ld.so.conf.d/libc.conf
    sudo ldconfig # updates the $PATH

If you do not have root privileges, you might update your .bashrc file instead::

    echo -e "\n# added by $(whoami) on $(date) to grant full access to shared libraries in the shell" >> $HOME/.bashrc
    echo 'export LD_LIBRARY_PATH="/usr/local/lib"' >> $HOME/.bashrc

.. highlight:: python

This will give you (and only you) access to all shared libraries in the Linux
shell, but the Python shell will still be unable to access them, since it does
not load the .bashrc file, even with the ``shell=True`` parameter. You will
have to alter the environment variables of the suprocesses started in your
Python session. There are two ways of doing that; first you may want to set an
environment variable on a case-by-case basis, using the ``env`` parameter of
the ``subprocess.call()`` and ``subprocess.Popen()`` methods::

    >>> import os
    >>> import subprocess
    >>> myenv=os.environ.copy()
    >>> myenv['LD_LIBRARY_PATH'] = "/usr/local/lib"
    >>> p = subprocess.call(["apbs", "--version"], shell=True)
    apbs: error while loading shared libraries: libmaloc.so: cannot open shared object file: No such file or directory
    >>> p = subprocess.call(["apbs", "--version"], shell=True, env=myenv)
    Version 1.4.1

Second, you might want to set this environment variable for all subprocesses
called in a Python interpreter::

    >>> import os
    >>> import subprocess
    >>> os.environ['LD_LIBRARY_PATH'] = "/usr/local/lib"
    >>> p = subprocess.call(["apbs", "--version"], shell=True, env=myenv)
    Version 1.4.1

The latter method is more appropriate in this situation, since Epitopsy does
not allow you to pass the environment argument when calling its built-in
subprocess functions.


