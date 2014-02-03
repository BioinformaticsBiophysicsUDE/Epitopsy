.. _install-index:

*******************
Installing Epitopsy
*******************

Epitopsy requires following softwares and modules to run:

    * `Biopython <http://biopython.org/wiki/Biopython>`_ (>= 1.6.0)
    * `Python <http://www.python.org/>`_ (>= 2.7.4)
    * `Cython <http://cython.org/>`_ (>= 0.18.0)
    * `Numpy <http://www.numpy.org/>`_ (>= 1.7.0)

If you have a local copy of `Anaconda <https://store.continuum.io/cshop/anaconda/>`_, then everything is already present.

Depending on your operating system, you may still need additional softwares to properly install Epitopsy.

Installation procedure
======================

.. highlight:: bash

For the impatients
------------------

Open a shell and type::

    cd Downloads/
    mkdir temp; cd temp/
    git clone https://code.google.com/p/epitopsy/
    cd epitopsy/
    python setup.py install
    conda package -u
    cd ..; rm -rf temp/

Step-by-step installation
-------------------------

You need `Git <http://git-scm.com/>`_ to be installed, if not, type in your
shell::

     sudo apt-get install git

Then move to a temporary directory and get a local copy of the installation
package from Google Code::

    cd Documents/
    mkdir temp
    cd temp/
    git clone https://code.google.com/p/epitopsy/

Move to the install directory and compile the Cython code::

    cd epitopsy/
    python setup.py install

Add the compiled package to your Anaconda development kit::

    conda package -u

Remove the installation files::

    cd ..
    rm -rf temp/

You are now ready to use Epitopsy in your Python interpreter:

.. sourcecode:: python

    import epitopsy

