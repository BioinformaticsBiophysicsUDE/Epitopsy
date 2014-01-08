***************************
Multiple Sequence Alignment
***************************

Biopython provides several modules to perform a Multiple Sequence Alignment
(MSA) using the command line versions of various MSA softwares. The module
:mod:`Bio.Align.Applications.ClustalwCommandline` is an example of module
based on the Clustal family of MSA softwares. It currently supports various
releases of Clustal, including ClustalX, ClustalW, ClustalW2 and Clustal
Omega. For non-challenging MSAs, running a local copy of the Clustal software
can be more efficient than using the online website.

All Clustal releases can be run from the shell, but when Biopython is used
extensively, it may prove useful to integrate such MSAs in the python shell
itself. Clustal can be invoked at any time using the subprocess module::

    >>> import os, subprocess
    >>> env = os.environ.copy()
    >>> env['LD_LIBRARY_PATH'] = "/usr/local/lib"
    >>> subprocess.call('clustalo -i ~/Downloads/test.fasta --infmt=fa -o ~/Downloads/test2.fasta --force',shell=True,env=env)

But this is clumsy, and Biopython already provides a module for this task::

    >>> from Bio.Align.Applications import ClustalwCommandline
    >>> file = "/home/grad/Downloads/test.fasta"
    >>> clustalo_cline = ClustalwCommandline("clustalo", infile=file)
    >>> print(clustalo_cline)
    clustalo -infile=/home/grad/Downloads/test.fasta
    >>> stdout, stderr = clustalo_cline()
    TODO

Installation of Clustal Omega
=============================

.. highlight:: bash

For the impatients
------------------

Open a shell and do the following::

    wget sourceforge.net/projects/argtable/files/argtable/argtable-2.13/argtable2-13.tar.gz
    tar xfz argtable2-13.tar.gz 
    cd argtable2-13/
    ./configure
    make
    sudo make install
    make clean; cd ..; rm -rf argtable2-13/ argtable2-13.tar.gz
    wget http://www.clustal.org/omega/clustal-omega-1.2.0.tar.gz
    tar -xfz clustal-omega-1.2.0.tar.gz 
    cd clustal-omega-1.2.0/
    ./configure
    make
    sudo make install
    make clean; cd ..; rm -rf clustal-omega-1.2.0/ clustal-omega-1.2.0.tar.gz

Step-by-step installation
-------------------------

First, we need to install Argtable::

    cd Downloads/
    wget sourceforge.net/projects/argtable/files/argtable/argtable-2.13/argtable2-13.tar.gz
    tar xvfz argtable2-13.tar.gz 
    cd argtable2-13/
    ./configure
    make
    make check
    sudo make install
    make clean

If you typed ``make check`` You should have got the following output::

    ./test_file.sh TESTS PASSED
    ----------------------------------
    PASS: test_file.sh
    ==================
    All 5 tests passed
    ==================

Check that following files were correctly created ``/usr/local/lib``::

    ls /usr/local/lib/
    libargtable2.a   libargtable2.so    libargtable2.so.0.1.8  
    libargtable2.la  libargtable2.so.0  [...]

You can now remove the install files and the archive::

    cd ..
    rm -rf argtable2-13/ argtable2-13.tar.gz

.. At this point, you will have to export LD_LIBRARY_PATH everytime you use
.. clustalo, provided that your libc.conf file links to the correct dynamic
.. library::
.. 
..     cat /etc/ld.so.conf.d/libc.conf
..     # libc default configuration
..     /usr/local/lib
..     LD_LIBRARY_PATH=/usr/local/lib
..     export LD_LIBRARY_PATH

At this point, you would have to export LD_LIBRARY_PATH everytime you use
Clustal Omega to gain access to ``libargtable2.so.0``::

    LD_LIBRARY_PATH=/usr/local/lib
    export LD_LIBRARY_PATH

Alternatively, you may add this export command into your local bashrc::

    echo -e '\n# added by <username> for libargtable2.so.0 >> ~/.bashrc
    echo 'export LD_LIBRARY_PATH="/usr/local/lib"' >> ~/.bashrc

Next, we need to install Clustal Omega::

    cd Downloads/
    wget http://www.clustal.org/omega/clustal-omega-1.2.0.tar.gz
    tar xvfz clustal-omega-1.2.0.tar.gz 
    cd clustal-omega-1.2.0/
    ./configure
    make
    make check
    sudo make install
    make installcheck
    make clean

You may get warnings of the type "warning: ISO C++ does not support
variable-length array types [-Wvla]", but these can be ignored. Now
check that everything ran correctly::

    which clustalo
    /usr/local/bin/clustalo
    clustalo --help
    Clustal Omega - 1.2.0 (AndreaGiacomo)
    [...]

You can remove the installation files and the archive::

    make clean
    cd ..
    rm -rf clustal-omega-1.2.0/ clustal-omega-1.2.0.tar.gz

.. highlight:: python


