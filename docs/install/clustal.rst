*********************************************
Clustal Omega --- Multiple Sequence Alignment
*********************************************

Biopython provides several modules to perform a Multiple Sequence Alignment
(MSA) using the command line versions of various MSA softwares. The module
:mod:`Bio.Align.Applications.ClustalwCommandline` is an example of module
based on the Clustal family of MSA softwares. It currently supports various
releases of Clustal, including ClustalX, ClustalW, ClustalW2 and Clustal
Omega. For non-challenging MSAs, running a local copy of the Clustal software
can be more efficient than using the online website.

Clustal can be invoked in the Python shell at any time using the subprocess
module::

    >>> import subprocess
    >>> seq = "/path/seq.fasta"
    >>> msa = "/path/msa.fasta"
    >>> subprocess.call('clustalo -i {0} --infmt=fa -o {1} --force'.format(seq, msa),shell=True,env=env)

Biopython also provides a module for this task, so you don't have to
remember the exact syntax of Clustal::

    >>> from Bio.Align.Applications import ClustalOmegaCommandline
    >>> import subprocess
    >>> seq = "/path/seq.fasta"
    >>> msa = "/path/msa.fasta"
    >>> clustalo_cline = ClustalOmegaCommandline(infile=seq, outfile=msa, verbose=True, auto=True, force=True)
    >>> stdout, stderr = clustalo_cline()
    >>> print stdout
    Using 4 threads
    Read 2 sequences (type: Protein) from /path/seq.fasta
    not more sequences (2) than cluster-size (100), turn off mBed
    Setting options automatically based on input sequence characteristics (might overwrite some of your options).
    Auto settings: Enabling mBed.
    Auto settings: Setting iteration to 1.
    Progressive alignment progress: 100 % (1 out of 1)
    Progressive alignment progress done. CPU time: 0.00u 0.00s 00:00:00.00 Elapsed: 00:00:00
    Iteration step 1 out of 1
    Computing new guide tree (iteration step 1148585808)
    Computing HMM from alignment
    Progressive alignment progress: 100 % (1 out of 1)
    Progressive alignment progress done. CPU time: 0.02u 0.00s 00:00:00.02 Elapsed: 00:00:00
    Alignment written to /path/msa.fasta
    0

If :file:`/usr/local/lib` isn't available in your 'LD_LIBRARY_PATH', add the
following code in the beginning of your Python script::

    >>> import os
    >>> env = os.environ.copy()
    >>> env['LD_LIBRARY_PATH'] = "/usr/local/lib"

.. highlight:: bash


Installation procedure
======================

For the impatients
------------------

Open a shell and type::

    cd $HOME/Downloads/
    mkdir argtable; cd $_
    wget -O - sourceforge.net/projects/argtable/files/argtable/argtable-2.13/argtable2-13.tar.gz | tar xfz - -C . --strip-components=1
    ./configure
    make
    sudo make install
    make clean; cd ..; rm -rf argtable/ 
    echo -e '\n# added by $(whoami) for libargtable2.so.0' >> ~/.bashrc
    echo 'export LD_LIBRARY_PATH="/usr/local/lib:$LD_LIBRARY_PATH"' >> ~/.bashrc
    mkdir clustal; cd $_
    wget -O - http://www.clustal.org/omega/clustal-omega-1.2.1.tar.gz | tar xfz - -C . --strip-components=1
    ./configure --prefix=$HOME/bin/clustal
    make
    make install
    make clean; cd ..; rm -rf clustal/
    ln -s $HOME/bin/clustal-1.2.1/bin/clustalo $HOME/bin/clustalo

Step-by-step installation
-------------------------

First, we need to install Argtable::

    cd $HOME/Downloads/
    wget sourceforge.net/projects/argtable/files/argtable/argtable-2.13/argtable2-13.tar.gz
    tar xfz argtable2-13.tar.gz 
    cd argtable2-13/
    ./configure
    make
    make check
    sudo make install

The ``make check`` command should print this message:

.. code-block:: none

    ./test_file.sh TESTS PASSED
    ---------------------------
    PASS: test_file.sh
    ==================
    All 5 tests passed
    ==================

Check with ``ls /usr/local/lib/libargtable*`` that following files were correctly created:

.. code-block:: none

    libargtable2.a   libargtable2.so    libargtable2.so.0.1.8  
    libargtable2.la  libargtable2.so.0

You can now remove the install files and the archive::

    make clean
    cd ..
    rm -rf argtable2-13/ argtable2-13.tar.gz

The libraries should now be sourced in the local .bashrc::

    echo -e '\n# added by $(whoami) for libargtable2.so.0' >> ~/.bashrc
    echo 'export LD_LIBRARY_PATH="/usr/local/lib:$LD_LIBRARY_PATH"' >> ~/.bashrc

Next, we need to install Clustal Omega::

    cd Downloads/
    wget http://www.clustal.org/omega/clustal-omega-1.2.1.tar.gz
    tar xfz clustal-omega-1.2.1.tar.gz 
    cd clustal-omega-1.2.1/
    ./configure --prefix=$HOME/bin/clustal
    make
    make check
    make install
    make installcheck
    ln -s $HOME/bin/clustal-1.2.1/bin/clustalo $HOME/bin/clustalo

You may get warnings of the type "warning: ISO C++ does not support
variable-length array types [-Wvla]", but these can be ignored. Now check that
the shell command ``clustalo --version`` outputs ``1.2.1``, and remove the
installation files and the archive::

    make clean
    cd ..
    rm -rf clustal-omega-1.2.1/ clustal-omega-1.2.1.tar.gz

