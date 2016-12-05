************************
Amber --- MD simulations
************************

Epitopsy cannot call Amber yet.


.. note::

    In order to run compile Amber12, following dependencies should be present on the system:

    * python2.X, gcc, gfortran, NetCDF, flex, zlib, libbz2, fftw-3.3, XBLAS, MTK++
    
    Most of them can be obtained *via* ``sudo apt-get install -y csh flex gfortran g++ xorg-dev zlib1g-dev libbz2-dev``

Note on AmberTools: AmberTools13 was developed for Amber12 and AmberTools1.5 for Amber11, but
AmberTools14 can be now be installed on any Amber starting from version 11 to 14.
Some official links:
install `AmberTools14 on Amber14 <http://jswails.wikidot.com/installing-amber14-and-ambertools14>`_,
install `AmberTools13 on Amber12 <http://jswails.wikidot.com/installing-amber12-and-ambertools-13>`_,
install `AmberTools14 together with Amber12+AmberTools13 or Amber11+AmberTools1.5
<http://jswails.wikidot.com/installing-ambertools-14-and-older-amber>`_.


Installation procedure
======================

.. highlight:: bash

Get your licensed version of Amber12 as a .tar.bz2 file and download the
AmberTools13 package from the official website.
Store them in :file:`$HOME/Downloads`. Open a shell and type::

    sudo apt-get install -y csh flex gfortran g++ xorg-dev zlib1g-dev libbz2-dev
    export AMBERHOME=$HOME/bin/amber12
    cd $HOME/Downloads
    tar xfj amber12.tar.bz2
    tar xfj AmberTools13.tar.bz2
    mv amber12 $AMBERHOME
    cd $AMBERHOME
    ./configure gnu # type "y" when prompted for downloading updates, also check that all dependencies were found
    make install
    # make test # optional
    echo -e "\n# added by $(whoami) on $(date) to source Amber12 and AmberTools13 binaries" >> $HOME/.bashrc
    echo 'export AMBERHOME=$HOME/bin/amber12' >> $HOME/.bashrc
    echo 'export PATH="$HOME/bin/amber12/bin/:$PATH"' >> $HOME/.bashrc
    echo 'export PATH="$HOME/bin/amber12/AmberTools/bin/:$PATH"' >> $HOME/.bashrc

The ``make test`` instruction may take some time depending on your machine.
The results are stored in :file:`$AMBERHOME/logs/`

The shell command ``$AMBERHOME/update_amber --version`` should output this:

.. code-block:: none

    Version is reported as <version>.<patches applied>
        AmberTools version 13.26
        Amber version 12.21

Amber Flowchart
===============

Below is a flowchart for a typical Amber simulation with a single protein.
Starting with a protein sequence file :file:`prot.fasta`, a .lib file is
generated in Xleap as well as a .rst file (atomic coordinates) and a .top
file (topology). The structure is minimized by sander, heated and used in the
production MD.

.. figure:: ../_static/figures/AMBER.*
   :target: ../_static/figures/AMBER.pdf
   :width: 556 px
   :height: 1448 px
   :scale: 80 %
   :alt: Flowchart listing all input and output files for Xleap and sander
   :align: center
   
   Legend:
   
   .pdb
       PDB coordinates
   .lib
       LIB coordinates
   .inpcrd
       input coordinates
   .prmtop
       topology file
   .in
       Amber instructions
   .frcmod
       forcefield modifications
   .mdcrd
       trajectory
   .rst
       RST coordinates
   .out
       log file


.. highlight:: python

