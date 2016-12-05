**********************************
Graph-Tool --- Graphs and Networks
**********************************

Graph-Tool 2.2.28 (`homepage <http://graph-tool.skewed.de/>`__, `doc <http://graph-tool.skewed.de/static/doc/index.html>`__, `download <http://graph-tool.skewed.de/download>`__).

Installation procedure
======================

.. highlight:: bash

For the impatients
------------------

Graph-Tool can be installed with a package. Get first the public key `98507F25
<http://pgp.skewed.de:11371/pks/lookup?op=get&search=0x612DEFB798507F25>`__ and
add it to your list::

    sudo apt-key add 0x612defb798507f25.html

The command ``apt-key list`` should output:

.. code-block:: none

    pub   4096R/98507F25 2013-10-17 [expires: 2018-10-16]
    uid                  Tiago de Paula Peixoto <tiago@skewed.de>
    uid                  Tiago de Paula Peixoto <tiago@itp.uni-bremen.de>
    sub   4096R/1A7ECE03 2013-10-17 [expires: 2018-10-16]
    sub   4096R/23F08CAF 2013-10-17 [expires: 2018-10-16]

If you skip this step, the command ``sudo apt-get update`` will later return
the following warning (this will not prevent you from installing the package):

.. code-block:: none

    W: GPG error: http://downloads.skewed.de saucy Release: The following
    signatures couldn't be verified because the public key is not available:
    NO_PUBKEY 92F371361A7ECE03

Proceed to the installation of the package::

    DISTRIBUTION=$(cat /etc/lsb-release | grep -P -o "(?<=CODENAME=)[a-z]+")
    echo "deb http://downloads.skewed.de/apt/$DISTRIBUTION $DISTRIBUTION universe" | sudo tee -a /etc/apt/sources.list
    echo "deb-src http://downloads.skewed.de/apt/$DISTRIBUTION $DISTRIBUTION universe" | sudo tee -a /etc/apt/sources.list
    sudo apt-get update
    sudo apt-get install python-graph-tool
    cp -r /usr/lib/python2.7/dist-packages/graph_tool/ $HOME/.anaconda/lib/python2.7/site-packages/graph_tool/

Step-by-step installation
-------------------------

From the source code
~~~~~~~~~~~~~~~~~~~~

.. warning:: Not tested.

.. seealso::

    You may follow the official Graph-Tool `compilation procedure
    <http://graph-tool.skewed.de/download#compilation>`__.

There are at least three dependencies to download and install:

* the Boost libraries v1.46 (`homepage <http://www.boost.org/>`__,
  `doc <http://www.boost.org/doc/libs/1_55_0/?view=categorized>`__,
  `download <http://www.boost.org/users/history/version_1_55_0.html>`__),
* expat (`homepage <http://expat.sourceforge.net/>`__,
  `download <http://sourceforge.net/projects/expat/>`__),
* CGAL v4.3 (`homepage <http://www.cgal.org/>`__,
  `doc <http://doc.cgal.org/latest/Manual/index.html>`__,
  `registration form <http://www.cgal.org/download.html>`__).

..
    https://gforge.inria.fr/frs/?group_id=52

Open a shell and type::

    cd $HOME/Downloads
    # get Boost
    wget http://downloads.sourceforge.net/project/boost/boost/1.55.0/boost_1_55_0.tar.gz
    tar xfz boost_1_55_0.tar.gz
    cd boost_1_55_0
    ./bootstrap.sh # --prefix=$HOME/bin/boost/
    ./b2
    ./b2 install
    cd ..
    rm -rf boost_1_55_0/  boost_1_55_0.tar.gz
    # get Expat
    wget http://downloads.sourceforge.net/project/expat/expat/2.1.0/expat-2.1.0.tar.gz
    # get CGAL
    wget https://gforge.inria.fr/frs/download.php/32999/CGAL-4.3-doc_html.tar.gz
    # get Graph-Tool
    wget http://downloads.skewed.de/graph-tool/graph-tool-2.2.28.tar.bz2

.. highlight:: python

