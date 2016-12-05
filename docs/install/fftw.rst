*************************************
FFTW --- Fast Fourier Transform Tools
*************************************

Some features of Epitopsy require a local copy of FFTW v3.3+ (`homepage
<http://www.fftw.org/>`__, `documentation <http://www.fftw.org/fftw3_doc/>`__)
in threaded, double precision to run (namely :mod:`FFTCorrelationScoring`)
as well as one of these Python wrappers:

    * pyFFTW v0.10.1 (`homepage <https://pypi.python.org/pypi/pyFFTW>`__,
      `documentation <http://hgomersall.github.io/pyFFTW/>`__)
    * anfft v0.2 (`homepage <https://code.google.com/p/anfft/>`__),
      last update in 2012
    * PyFFTW3 v0.2.1 (`homepage <https://pypi.python.org/pypi/PyFFTW3>`__),
      last update in 2010
    


.. highlight:: bash

Installation procedure
======================


Install the FFTW v3 package and pyFFTW::

    sudo apt-get install libfftw3-dev
    sudo pip install pyfftw

If package `libfftw3-dev` is unavailable for your distribution, compile it
from sources::

    version="3.3.4"
    mkdir ~/Downloads/fftw; cd $_
    wget -O - http://www.fftw.org/fftw-${version}.tar.gz | tar xfz - -C . --strip-components=1
    ./configure --enable-shared --enable-threads --enable-long-double #--prefix=~/bin/fftw-${version}
    make -j$(nproc)
    sudo make install
    make clean; cd ..; rm -rf fftw
    echo -e "\n# added by $(whoami) on $(date) to source FFTW3" >> ~/.bashrc
    echo 'export FFTW_LIBRARY_PATH=/usr/local/lib' >> ~/.bashrc
    echo 'export FFTW_PATH=/usr/local/lib' >> ~/.bashrc

If pip crashes while compiling pyFFTW, make sure that FFTW was
compiled with `--enable-shared --enable-threads --enable-long-double`.
Perform a manual compilation to get complete error messages::

    version="0.10.1"
    mkdir ~/Downloads/pyfftw; cd $_
    wget -O - https://pypi.python.org/packages/source/p/pyFFTW/pyFFTW-${version}.tar.gz | tar xfz - -C . --strip-components=1
    python setup.py build_ext --inplace # don't install, just read error message

If `--enable-shared` was not set when compiling FFTW:

.. code-block:: none

    /usr/bin/ld: //usr/local/lib/libfftw3.a(export-wisdom.o): relocation R_X86_64_32
    against `.text' can not be used when making a shared object; recompile with -fPIC

If `--enable-threads` or `--enable-long-double` was not set when compiling
FFTW:

.. code-block:: none

    x86_64-linux-gnu-gcc -pthread -shared -Wl,-O1 -Wl,-Bsymbolic-functions -Wl,
    -Bsymbolic-functions -Wl,-z,relro -fno-strict-aliasing -DNDEBUG -g -fwrapv
    -O2 -Wall -Wstrict-prototypes -D_FORTIFY_SOURCE=2 -g -fstack-protector
    --param=ssp-buffer-size=4 -Wformat -Werror=format-security
    build/temp.linux-x86_64-2.7/home/grad/Downloads/pyfftw/pyfftw/pyfftw.o
    -lfftw3 -lfftw3f -lfftw3l -lfftw3_threads -lfftw3f_threads -lfftw3l_threads
    -o build/lib.linux-x86_64-2.7/pyfftw/pyfftw.so
    /usr/bin/ld: cannot find -lfftw3f           <-- missing --enable-long-double
    /usr/bin/ld: cannot find -lfftw3f_threads   <--                  "
    /usr/bin/ld: cannot find -lfftw3l           <--                  "
    /usr/bin/ld: cannot find -lfftw3l_threads   <--                  "
    /usr/bin/ld: cannot find -lfftw3_threads    <-- missing --enable-threads

An alternative to pyFFTW is anfft::

    mkdir ~/Downloads/anfft; cd $_
    wget -O - https://anfft.googlecode.com/files/anfft-0.2.tar.gz | tar xfz - -C . --strip-components=1
    sudo python setup.py install
    python -c "import anfft" # raises an error and creates ~/.local/.anfft
    python -c "import anfft" # raises no error
    cd ..; sudo rm -rf anfft

If pyFFTW or anfft freeze while processing a FFT, make sure environment
variables $FFTW_PATH and $FFTW_LIBRARY_PATH are correctly set.


