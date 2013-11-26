Requirements:
The following software packages and versions are required:
python >= 2.7.4 (python3.* has to be tested)
+ python-dev (headers are required for compilation)
cython >= 0.18
numpy >= 1.7.0
biopython >= 1.6.0
anfft (https://code.google.com/p/anfft/) and include it in you pythonpath 

use pip to install python packages

    sudo pip install numpy
    sudo pip install cython
    sudo pip install scipy
    sudo pip install pandas
    sudo pip install biopython 

    install numba (this is a little difficult, but this is a really cool package)

Alternatively install the anaconda package, where most packages are already included:

1. Get right version from: http://continuum.io/downloads 
2. $ bash <downloaded file>
3. conda install biopython 
4. Go to epitopsy directory
5. $ python setup.py install
