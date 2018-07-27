
# Epitopsy

Python package for FFT-driven protein-ligand energy scanning and Direct Coupling Analysis.

## Getting started

### Dependencies

* python >= 2.7.4 (python3.* has not been tested)
* python-dev (headers are required for compilation)
* cython >= 0.28.2
* numpy >= 1.7.0
* biopython >= 1.6.0
* scipy >= 0.19
* sphinx >= 1.7
* [anfft](https://code.google.com/archive/p/anfft/) 2.0 (a manual installation is required)

Command line:

```bash
pip2 install numpy
pip2 install cython
pip2 install scipy
pip2 install biopython
pip2 install sphinx
apt-get install python2.7-dev
```

### Installation

```bash
git clone https://github.com/BioinformaticsBiophysicsUDE/Epitopsy
cd Epitopsy
python2 setup.py install
```

### Unit testing

Tests are carried out within the `unittest` framework:

```bash
cd Epitopsy/epitopsy/epitopsy_unittests/
python2 unittest_suite.py
```

### Documentation

Sphinx can generate the documentation:

```bash
cd Epitopsy/docs
make html
```

Then open `Epitopsy/docs/build/html/index.html` in a web browser.

## Citation

* J.-N. Grad, A. Gigante, C. Wilms, J. N. Dybowski, L. Ohl, C. Ottmann, C. Schmuck and D. Hoffmann, Locating Large, Flexible Ligands on Proteins, *J. Chem. Inf. Model.* **2018**, 58(2), 315-327. DOI: [10.1021/acs.jcim.7b00413](https://doi.org/10.1021/acs.jcim.7b00413)

## Authors

* **Jan Nikolaj Dybowski** - *initial work*
* **Christoph Wilms** - *DCA and FFT modules, porting to Python*
* **Ludwig Ohl** - *DCA module expansion*
* **Jean-Noël Grad** - *FFT module expansion*

## Built with

* [Biopython](https://biopython.org/wiki/Biopython) for protein parsing
* [Cython](http://cython.org/) for faster execution
* [NumPy](https://www.numpy.org/) for most operations on arrays
* [ANFFT](https://code.google.com/archive/p/anfft/) for FFT operations
* [Sphinx](http://sphinx-doc.org/) for documentation

## License

The code is released under the [GNU LGPL v3](https://www.gnu.org/licenses/lgpl.html), unless stated otherwise:

* `Epitopsy/epitopsy/Structure.py`, derived from [Biopython](https://biopython.org/) and released under the [Biopython license](https://github.com/biopython/biopython/blob/master/LICENSE.rst)
* `Epitopsy/docs/_static/copybutton.js`, derived from [CPython](https://github.com/python/cpython) and released under the [Python Software Foundation License](https://github.com/python/cpython/blob/master/LICENSE)

