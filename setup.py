from setuptools import setup
from Cython.Build import cythonize
from distutils.extension import Extension
from distutils.sysconfig import get_python_inc
import numpy

# setuptools DWIM monkey-patch madness
# http://mail.python.org/pipermail/distutils-sig/2007-September/thread.html#8204
import sys
if 'setuptools.extension' in sys.modules:
    m = sys.modules['setuptools.extension']
    m.Extension.__dict__ = m._Extension.__dict__



setup(
    name='epitopsy',
    version='0.1.0',
    author=['Christoph Wilms', 'Ludwig Ohl'],
    author_email='Ludwig.Ohl@uni-due.de',
    packages=[
		'epitopsy', 
		'epitopsy.tools', 
		'epitopsy.Optimization',
		'epitopsy.functions',
		'epitopsy.cython',
		'epitopsy.scoring',
		'epitopsy.result'
	],
    license='LICENSE.txt',
    description='Bioinforamtics Toolkit.',
    long_description=open('README.txt').read(),
    install_requires=[
	"python >= 2.7.4",
        "numpy >= 1.7.0",
        "cython >= 0.18.0",
	"biopython >= 1.6.0"
    ],
    setup_requires=[
            'setuptools_cython',
    ],
    ext_modules = [
        Extension('epitopsy.cython', ["epitopsy/cython/CalculateInteractionEnergyCython.pyx", "epitopsy/cython/dx_cython.pyx", "epitopsy/cython/pdb_cython.pyx", "epitopsy/cython/FFTCorrelationScoringCython.pyx", "epitopsy/cython/optimizer_cython.pyx", "epitopsy/cython/fast_dca.pyx", "epitopsy/cython/fast_dca_float32.pyx"], include_dirs= [numpy.get_include(), get_python_inc()])
    ],
)
