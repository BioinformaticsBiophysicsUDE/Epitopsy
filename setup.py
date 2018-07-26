#from setuptools import setup
from distutils.core import setup
from Cython.Build import cythonize
from Cython.Distutils import build_ext
from distutils.extension import Extension
#from distutils.sysconfig import get_python_inc
from numpy import get_include

# get long description from README file
try:
    import pypandoc
    longdesc = pypandoc.convert('README.md', u'rst', format=u'markdown_github')
except ImportError:
    with open('README.md') as f:
        longdesc = f.read()

# setuptools DWIM monkey-patch madness
# http://mail.python.org/pipermail/distutils-sig/2007-September/thread.html#8204
import sys
if 'setuptools.extension' in sys.modules:
    m = sys.modules['setuptools.extension']
    m.Extension.__dict__ = m._Extension.__dict__

### Python and Numpy Header paths
cyHeaders = [get_include()]

### Cython Files that need compiling
cyFileNames = [
    'cython/interaction_explicit_sampling',
    'cython/dx_cython', 
    'cython/pdb_cython', 
    'cython/FFTCorrelationScoringCython', 
    'cython/optimizer_cython', 
    'cython/fast_dca', 
    'cython/fast_dca_float32'
]
XtraCompileArgs = ['-shared','-pthread','-fPIC','-fwrapv','-O2','-Wall','-fopenmp','-fno-strict-aliasing']
XtraLinkArgs =    ['-shared','-pthread','-fPIC','-fwrapv','-O2','-Wall','-fopenmp','-fno-strict-aliasing']
cythonExtensions = [Extension('epitopsy/'+cyFile, ['epitopsy/'+cyFile+'.pyx'],
                    include_dirs=cyHeaders, extra_compile_args=XtraCompileArgs,
                    extra_link_args=XtraLinkArgs) for cyFile in cyFileNames]


setup(
    name='epitopsy',
    version='0.1.0',
    author=('Christoph Wilms', 'Ludwig Ohl', 'Jean-Noel Grad'),
    author_email='Ludwig.Ohl@uni-due.de',
    packages=(
        'epitopsy',
        'epitopsy.tools',
        'epitopsy.Optimization',
        'epitopsy.functions',
        'epitopsy.cython',
        'epitopsy.EnergyGrid'
    ),
    license='LICENSE.txt',
    description='Bioinformatics Toolkit.',
    long_description=longdesc,
    url='https://github.com/BioinformaticsBiophysicsUDE/Epitopsy',
    requires=(
        'python (>= 2.7.4)',
        'numpy (>= 1.7.0)',
        'cython (>= 0.18.0)',
        'biopython (>= 1.6.0)'
    ),
    ext_modules = cythonize(cythonExtensions),
    classifiers=(
        'Programming Language :: Python :: 2',
    ),
)
