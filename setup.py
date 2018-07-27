#from setuptools import setup
from setuptools import setup, Extension
from numpy import get_include

# get long description from README file
try:
    import pypandoc
    longdesc = pypandoc.convert('README.md', u'rst', format=u'markdown_github')
    longdesc_type = 'text/x-rst'
except ImportError:
    with open('README.md') as f:
        longdesc = f.read()
    longdesc_type = 'text/markdown'

### Python and Numpy Header paths
cyHeaders = [get_include()]

### Cython Files that need compiling
cyFileNames = [
    'epitopsy/cython/interaction_explicit_sampling',
    'epitopsy/cython/dx_cython',
    'epitopsy/cython/pdb_cython',
    'epitopsy/cython/FFTCorrelationScoringCython',
    'epitopsy/cython/optimizer_cython',
    'epitopsy/cython/fast_dca',
    'epitopsy/cython/fast_dca_float32'
]
XtraCompileArgs = ['-shared','-pthread','-fPIC','-fwrapv','-O2','-Wall','-fopenmp','-fno-strict-aliasing']
XtraLinkArgs =    ['-shared','-pthread','-fPIC','-fwrapv','-O2','-Wall','-fopenmp','-fno-strict-aliasing']
cythonExtensions = [Extension(cyFile.replace('/', '.'), [cyFile + '.pyx'],
                    include_dirs=cyHeaders, extra_compile_args=XtraCompileArgs,
                    extra_link_args=XtraLinkArgs, language='c')
                    for cyFile in cyFileNames]

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
    long_description_content_type=longdesc_type,
    url='https://github.com/BioinformaticsBiophysicsUDE/Epitopsy',
    requires=(
        'python (>= 2.7.4)',
        'numpy (>= 1.7.0)',
        'cython (>= 0.28.2)',
        'biopython (>= 1.6.0)'
    ),
    ext_modules = cythonExtensions,
    classifiers=(
        'Programming Language :: Python :: 2',
    ),
)
