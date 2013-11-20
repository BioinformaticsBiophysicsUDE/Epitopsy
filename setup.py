from distutils.core import setup

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
)
