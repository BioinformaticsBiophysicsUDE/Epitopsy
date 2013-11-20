'''
Created on Sep 30, 2011
Edited on Nov 14, 2013
@author: Christoph Wilms, Ludwig Ohl
'''
import os
import sys
from subprocess import call
import numpy
from distutils.sysconfig import get_python_inc

if __name__ == '__main__':
    """
    This scripts builds all extensions that have been built with cython.
    """
    python_version = '2' # 2 or 3

    # this list contains all files
    file_list = []

    file_list.append('cython/CalculateInteractionEnergyCython')

    file_list.append('cython/dx_cython')

    file_list.append('cython/pdb_cython')

    file_list.append('cython/FFTCorrelationScoringCython')

    file_list.append('cython/optimizer_cython')

    file_list.append('cython/fast_dca')

    file_list.append('cython/fast_dca_float32')

#    file_list.append('cython/test_dca')


    # compiler  + options:
    compiler = ['gcc']
    compiler.append('-shared')
    compiler.append('-pthread')
    compiler.append('-fPIC')
    compiler.append('-fwrapv')
    compiler.append('-O2')
    compiler.append('-Wall')
    compiler.append('-fopenmp')
    compiler.append('-fno-strict-aliasing')
    if python_version == '2':
        python_headers_usr = os.path.join("/usr", "include","python2.7")
        numpy_headers_usr = os.path.join("/usr", "include","python2.7","numpy")

        python_headers_cluster = get_python_inc()
	print 'Using python Headers from ', python_headers_cluster
        numpy_headers_cluster = numpy.get_include()
	print 'Using numpy Headers from ', numpy_headers_cluster

        if(os.path.exists(python_headers_usr)
                and os.path.exists(numpy_headers_usr)):
            python_headers = python_headers_usr
            numpy_headers = numpy_headers_usr

        elif (os.path.exists(python_headers_cluster)
                and os.path.exists(numpy_headers_cluster)):
            python_headers = python_headers_cluster
            numpy_headers = numpy_headers_cluster

        else:
            raise AttributeError("Could not find python and numpy header files!")

#        python_headers_anaconda = "{0}".format(os.path.join(os.getenv("HOME"),
#                                      "Programs",
#                                      "anaconda",
#                                      "include",
#                                      "python2.7"))
#        numpy_headers_anaconda = "-I{0}".format(os.path.join(os.getenv("HOME"),
#                                       "Programs",
#                                       "anaconda",
#                                       "pkgs",
#                                       "numpy-1.7.0-py27_0",
#                                       "lib",
#                                       "python2.7",
#                                       "site-packages",
#                                       "numpy",
#                                       "core",
#                                       "include"))



    else:
        raise ValueError('Unkwon python version: {0}'.format(python_version))

    compiler.append('-I{0}'.format(python_headers))
    compiler.append('-I{0}'.format(numpy_headers))
    compiler.append('-o')
#    compiler.append

    wd = os.getcwd()
    for file_path in file_list:
        sub_wd = os.path.split(file_path)[0]
        sub_path = os.path.split(file_path)[1]
        os.chdir(sub_wd)
        pyx_path = "{0}.pyx".format(sub_path)
        c_path = "{0}.c".format(sub_path)
        so_path = "{0}.so".format(sub_path)
        if os.path.exists(c_path):
            os.remove(c_path)
        if os.path.exists(so_path):
            os.remove(so_path)
        print('#### compiling {0}'.format(sub_path))
        call(['cython', "-a", pyx_path])
        call(compiler + [so_path, c_path])

        os.chdir(wd)
