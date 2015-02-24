'''
Created on Sep 30, 2011
Edited on Nov 14, 2013
@author: Christoph Wilms, Ludwig Ohl
'''
import os
from subprocess import call
import numpy
from distutils.sysconfig import get_python_inc

import epitopsy.cython

if __name__ == '__main__':
    """
    This scripts builds all extensions that have been built with cython.
    """
    python_version = '2'  # 2 or 3

    cython_folder = epitopsy.cython.__path__[0]
    file_list = filter(lambda x: x.endswith(".pyx"),
                       os.listdir(cython_folder))

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

    python_headers = get_python_inc()
    numpy_headers = numpy.get_include()

    compiler.append('-I{0}'.format(python_headers))
    compiler.append('-I{0}'.format(numpy_headers))
    compiler.append('-o')

    wd = os.getcwd()
    os.chdir(cython_folder)
    for file_path in file_list:
        file_template = file_path.rstrip(".pyx")
        pyx_path = file_path
        c_path = "{0}.c".format(file_template)
        so_path = "{0}.so".format(file_template)
        if os.path.exists(c_path):
            os.remove(c_path)
        if os.path.exists(so_path):
            os.remove(so_path)
        print('#### compiling {0}'.format(file_template))
        call(['cython', "-a", pyx_path])
        call(compiler + [so_path, c_path])
