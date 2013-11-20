# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 08:12:34 2013

@author: Christoph Wilms
"""
import os
import sys
import subprocess

import epitopsy.tools
from epitopsy.cython.fast_dca_float32 import dca_new, dca_old, mi


def run_dca_old(aln_path, di_path, identity_threshold):
    algo = "dca_old"
    path = epitopsy.tools.__path__[0]
    script_path = os.path.join(path, "run_dca.py")
    p = subprocess.call(["python", script_path,
                         "algorithm={0}".format(algo),
                         "in={0}".format(aln_path),
                         "out={0}".format(di_path),
                         "x={0}".format(identity_threshold)])

def run_dca_new(aln_path, di_path, identity_threshold):
    algo = "dca_new"
    path = epitopsy.tools.__path__[0]
    script_path = os.path.join(path, "run_dca.py")
    p = subprocess.call(["python", script_path,
                         "algorithm={0}".format(algo),
                         "in={0}".format(aln_path),
                         "out={0}".format(di_path),
                         "x={0}".format(identity_threshold)])

def run_mi(aln_path, di_path):
    algo = "mi"
    path = epitopsy.tools.__path__[0]
    script_path = os.path.join(path, "run_dca.py")
    p = subprocess.call(["python", script_path,
                         "algorithm={0}".format(algo),
                         "in={0}".format(aln_path),
                         "out={0}".format(di_path)])


## settings
if __name__ == '__main__':
    algorithm = None
    aln_path = None
    di_path = None
    identity_threshold = None
    key_algorithm = "algorithm"
    key_in = "in"
    key_out = "out"
    key_x = "x"
    usage = "python run_dca.py {algo}=<dca_new or dca_old or mi> {aln}=<input_aln_path> {out}=<output_path> {x}=<identity_threshold>".format(
        algo=key_algorithm,
        aln=key_in,
        out=key_out,
        x=key_x)

    for item in sys.argv[1:]:
        item_key = item.split('=')[0]
        item_value = item.split('=')[1]
        if item_key == key_algorithm:
            algorithm = item_value
        elif item_key == key_in:
            aln_path = item_value
        elif item_key == key_out:
            di_path = item_value
        elif item_key == key_x:
            identity_threshold = float(item_value)

    if(algorithm is None or aln_path is None or di_path is None):
        print(usage)
        sys.exit(1)

    if identity_threshold is None:
        if algorithm != "mi":
            print(usage)
            sys.exit(1)

    if algorithm == "dca_new":
        dca_new(aln_path, di_path, identity_threshold=identity_threshold)
    elif algorithm == "dca_old":
        dca_old(aln_path, di_path, identity_threshold=identity_threshold)
    elif algorithm == "mi":
        mi(aln_path, di_path)
    else:
        raise ValueError("Unkown input for key word 'algorithm', can be 'dca_new', 'dca_old' or 'mi'!")
