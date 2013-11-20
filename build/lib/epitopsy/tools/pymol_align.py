import __main__
# quiet and no gui
__main__.pymol_argv = ['pymol','-qc']

import os
import sys
import time
import pymol

pymol.finish_launching()


## settings
usage = "python pymol_align.py <ref.pdb> <query_0.pdb> <...> <query_N.pdb>"

if len(sys.argv[1:]) < 2:
    raise AttributeError(usage)

ref_name = "ref"
ref_path = sys.argv[1]

pymol.cmd.load(ref_path, ref_name)
for i,query_item in enumerate(sys.argv[2:]):
    query_name = "query_{0}".format(i)
    pymol.cmd.load(query_item, query_name)
    pymol.cmd.align(query_name, ref_name)
    pymol.cmd.save(query_item, query_name)
    pymol.cmd.delete(query_name)

pymol.cmd.quit()

