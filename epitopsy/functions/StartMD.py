'''
Created on May 30, 2012

@author: chris
'''

import sys
import time
from epitopsy.tools.Wrapper import Gromacs

if __name__ == '__main__':
    if len(sys.argv[1:]) != 5:
        print('Missing at least one of the arguments:\n<pdb_path> <simulation_time> <temperature> <new_dir_template> <seeds>')
        sys.exit(1)
    
    pdb_path = sys.argv[1]
    md_time = sys.argv[2]
    temp = sys.argv[3]
    new_dir = sys.argv[4]
    seeds = int(sys.argv[5])
    
    print('PDB Structure:\t{0}'.format(pdb_path))
    print('Simulation time:\t{0} ns'.format(md_time))
    print('Temperature:\t{0} K'.format(temp))
    print('Directory:\t{0}_0:{1}'.format(new_dir, seeds - 1))
    print('#Seeds:\t{0}'.format(seeds))
    
    time.sleep(5)
    
    for i in range(seeds):
        md_dir = '{0}_{1}'.format(new_dir, i)
        if not os.path.exists(md_dir):
            gromacs = Gromacs(pdb_path, md_time, temp, md_dir)
            gromacs.run()
        else:
            continue
