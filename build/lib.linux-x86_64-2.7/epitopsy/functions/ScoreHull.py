'''
Created on Sep 29, 2011

@author: Christoph Wilms
'''

import sys
import os
import time
import numpy as np

from epitopsy.Structure import PDBFile
from epitopsy.APBS import APBSWrapper
from epitopsy.tools.UtilityClasses import ESP_Profile_Manager

class ScoreHull(object):
    '''
    classdocs
    '''


    def __init__(self, pdb_name, hull_path, box_dim, box_center,
                 mesh_size = [1, 1, 1], save_hull = False,
                 atom_coord = False, clean_up = False,
                 score_method = 'diff-vec' , print_result = False):
        '''
        This class scores the given profile for the given pdb structure.
        To compare structures it is necessary to provide the target structure, 
        so that the infile with the apbs parameters can be constructed from 
        this structure, so that the offset is the same. Furthermore it is
        necessary to provide the '--extend=<number>' parameter, that has been
        used with the target structure. If none has been used, there is no need
        to use it here either!
        If the option '--atom-coord=yes' is given, the hull is saved as 
        atomic coordinates instead of the array indices.
        The hull of the current protein can be saved if neccessary by
        '--out=<path>'.
        
        Usage:
            python ScoreHull.py --hullpdb=<path> --hull=<path> --target=<path> 
                                        
        optional arguments:
                --meshsize=<number>          default=1
                --out=<path>                 default=no
                --atom-coord=yes             default=no
                --clean-up=yes               default=no 
                --score=diff-vec             default=diff-vec
                                             abs-vec
                       =no
        '''
        score = None
        self.print_result = print_result
        
        pdb = PDBFile(pdb_name)
        
        apbs = APBSWrapper("/usr/bin/apbs", "/usr/bin/pdb2pqr")

        box_dict = apbs.get_dxbox(pdb_path = pdb_name,
                                  mesh_size = mesh_size,
                                  pqr_path = None, box_dim = box_dim,
                                  box_center = box_center, extend = None,
                                  box_type = ['esp'], cubic_box = None,
                                  close_boundaries = False)
            
        esp_box = box_dict['esp']
                                  
        hull_coord, target_phi = ESP_Profile_Manager().read_hull_coordinates(hull_path, True)
        
        try:
            current_phi = esp_box.box[hull_coord]
        except:
            print('Given hull coordinates extend the current grid size! Existing ...')
            sys.exit(1)
        
        # clean up
        if clean_up == True:
            os.remove(pdb.structure_path.replace('.pdb', '.pqr'))
            os.remove(pdb.structure_path.replace('.pdb', '.in'))
            #os.remove(pdb.structure_path.replace('.pdb', '.propka'))
            os.remove(pdb.structure_path.replace('.pdb', '_esp-PE0.dx'))
            #os.remove(pdb.structure_path.replace('.pdb', '_vdw-PE0.dx'))
        
        if score_method != 'no':
            if score_method == 'diff-vec':
                self.score = self.score_diff_vec(target_phi, current_phi)
            elif score_method == 'abs-vec':
                self.score = self.score_abs_vec(target_phi, current_phi)
                
    def score_diff_vec(self, vec1, vec2):      
        # works best
        diff_vec = vec1 - vec2
        score_diff = np.sqrt(np.dot(diff_vec, diff_vec))
        if self.print_result is True:
            print(score_diff)
        
        return score_diff
    
    def score_abs_vec(self, vec1, vec2):
        vec_abs_1 = np.abs(vec1)
        vec_abs_2 = np.abs(vec2)
        abs_vec = vec_abs_1 * vec_abs_2
        score_abs = np.sum(abs_vec) / len(abs_vec)
        if self.print_result is True:
            print(score_abs)
            
        return score_abs 

    def get_score(self):
        return self.score
        
if __name__ == '__main__':
    arguments = {}
    for item in sys.argv[1:]:
        arguments[item.split('=')[0]] = item.split('=')[1]
    if '--hullpdb' in arguments:
        pdb_name = arguments['--hullpdb']
    else:
        print("Usage:\npython ScoreHull.py --hullpdb=<path> --hull=<path> --target=<path>")
        sys.exit(1)
    if '--target' in arguments:
        target_pdb_name = arguments['--target']
    else:
        print("Usage:\npython ScoreHull.py --hullpdb=<path> --hull=<path> --target=<path>")
        sys.exit(1)
    if '--out' in arguments:
        save_hull = True
        out_file_name = arguments['--out']
    else:
        save_hull = False
    if '--extend' in arguments:
        extend = int(arguments['--extend'])
    else:
        extend = 6
    if '--hull' in arguments:
        hull_path = arguments['--hull']
    else:
        print("Usage:\npython ScoreHull.py --hullpdb=<path> --hull=<path> --target=<path>")
        sys.exit(1)
    if '--atom-coord' in arguments:
        if arguments['--atom-coord'] == 'yes' or arguments['--atom-coord'] == 'y': 
            atom_coord = True
        else:
            atom_coord = False
    else:
        atom_coord = False
        
    if '--score' in arguments:
        if arguments['--score'] == 'diff-vec':
            score_method = arguments['--score']
        elif arguments['--score'] == 'no' or arguments['--score'] == 'n':
            score_method = 'no'
        else:
            print('Unrecognized scoring method! Existing ...')
            sys.exit(1)
    else:
        score_method = 'diff-vec'
    if '--meshsize' in arguments:
        m = float(arguments['--meshsize'])
    else:
        m = 1
    if '--clean-up' in arguments:
        if arguments['--clean-up'] == 'yes' or arguments['--clean-up'] == 'y': 
            clean_up = True
        else:
            clean_up = False
    else:
        clean_up = False
    ScoreHull(pdb_name, hull_path, target_pdb_name, [m, m, m],
              extend, save_hull, atom_coord, clean_up, score_method, True)
