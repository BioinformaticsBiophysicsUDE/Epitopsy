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

class ScoreHullEnsemble(object):
    '''
    classdocs
    '''


    def __init__(self, pdb_name_list, hull_path, box_dim, box_center,
                 mesh_size = [1, 1, 1], save_hull = False,
                 atom_coord = False, clean_up = False,
                 score_method = 'diff-vec' , print_result = False):
        '''
        TODO: This text needs an update!!!
        
        
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
        np.seterr('raise')
        self.number_of_calculations = len(pdb_name_list)
        
        hull_coord, target_phi = ESP_Profile_Manager().read_hull_coordinates(hull_path, True)
        
        mean_phi = np.zeros(target_phi.shape)
        
        for pdb_name in pdb_name_list:

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
            
            try:
                current_phi = esp_box.box[hull_coord]
            except:
                print('Given hull coordinates extend the current grid size! Existing ...')
                sys.exit(1)
            
            mean_phi = mean_phi + current_phi
            # clean up
            if clean_up == True:
                os.remove(pdb.structure_path.replace('.pdb', '.pqr'))
                os.remove(pdb.structure_path.replace('.pdb', '.in'))
                #os.remove(pdb.structure_path.replace('.pdb', '.propka'))
                os.remove(pdb.structure_path.replace('.pdb', '_esp-PE0.dx'))
                #os.remove(pdb.structure_path.replace('.pdb', '_vdw-PE0.dx'))
            
        # normalize 
        mean_phi = mean_phi / self.number_of_calculations
        
        if score_method != 'no':
            if score_method == 'diff-vec':
                self.score = self.score_diff_vec(target_phi, mean_phi)
            elif score_method == 'abs-vec':
                self.score = self.score_abs_vec(target_phi, mean_phi)
                
    def score_diff_vec(self, vec1, vec2):      
        # works best
        diff_vec = vec1 - vec2
        score_diff = np.sqrt(np.dot(diff_vec, diff_vec))
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
        # normalize
        return self.score 
        
