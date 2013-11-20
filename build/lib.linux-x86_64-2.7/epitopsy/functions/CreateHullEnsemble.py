'''
Created on Sep 29, 2011

@author: Christoph Wilms
'''

import sys
import os
import time
import numpy as np
from subprocess import call

from epitopsy.Structure import PDBFile
from epitopsy.APBS import APBSWrapper, InFile
from epitopsy.DXFile import DXReader
from epitopsy.tools.UtilityClasses import ESP_Profile_Manager


def Create_Hull_Ensemble(pdb_name_list, out_file_name, box_dim, box_center,
                         mesh_size, extend_surface = 6,
                         clean_up = False):
    '''
    This class calculates the electrostatic hull for one given protein 
    structure and saves it under the given path. 
    The hull can be extended by the command '--extend=<number>' to account
    for structural changes if it is compared to some other profile.
    If the option '--atom-coord=yes' is given, the hull is saved with 
    atomic coordinates instead of the array indices.
    
    Usage:
        python CreateHull.py --pdb=<path> --out=<path>
                                    
    optional arguments:
            --extend=<number_in_A>       default=6 (Nikos suggestion)
            --atom-coord=yes             default=no
            --clean-up=yes               default=no 
    '''
    average_box = None
    starttime = time.time()
    
    for i, pdb_name in enumerate(pdb_name_list):
        # account for meshsize
        extend_surface = np.ceil(extend_surface / mesh_size[0])
        
        pdb = PDBFile(pdb_name)
        
        apbs = APBSWrapper("/usr/bin/apbs", "/usr/bin/pdb2pqr")
    
        box_dict = apbs.get_dxbox(pdb_path = pdb_name,
                                  mesh_size = mesh_size,
                                  pqr_path = None, box_dim = box_dim,
                                  box_center = box_center, extend = None,
                                  box_type = ['esp', 'vdw'], cubic_box = None,
                                  close_boundaries = False)
            
        esp_box = box_dict['esp']
        
        vdw_box = box_dict['vdw']
        
        if i == 0:
            if extend_surface > 0:
                # mesh_size[0] = mesh_size[1] = mesh_size[2]
                vdw_box.extendSurface(extend_surface)
            
            vdw_box.flood()
            
            vdw_box.find_surface()
            
            average_box = esp_box.box
            # get phi_values
            hull_coordinates = np.nonzero(vdw_box.box == vdw_box.score_of_surface)
        else:
            average_box = average_box + esp_box.box
        
        # clean up
        if clean_up == True:
            os.remove(pdb.structure_path.replace('.pdb', '.pqr'))
            os.remove(pdb.structure_path.replace('.pdb', '.in'))
            #os.remove(pdb.structure_path.replace('.pdb', '.propka'))
            os.remove(pdb_name.replace('.pdb', "_esp-PE0.dx"))
            os.remove(pdb_name.replace('.pdb', "_vdw-PE0.dx"))
    
    # mean energy
    average_box = average_box / len(pdb_name_list)
    
    ESP_Profile_Manager().write_profile(average_box, hull_coordinates,
                                        1, out_file_name, False) 
    
    
    
    time_for_execution = time.time() - starttime
    print('time spent: {0}s'.format(time_for_execution))
