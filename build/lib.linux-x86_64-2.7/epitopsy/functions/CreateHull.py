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

def Create_Hull(pdb_name, out_file_name, mesh_size, extend = 6,
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
    starttime = time.time()
    
    pdb = PDBFile(pdb_name)
    
    apbs = APBSWrapper("/usr/bin/apbs", "/usr/bin/pdb2pqr")

    box_dict = apbs.get_dxbox(pdb_path = pdb_name,
                              mesh_size = mesh_size,
                              pqr_path = None, box_dim = None,
                              box_center = None, extend = extend,
                              box_type = ['esp', 'vdw'], cubic_box = None,
                              close_boundaries = False)
            
    esp_box = box_dict['esp']
    
    vdw_box = box_dict['vdw']
    
    # account for meshsize
    extend = np.ceil(extend / mesh_size[0])
    
    if extend > 0:
        # mesh_size[0] = mesh_size[1] = mesh_size[2]
        vdw_box.extendSurface(extend)
    
    vdw_box.flood()
    
    vdw_box.find_surface()
    
    # get phi_values
    hull_coordinates = np.nonzero(vdw_box.box == vdw_box.score_of_surface)
    
    ESP_Profile_Manager().write_profile(esp_box.box, hull_coordinates,
                                        1, out_file_name, False) 
    
    # clean up
    if clean_up == True:
        os.remove(pdb.structure_path.replace('.pdb', '.pqr'))
        os.remove(pdb.structure_path.replace('.pdb', '.in'))
        #os.remove(pdb.structure_path.replace('.pdb', '.propka'))
    
    time_for_execution = time.time() - starttime
    #print('time spent: {0}s'.format(time_for_execution))

class CreateHull2(object):
    '''
    classdocs
    '''
    warn_usage = 'Usage:\npython CreateHull.py --pdb_top=<path> --pdb_bottom --out=<path>'  

    def __init__(self, args):
        '''
        This class calculates the electrostatic hull for one given protein 
        structure and saves it under the given path. 
        The hull can be extended by the command '--extend=<number>' to account
        for structural changes if it is compared to some other profile.
        If the option '--atom-coord=yes' is given, the hull is saved with 
        atomic coordinates instead of the array indices.
        
        Usage:
            python CreateHull.py --pdb_top=<path> --pdb_bottom --out=<path>
                                        
        optional arguments:
                --meshsize=<number>          default=1
                --extend=<number_in_A>       default=6 (Nikos suggestion)
                --atom-coord=yes             default=no
                --clean-up=yes               default=no 
        '''
        arguments = {}
        for item in args:
            arguments[item.split('=')[0]] = item.split('=')[1]
        if '--pdb_top' in arguments:
            pdb_top_name = arguments['--pdb_top']
        else:
            print(self.warn_usage)
            sys.exit(1)
        if '--pdb_bottom' in arguments:
            pdb_bottom_name = arguments['--pdb_bottom']
        else:
            print(self.warn_usage)
            sys.exit(1)
        if '--out' in arguments:
            out_file_name = arguments['--out']
        else:
            print(self.warn_usage)
            sys.exit(1)
        if '--extend' in arguments:
            # extend hull by 6 Angstroem
            extend = int(arguments['--extend'])
        else:
            extend = 6
        if '--atom-coord' in arguments:
            if arguments['--atom-coord'] == 'yes' or arguments['--atom-coord'] == 'y': 
                atom_coord = True
            else:
                atom_coord = False
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
        
        mesh_size = [m, m, m]
        
        # account for meshsize
        extend = int(extend / mesh_size[0])
        
        starttime = time.time()
        
        pdb_top = PDBFile(pdb_top_name)
        pdb_bottom = PDBFile(pdb_bottom_name)
        
        # run apbs
        self.run_apbs(pdb_top, extend, mesh_size)
        self.run_apbs(pdb_bottom, extend, mesh_size)
        
        esp_top_box = DXReader().parse(pdb_top_name[:-4] + "_esp-PE0.dx", DXReader().ESP,
                                  mesh_size)
        
        vdw_top_box = DXReader().parse(pdb_top_name[:-4] + "_vdw-PE0.dx", DXReader().VDW,
                                  mesh_size)
        
        esp_bottom_box = DXReader().parse(pdb_bottom_name[:-4] + "_esp-PE0.dx",
                                          DXReader().ESP)
        
        vdw_bottom_box = DXReader().parse(pdb_bottom_name[:-4] + "_vdw-PE0.dx",
                                          DXReader().VDW)
        
        # take top of esp_top_box and put it on the bottom of esp_bottom_box
        # cut x 
        # top: modelneg_repos
        # bottom: 1LP1_repos
        box_x = np.zeros(esp_bottom_box.box.shape)
        box_x[0:35][:, :] = esp_bottom_box.box[0:35][:, :]
        box_x[35:][:, :] = esp_top_box.box[35:][:, :]
        
        box_x1 = esp_bottom_box.box[0:35][:, :]
        box_x2 = esp_top_box.box[35:][:, :]
        
        if extend > 0:
            # mesh_size[0] = mesh_size[1] = mesh_size[2]
            vdw_top_box.extendSurface(extend / mesh_size[0])
        
        vdw_top_box.flood()
        
        vdw_top_box.find_surface()
        
        # get phi_values
        hull_coordinates = np.nonzero(vdw_top_box.box == vdw_top_box.score_of_surface)
        
        ESP_Profile_Manager().write_profile(box_x, hull_coordinates,
                                            1, out_file_name, False) 
        
        # clean up
        if clean_up == True:
            os.remove(pdb_top.PDBFilename[:-4] + '.pqr')
            os.remove(pdb_top.PDBFilename[:-4] + '.in')
        
        time_for_execution = time.time() - starttime
        #print('time spent: {0}s'.format(time_for_execution))

    def run_apbs(self, pdb, extend, mesh_size):
        apbs = APBSWrapper("/usr/bin/apbs", "/usr/bin/pdb2pqr")

        apbs.runPDB2PQR(pdb.structure_path,
                        pdb.structure_path[:-4] + '.pqr')
                
        # set padding to extend
        padding = extend
        
        template_in = InFile("", "", "", mesh_size)
        
        # create template InFile 
        template_in.generateFromPDB(pdb, padding, True)
        
        template_in.setPQRFilePath(pdb.structure_path[:-4] + '.pqr')
        template_in.setOutSurfaceDXPath(pdb.structure_path[:-4] + "_vdw")
        template_in.setOutPotentialDXPath(pdb.structure_path[:-4] + "_esp")
        apbs.runAPBS(template_in, pdb.structure_path[:-4] + '.in')

if __name__ == '__main__':
    arguments = {}
    for item in sys.argv[1:]:
        arguments[item.split('=')[0]] = item.split('=')[1]
    if '--pdb' in arguments:
        pdb_name = arguments['--pdb']
    else:
        print("Usage:\n\tpython CalculateESPProfiles --pdb=<path> --out=<path>")
        sys.exit(1)
    if '--out' in arguments:
        out_file_name = arguments['--out']
    else:
        print("Usage:\n\tpython CalculateESPProfiles --pdb=<path> --out=<path>")
        sys.exit(1)
    if '--extend' in arguments:
        # extend hull by 6 Angstroem
        extend = int(arguments['--extend'])
    else:
        extend = 6
    if '--atom-coord' in arguments:
        if arguments['--atom-coord'] == 'yes' or arguments['--atom-coord'] == 'y': 
            atom_coord = True
        else:
            atom_coord = False
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
    Create_Hull(pdb_name, out_file_name, [m, m, m], extend, clean_up)
