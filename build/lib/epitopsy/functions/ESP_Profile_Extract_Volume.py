'''
Created on 21.11.2011

@author: chris
'''

import os
import sys
import numpy as np

from epitopsy.pdb.PDBFile import PDBFile
from epitopsy.apbs.APBSWrapper import APBSWrapper
from epitopsy.apbs.InFile import InFile
from epitopsy.grid.Grid import Grid
from epitopsy.dx.DXBox import DXBox
from epitopsy.dx.DXReader import DXReader
from epitopsy.dx.DXWriter import DXWriter

class ESP_Profile_Extract_Volume(object):
    '''
    classdocs
    '''


    def __init__(self, args):
        '''
        Constructor
        '''
        arguments = {}
        for item in args:
            arguments[item.split('=')[0]] = item.split('=')[1]
        
        if '--pdb' not in arguments:
            raise AttributeError(self.usage_warning)
        else:
            pdb_path = arguments['--pdb']
        
        if '--coord' in arguments:
            # remove brackets and spaces etc.
            coord_str = arguments['--coord']
            coord_str = coord_str.lstrip('[')
            coord_str = coord_str.rstrip(']')
            coord_str = coord_str.strip()
            
            # only works with lists
            coord_str_list = coord_str.split(',')
                
            # transform str to numbers
            for i, item in enumerate(coord_str_list):
                coord_str_list[i] = float(item)
            
            # make an array from the list
            coord = np.array(coord_str_list)
            
        else:
            raise AttributeError(self.usage_warning)
        
        if '--meshsize' in arguments:
            mesh = float(arguments['--meshsize'])
            mesh_size = [mesh, mesh, mesh]
        else:
            mesh = 1.0
            mesh_size = [mesh, mesh, mesh]
        
        if '--dx' in arguments:
            if arguments['--dx'] == 'True':
                if(os.path.exists('{0}{1}'.format(pdb_path.replace('.pdb', ''),
                                                  '_esp-PE0.dx')) and
                   os.path.exists('{0}{1}'.format(pdb_path.replace('.pdb', ''),
                                                  '_vdw-PE0.dx'))):
                    apbs_calculation_done = True
                else:
                    apbs_calculation_done = False
            else:
                apbs_calculation_done = False
        
        # in the following everything is divided by mesh to scale everything
        # accordingly
        
        if '--sphere-radius' in arguments:
            sphere_radius = int(int(arguments['--sphere-radius']) / mesh)
        else:
            sphere_radius = int(5 / mesh)
        
        # filename for the extracted consensus volume
        extract_dx_filename = '{0}{1}'.format(pdb_path.replace('.pdb', ''),
                                             '_extract_vol.dx')
        
        #=======================================================================
        # run apbs if necessary
        #=======================================================================
        if apbs_calculation_done is False:
            # load the structure
            pdb = PDBFile(pdb_path)
            
            #=======================================================================
            # run esp calculations
            #=======================================================================
            apbs = APBSWrapper("apbs", "pdb2pqr")
            template_in = InFile("", "", "", mesh_size)
            
            # padd the box by 2, to ensure, that the algorithm
            # will not touch the borders of the box
            padding = 2
            
            template_in.generateFromPDB(pdb, padding, True)
            
            # run pdb2pqr
            apbs.runPDB2PQR(pdb_path, pdb_path[:-4] + '.pqr')
            
            template_in.setPQRFilePath(pdb_path[:-4] + '.pqr')
            template_in.setOutSurfaceDXPath(pdb_path[:-4] + "_vdw")
            template_in.setOutPotentialDXPath(pdb_path[:-4] + "_esp")
            apbs.runAPBS(template_in, pdb_path[:-4] + '.in')
          
        #=======================================================================
        # read esp grid
        #=======================================================================
        dx_esp = DXReader().parse(pdb_path[:-4] + "_esp-PE0.dx",
                                  DXReader().ESP, mesh_size)
        espbox = dx_esp.getBox()
        
        # transform the coord_real to the box space
        coord_box = espbox.transform_real_to_box_space(coord)
        
        # make a sphere at (0,0,0)
        box_indices = ([], [], [])
        esp_indices = ([], [], [])

        for x in range(-sphere_radius, sphere_radius + 1):
            for y in range(-sphere_radius, sphere_radius + 1):
                for z in range(-sphere_radius, sphere_radius + 1):
                    if np.sqrt(x ** 2 + y ** 2 + z ** 2) < sphere_radius:
                        box_indices[0].append(x + sphere_radius)
                        box_indices[1].append(y + sphere_radius)
                        box_indices[2].append(z + sphere_radius)
                        esp_indices[0].append(x + coord_box[0])
                        esp_indices[1].append(y + coord_box[1])
                        esp_indices[2].append(z + coord_box[2])
                        
        
        # new shape is reduced to the sphere, +1 because of (0,0,0)
        new_shape = [2 * sphere_radius + 1, 2 * sphere_radius + 1,
                     2 * sphere_radius + 1]
                
        new_box = np.zeros(new_shape)
        # add sphere_radius, so that the sphere is centered at the middle of
        # array
        # add coord_box, this is where the P is
        new_box[box_indices] = espbox.box[esp_indices]
        
        DXWriter().writeBox(extract_dx_filename, DXBox(new_box,
                                                       new_shape,
                                                       mesh_size,
                                                       [-sphere_radius, -sphere_radius, -sphere_radius]))
        print('wrote {0}'.format(extract_dx_filename))
        
if __name__ == '__main__':
    ESP_Profile_Extract_Consensus_Volume(sys.argv[1:])
