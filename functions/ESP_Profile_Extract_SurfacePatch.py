'''
Created on 04.11.2011

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

class ESP_Profile_Extract_SurfacePatch(object):
    '''
    This class extracts the surface patch, which lies inside the volume of a 
    sphere around the given coordinates, in form a '.dx' file, which then
    contains the esp-values on the vdw surface.
    
    The extracted '.dx' file will be saved with the extension: 
    '*_extract_sp.dx', where '*' stands for the pdb_name without the'.pdb'.
    
    Usage:
        python ESP_Profile_Extract_SurfacePatch.py --pdb=<path> --coord=<[x,y,z]> --sphere-radius=<number>
        
    Additional options:
        --dx=True    if this if True, it looks for the file:
            '*_vdw-PE0.dx' and '*_esp-PE0.dx', where '*' is the pdb_name.
        --meshsize=<meshsize_in_x_y_z>
    '''
    
    usage_warning = 'python ESP_Profile_Extract_SurfacePatch.py --pdb=<path> --coord=<[x,y,z]> --sphere-radius=<number>'
    
    def __init__(self, args):
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
        
        # filename for the extracted surface patch
        extract_dx_filename = '{0}{1}'.format(pdb_path.replace('.pdb', ''),
                                             '_extract_sp.dx')
        
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
        # read esp and vdw grids
        #=======================================================================
        dx_esp = DXReader().parse(pdb_path[:-4] + "_esp-PE0.dx",
                                  DXReader().ESP, mesh_size)
        espbox = dx_esp.getBox()
        dx_vdw = DXReader().parse(pdb_path[:-4] + "_vdw-PE0.dx",
                                  DXReader().VDW, mesh_size)
        vdwbox = dx_vdw.getBox()
        
        # flood the vdw box, to remove holes etc.
        vdwbox.flood()
        
        # find surface
        vdwbox.find_surface()
        
        # collect all points on the surface
        surf_x, surf_y, surf_z = np.nonzero(vdwbox.box == vdwbox.score_of_surface)
        
        # transform the coord_real to the box space
        coord_box = vdwbox.transform_real_to_box_space(coord)
        
        # subtract coord_box from all surface points
        diff_x = surf_x - coord_box[0]
        diff_y = surf_y - coord_box[1]
        diff_z = surf_z - coord_box[2]
        
        # calculate length for each vector
        distance = np.sqrt(diff_x * diff_x + diff_y * diff_y + diff_z * diff_z)
        
        # find all elements which distance is smaller than sphere_radius, with
        # these elements the surface points in the sphere can be accessed
        sphere_points = np.nonzero(distance <= sphere_radius)[0]
        
        # coordinates, which lie inside the sphere
        box_x = surf_x[sphere_points]
        box_y = surf_y[sphere_points]
        box_z = surf_z[sphere_points]
        
        new_box = np.zeros(espbox.box.shape)
        new_box[box_x, box_y, box_z] = espbox.box[box_x, box_y, box_z]
        
        DXWriter().writeBox(extract_dx_filename, DXBox(new_box,
                                                       espbox.getDimensions(),
                                                       espbox.getMeshSize(),
                                                       espbox.getOffset()))
        
if __name__ == '__main__':
    ESP_Profile_Extract_SurfacePatch(sys.argv[1:])
