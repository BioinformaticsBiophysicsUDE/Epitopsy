'''
Created on 03.11.2011

@author: chris
'''

import sys
import numpy as np

from epitopsy.pdb.PDBFile import PDBFile
from epitopsy.apbs.APBSWrapper import APBSWrapper
from epitopsy.apbs.InFile import InFile
from epitopsy.grid.Grid import Grid
from epitopsy.dx.DXBox import DXBox
from epitopsy.dx.DXReader import DXReader
from epitopsy.dx.DXWriter import DXWriter

class ESP_Profile_SurfacePatch(object):
    '''
    This class calculates the ESP-Profie for a patch on the surface of the 
    given protein. The patch is calculated like this:
        1) the vdw surface of the protein is extended by the supplied value
        2) a sphere with the given value is created around the supplied
            coordinates
        3) all points on the protein surface which lie inside the sphere
            make up the surface patch
    
    This patch is then written to a file for further analysis.
    
    The input pdb should only contain the proteins, from which the esp hull
    is of interest. If one wants to know the esp-profile of a binding site
    for a ligand, that ligand has to be removed!!!
    
    Parameters of 'extend' and 'sphere radius' are in Angstroem.
    
    Usage:
        python ESP_Profile_SurfacePatch.py --pdb=<path> --out=<path> --coord=<[x,y,z] or [x  y  z]> --extend=<number> --sphere-radius=<number>
        
    Notice:
        * the best way to supply the '--coord' is to use 'str(list(coord))'
            in a python call (real space)
    
    Additional options:
        --meshsize=<meshsize_in_x_y_z>
    
    '''
    usage_warning = 'python ESP_Profile_SurfacePatch.py --pdb=<path> --out=<path> --coord=<[x,y,z]> --extend=<number> --sphere-radius=<number>'
    
    
    def __init__(self, args):
        arguments = {}
        for item in args:
            arguments[item.split('=')[0]] = item.split('=')[1]
        
        if '--pdb' not in arguments:
            raise AttributeError(self.usage_warning)
        else:
            pdb_path = arguments['--pdb']
        
        if '--out' not in arguments:
            raise AttributeError(self.usage_warning)
        else:
            result_path = arguments['--out']
        
        
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
        
        # accordingly
        # in the following everything is divided by mesh to scale everything
        
        if '--extend' in arguments:
            extend_surface = int(int(arguments['--extend']) / mesh)
        else:
            extend_surface = int(5 / mesh)
            
        if '--sphere-radius' in arguments:
            sphere_radius = int(int(arguments['--sphere-radius']) / mesh)
        else:
            sphere_radius = int(5 / mesh)
        
        
        # load the structure
        pdb = PDBFile(pdb_path)
        
        #=======================================================================
        # run esp calculations
        #=======================================================================
        apbs = APBSWrapper("apbs", "pdb2pqr")
        template_in = InFile("", "", "", mesh_size)
        
        # padd the box by extend_surface + 8, to ensure, that the algorithm
        # will not touch the borders of the box
        padding = int(extend_surface + np.ceil(2 / mesh_size[0]))
        
        # there were some errors with very lengthy proteins
        if padding < 9:
            padding = 9
        
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
        
        # extend the vdw box
        vdwbox.extendSurface(extend_surface)
        
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
        
        # get esp-values
        esp_values = espbox.box[box_x, box_y, box_z]
        
        f = open(result_path, 'a')
        f.write(pdb_path)
        for item in esp_values:
            f.write(' ' + str(item))
        f.write('\n')
        f.close()


        
if __name__ == '__main__':
    ESP_Profile_SurfacePatch(sys.argv[1:])
