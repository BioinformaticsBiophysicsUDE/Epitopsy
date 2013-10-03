"""
Created on Aug 2, 2011

@author: chris
"""

from subprocess import call
import sys
import os
import time
from numpy import nonzero, ceil


from epitopsy.pdb.PDBFile import PDBFile
from epitopsy.apbs.InFile import InFile
from epitopsy.apbs.APBSWrapper import APBSWrapper
from epitopsy.grid.Grid import Grid
from epitopsy.dx.DXBox import DXBox
from epitopsy.dx.DXReader import DXReader
from epitopsy.dx.DXWriter import DXWriter
from epitopsy.tools.UtilityClasses import ESP_Profile_Manager

class CalculateESPProfiles:
    """
    classdocs
    """


    def __init__(self, args):
        """
        This class reads pdb-names from a given list and calculates the 
        electrostatic potential for each structure. Then it extends the surface
        until all structures are included. Therefore the structures should be 
        already aligned!
        Once the optimal shell is found it extracts the electrostatic potential
        at these points for all structures and writes them to a file.
        The template for the optimal shell is the largest structure. 
        
        With the option '--no-clean-up' all '.pqr' and  '.in' files are kept, 
        otherwise they will be deleted.
        
        Usage:
            python CalculateESPProfiles --list=<path> --out=<path>
            
        optional arguments:
                --extend=
                    biggest    ->             extends the hull of the biggest
                                                structure until all structure
                                                fit into this hull
                    all        -> default!    finds a hull that covers all 
                                                structures and extends it by 
                                                one mesh unit
                --hull-coordinates=<path>
                --meshsize=<number>
                --clean-up=no
                
        """
        arguments = {}
        for item in args:
            arguments[item.split('=')[0]] = item.split('=')[1]
        if '--list' in arguments:
            pdb_list_name = arguments['--list']
        else:
            print("Usage:\n\tpython CalculateESPProfiles --list=<path> --out=<path>")
            sys.exit(1)
        if '--out' in arguments:
            out_file_name = arguments['--out']
        else:
            print("Usage:\n\tpython CalculateESPProfiles --list=<path> --out=<path>")
            sys.exit(1)
        if '--extend' in arguments:
            if (arguments['--extend'] != 'all' or 
                arguments['--extend'] != 'biggest'):
                print("Option '--extend' can either be all or biggest!")
                sys.exit(1) 
            extend_mode = arguments['--extend']
        else:
            extend_mode = 'all'
        if '--hull-coordinates' in arguments:
            hull_coordinates_exist = True
            hull_coordinates_file = arguments['--hull-coordinates']
        else:
            hull_coordinates_exist = False
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
        starttime = time.time()
        
        # TODO: do structural alignment
        
        # read and load pdb structures        
        pdb_list = self.read_pdb_list(pdb_list_name)
        
        apbs = APBSWrapper("/usr/bin/apbs", "/usr/bin/pdb2pqr")
        
        # find the biggest structure, this is then used as a "template"
        # calculate pqr files as well
        max_diameter = 0
        max_diameter_id = 0
        for i, pdb in enumerate(pdb_list):
            pdb_diameter = pdb.determineMaxDiameter()
            if pdb_diameter > max_diameter:
                max_diameter = pdb_diameter
                max_diameter_id = i
            apbs.runPDB2PQR(pdb.PDBFilename,
                            pdb.PDBFilename[:-4] + '.pqr')
                
        # set padding to 1/2 of max diameter
        padding = int(ceil(max_diameter))
        
        template_in = InFile("", "", "", mesh_size)
        
        # create template InFile for the biggest structure
        template_in.generateFromPDB(pdb_list[max_diameter_id], padding, True)
        
        # run apbs for all structures
        for pdb in pdb_list:
            template_in.setPQRFilePath(pdb.PDBFilename[:-4] + '.pqr')
            template_in.setOutSurfaceDXPath(pdb.PDBFilename[:-4] + "_vdw")
            template_in.setOutPotentialDXPath(pdb.PDBFilename[:-4] + "_esp")
            apbs.runAPBS(template_in, pdb.PDBFilename[:-4] + '.in')
        
        # clean up
        if clean_up == True:
            call(['rm', '*.pqr'])
            call(['rm', '*.in'])
        
        if hull_coordinates_exist == False:
            # load vdw_grid of the biggest hull
            big_box = DXReader().parse(
                pdb_list[max_diameter_id].PDBFilename[:-4] + '_vdw-PE0.dx',
                DXReader().VDW, mesh_size).getBox()
            template_box = big_box.box.copy()
            
            """
            Iterate over all pdb's and multiply the template_box with each
            hull. Then template_box is 0 at every point where a protein is and
            the hull of the biggest structure can be extended until all 0
            are inside the hull.
            """
            for i, pdb in enumerate(pdb_list):
                if i != max_diameter_id:
                    current_box = DXReader().parse(
                        pdb.PDBFilename[:-4] + '_vdw-PE0.dx',
                        DXReader().VDW, mesh_size).getBox()
                    template_box = template_box * current_box.box
            
            #
            # debug
            #
            DXWriter().writeBox('testbox.dx', DXBox(template_box,
                                                    big_box.getDimensions(),
                                                    big_box.getMeshSize(),
                                                    big_box.getOffset()))    
            
            if extend_mode == 'all':
                self.extend_all(big_box, template_box)
            elif extend_mode == 'biggest':
                self.extend_biggest(big_box, template_box, pdb_list,
                                    max_diameter_id)
            else: 
                print("Unknown extension mode {0}! Existing!".format(extend_mode))
                sys.exit(1)
            
            #
            # debug
            #
            DXWriter().writeBox('hull_box.dx', big_box)
            
            # extract the hull
            big_box.find_surface()
            hull_coordinates = nonzero(big_box.box == big_box.score_of_surface)
            ESP_Profile_Manager().write_hull_coordinates(hull_coordinates)
        else:
            hull_coordinates = ESP_Profile_Manager().read_hull_coordinates(hull_coordinates_file)
        
        # check if out_file_name exists, if so, then set append to False
        if os.path.exists(out_file_name):
            result_append = False
        else:
            result_append = True
        
        for pdb in pdb_list:
            dx_box = DXReader().parse(pdb.PDBFilename[:-4] + '_esp-PE0.dx',
                                      DXReader().ESP, mesh_size).getBox()
            ESP_Profile_Manager().write_profile(dx_box.box, hull_coordinates,
                pdb.PDBFilename[:-4], out_file_name, result_append)
            result_append = True
        
        print("done after {0}s".format(time.time() - starttime))
        
    def read_pdb_list(self, pdb_list_name):
        """
        Reads pdb-filenames from a given list and loads them as pdb structures.
        """
        pdb_list = []
        infile = open(pdb_list_name, 'r')
        for line in infile:
            line = line.rstrip()
            pdb_list.append(PDBFile(line))
        infile.close()
        return pdb_list 
    
    def extend_all(self, big_box, template_box, extend_by = 1):
        """
        This method works like this:
        template_box is 0 on all points, where a protein is and so all
        points are already inside the 'hull'. Then the hull is extended by a
        given unit and thats it.
        The hull sits much more closely around all structures.
        """
        big_box.box = template_box
        
        # flood big_box
        big_box.flood()
        
        # extend by
        big_box.extendSurface(extend_by)
        
    def extend_biggest(self, big_box, template_box, pdb_list, max_diameter_id):
        """
        Flip template_box 1->0, 0->1, then multiplication with the 
        biggest hull should give all 0, if not, then the hull needs to
        be extended.                
        """
        template_box = (template_box - 1) * (-1)
        
        # flood big_box
        big_box.flood()
        
        extend_counter = 0
        hull_found = False
        while hull_found is False:
            if 1 not in (template_box * big_box.box):
                hull_found = True
                break
            big_box.extendSurface(1)
            extend_counter = extend_counter + 1
            
        # extend surface 1 more, so that it is 1 point away from the border
        big_box.extendSurface(1)
        extend_counter = extend_counter + 1
        print("Extended largest pdb {0} by {1} units.".format(
            pdb_list[max_diameter_id].PDBFilename, extend_counter))
    
        
if __name__ == '__main__':
    CalculateESPProfiles(sys.argv[1:])
