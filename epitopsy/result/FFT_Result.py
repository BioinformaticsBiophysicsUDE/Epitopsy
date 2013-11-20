'''
Created on 14.11.2011

@author: chris
'''

import os
import sys
import numpy as np
import operator

from epitopsy.Structure import PDBFile

class FFT_Result(object):
    '''
    This class stores the results from FFT Correlation Docking.
    It stores them in a list, which contains a dictionary for every entry.
    
    The format looks like this:
    
    x_coord y_coord z_coord phi theta psi score.
    '''


    def __init__(self):
        '''
        some stuff
        '''
        self.results = []
    
    def _add_score(self, x, y, z, phi, theta, psi, score):
        '''
        This method adds a score to the list.
        '''
        new_element = {}
        new_element['x'] = x
        new_element['y'] = y
        new_element['z'] = z
        new_element['phi'] = phi
        new_element['theta'] = theta
        new_element['psi'] = psi
        new_element['score'] = score
        self.results.append(new_element)
        
    def _sort(self):
        '''
        Sort the results by their score, with the highest beeing the first.
        '''
        self.results.sort(key = operator.itemgetter('score'), reverse = True)
    
    def find_scores(self, box, phi, theta, psi, number_of_elements,
                    dxbox = None):
        '''
        This method finds the highest scoring elements from the box array. It
        needs the dxbox to transform the coordinates to real space, if none is
        given, it will store the box indices.
        '''
        # make a copy of the box
        mod_box = box.copy()
        
        for i in range(0, number_of_elements):
            # find score
            score, [i, j, k] = self._return_highest_score(mod_box)
            # transform indices [x, y, z] to real space
            if dxbox is not None:
                # transform the fft-shift
                [xdim, ydim, zdim] = dxbox.box_dim
                if i < xdim / 2:
                    x = i
                else:
                    x = i - xdim
                if j < ydim / 2:
                    y = j
                else:
                    y = j - ydim
                if k < zdim / 2:
                    z = k
                else:
                    z = k - zdim
                x = xdim / 2 - x - 1
                y = ydim / 2 - y - 1
                z = zdim / 2 - z - 1
                indices = [x, y, z]
                indices = dxbox.transform_box_to_real_space(indices) 
            
            # add score
            self._add_score(x = indices[0], y = indices[1], z = indices[2],
                            phi = phi, theta = theta, psi = psi, score = score)
            
        
        
    def _return_highest_score(self, mod_box):
        '''
        This method returns the highest score + the box coordinates of this 
        score. For further calls, it sets the highest score to the lowest, so 
        that the next highest score can be returned the next time.
        '''
        # find highest score
        highest_score = np.amax(mod_box)
        # find coordinates
        coord = np.nonzero(mod_box == highest_score)
        # save coordinates
        box_indices = [coord[0][0], coord[1][0], coord[2][0]]
        # set found score to the lowest score
        mod_box[box_indices[0], box_indices[1], box_indices[2]] = np.amin(mod_box)
        return highest_score, box_indices
    
    def write_to_disk(self, filename):
        '''
        Write results to disk. Before doing so, the list will be sorted by the
        score.
        '''
        # Make sure the results are sorted.
        self._sort()
        
        with open(filename, 'w') as f:
            for item in self.results:
                w_line = ('{0} {1} {2} {3} {4} {5} {6}\n'
                          .format(item['x'], item['y'], item['z'],
                                  item['phi'], item['theta'], item['psi'],
                                  item['score']))
                f.write(w_line)
        
    def read_from_filename(self, filename):
        '''
        Read the results from disk.
        '''
        with open(filename, 'r') as f:
            for line in f:
                if line.rstrip() != '':
                    content = line.rstrip().split(' ')
                    x = float(content[0])
                    y = float(content[1])
                    z = float(content[2])
                    phi = float(content[3])
                    theta = float(content[4])
                    psi = float(content[5])
                    score = float(content[6])
                    self._add_score(x = x, y = y, z = z,
                                    phi = phi, theta = theta, psi = psi,
                                    score = score)
    
    def make_pdb_results(self, pdb_path, storage_dir, new_name, best_x):
        '''
        This method takes the 'best_x' results and moves the pdb structure 
        to these positions. It stores the moved and rotated structures in the
        specified directory. The name of the pdb works like this:
            'new_name' + '_number.pdb', with number indexing the results from 
            0 as the best to the end.
            
        This method assuemes, that the given coordinates are in real space!!!
        '''
        # check if there are enough results
        if best_x > len(self.results):
            print('Not enough results found! Writing all results ...')
            best_x = len(self.results)
            
        # sort the list
        self._sort()
        
        # check if storage dir exists
        if not os.path.exists(storage_dir):
            os.mkdir(storage_dir)
        
        for i in range(0, best_x):
            x = self.results[i]['x']
            y = self.results[i]['y']
            z = self.results[i]['z']
            phi = self.results[i]['phi']
            theta = self.results[i]['theta']
            psi = self.results[i]['psi']
            # load pdb, its easier to load it every time again, than to create
            # copys of the pdb structure
            pdb = PDBFile(pdb_path)
            pdb.translate_origin_and_rotate(phi, theta, psi)
            pdb.translate([x, y, z])
            # make filename
            new_filename = '{0}_{1}.pdb'.format(new_name, i)
            # update filename
            new_filename = os.path.join(storage_dir, new_filename)
            pdb.save_to_file(new_filename)
    
