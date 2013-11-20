'''
Created on Jun 29,  2011

@author: chris
'''

from numpy import ndarray, zeros

from epitopsy.DXFile import DXBox

class Scoring:
    '''
    This class is in the process of an update of all functions!
    '''
    def __init__(self, box = None):
        if isinstance(box, DXBox):
            self.box = box.box.copy()
            self.box_size_x = box.box.shape[0]
            self.box_size_y = box.box.shape[1]
            self.box_size_z = box.box.shape[2]
            self.dimensions = [self.box_size_x, self.box_size_y, self.box_size_z]
        elif isinstance(box, ndarray):
            self.box = box.copy()
            self.box_size_x = box.shape[0]
            self.box_size_y = box.shape[1]
            self.box_size_z = box.shape[2]
            self.dimensions = [self.box_size_x, self.box_size_y, self.box_size_z]
        else:
            # do nothing
            pass
    
    def calculate_fft_translation(self, x, y = None, z = None):
        """
        Calculates the FFT-suggested optimal translation.
        """
        if isinstance(x, ndarray):
            z = x[2]
            y = x[1]
            x = x[0]
            
        trans = zeros(3, dtype = int)
        
        if(x < self.box_size_x / 2):
            dx = int(x)
        else:
            dx = int(x - self.box_size_x)
        if(y < self.box_size_y / 2):
            dy = int(y)
        else:
            dy = int(y - self.box_size_y)
        if(z < self.box_size_z / 2):
            dz = int(z)
        else:
            dz = int(z - self.box_size_z)
        trans[0] = dx
        trans[1] = dy
        trans[2] = dz
        return trans
        
    def getBoxSize(self): 
        return [self.box_size_x, self.box_size_y, self.box_size_z]
    
    def clone(self): 
        return Scoring(self.box.copy())

