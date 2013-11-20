'''
Created on Jul 4, 2011

@author: chris
'''

from Scoring import Scoring
from epitopsy.dx.ESPBox import ESPBox
from epitopsy.result.PotentialResult import PotentialResult

from numpy import abs

class PotentialScoring(Scoring):
    '''
    classdocs
    '''


    def __init__(self, refBox):
        '''
        DXBox refBox
        '''
        Scoring.__init__(self, refBox)
    def score(self, compBox):
        '''
        ESPBox compBox
        
        This method does a node per node comparison of two equally dimensioned
        grids.
        It returns the sum of absolute distances over all nodes.
        '''
        diff = 0.0
        counter = 0
        if not(compBox.hasSameSizeAs(self.getBoxSize())):
            raise AttributeError("Grid dimensions differ. Cannot complete potential scoring.")
        for x in range(compBox.getDimensionX()):
            for y in range(compBox.getDimensionY()):
                for z in range(compBox.getDimensionZ()):
                    if not( self.refBox.getElement(x, y, z) == ESPBox().ESP_WIPEOUT_VALUE or
                            compBox.getElement(x, y, z) == ESPBox().ESP_WIPEOUT_VALUE):
                        '''
                        get diff^2, or not, code for that was commented...
                        '''
                        localDiff = abs(self.refBox.getElement(x, y, z) - 
                                        compBox.getElement(x, y ,z))  
                        diff = diff - localDiff
                        counter = counter + 1
        result = PotentialResult()
        if counter > 0:
            result.setScore(diff / counter)
        else:
            result.setScore(ESPBox().ESP_WIPEOUT_VALUE)
        return result