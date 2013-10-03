'''
Created on Jun 29, 2011

@author: chris
'''

from ResultSet import ResultSet
from SurfaceResult import SurfaceResult

import sys

class ResultReader:
    '''
    classdocs
    '''
    def __init__(self):
        '''
        Constructor
        '''
        pass
    def parse(self,filename,header):
        rs = ResultSet()
        try:
            nLine = 0
            reader = open(filename,'r')
            '''
            Read header, if requested.
            Don't know how that should work...
            '''
            if header == True:
                line = reader.next()
            for line in reader:
                nLine = nLine+1
                components = line.split('\t')
                if len(components) != 8:
                    raise NameError('Malformed result expression at line {0} in file {1}!'.format(nLine,filename))
                r = SurfaceResult()
                r.setId(int(components[0]))
                r.setTheta(float(components[1]))
                r.setPhi(float(components[2]));
                r.setPsi(float(components[3]));
                r.setScore(float(components[4]));
                r.setCoordinateTranslationX(float(components[5]));
                r.setCoordinateTranslationY(float(components[6]));
                r.setCoordinateTranslationZ(float(components[7]));
                rs.append(r);
            reader.close
        except NameError as warning:
            print(warning)
            sys.exit(1)
        except Exception as warning:
            print(warning)
            sys.exit(1)
        return rs