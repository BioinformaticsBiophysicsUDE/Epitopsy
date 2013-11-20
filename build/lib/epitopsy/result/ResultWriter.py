'''
Created on Jun 30, 2011

@author: chris
'''

import sys
from ResultSet import ResultSet, Result

'''
Imports at the bottom ... circular imports ... what a mess!
'''

class ResultWriter:
    '''
    Again I copied Niko's class, but I had to change the function print as this 
    is a keyword in python.
    print(...) -> printResult(...)
    '''
    

    def __init__(self):
        '''
        Constructor
        '''
        pass
    '''
    Writes a given ResultSet to specified file, overwriting old results.
    '''
    def write(self, filename, rs, append = False):
        try:
            if append != True:
                self.createResultFile(filename)
            writer = open(filename, 'a')
            if isinstance(rs, ResultSet):
                for r in rs.results:
                    writer.write(self.printResult(r))
            elif isinstance(rs, Result):
                writer.write(self.printResult(rs))
            else:
                raise NameError('Given argument to ResultWriter.write() is neither a ResultSet nor a Result!')
            writer.close()
        except NameError as warning:
            print(warning)
            sys.exit(1)
        except Exception as warning:
            print(warning)
            print("ERROR in ResultWriter.write: Failed to write to {0}!".format(filename))
            sys.exit(1)
    '''
    Prints Result as formatted String to given PrintWriter Object.
    '''
    def printResult(self, r):
        # TODO: check if this formatting is correct!
        str2write = "{0}\t{1:.2f}\t{2:.2f}\t{3:.2f}\t{4:.4f}\t{5:.3f}\t{6:.3f}\t{7:.3f}\n".format(r.getId(),
                                                                                                r.getTheta() + r.getThetaOffset(),
                                                                                                r.getPhi() + r.getPhiOffset(),
                                                                                                r.getPsi() + r.getPsiOffset(),
                                                                                                r.getScore(),
                                                                                                r.getCoordinateTranslationX(),
                                                                                                r.getCoordinateTranslationY(),
                                                                                                r.getCoordinateTranslationZ())
        return str2write
    def append(self, filename, rs):
        self.write(filename, rs, True)
    def createResultFile(self, path, psiAxis = False):
        writer = open(path, 'w')
        if psiAxis != False:
            preString = 'Psi-Axis  =  {0}/{1}/{2}\n'.format(psiAxis[0], psiAxis[1], psiAxis[2])
            writer.write(preString)
        str2write = 'id\t' + 'theta\t' + 'phi\t' + 'psi\t' + 'score\t' + 'translation_x\t'
        str2write = str2write + 'translation_y\t' + 'translation_z\t\n'
        writer.write(str2write)
        writer.close()


