'''
Created on Jun 29,  2011

@author: chris
'''

import sys
import numpy as np

from epitopsy.cython import FFTCorrelationScoringCython
from epitopsy.DXFile import DXBox

try:
    import anfft
    __anfft__ = True
except:
    __anfft__ = False


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
        elif isinstance(box, np.ndarray):
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
        if isinstance(x, np.ndarray):
            z = x[2]
            y = x[1]
            x = x[0]
            
        trans = np.zeros(3, dtype = int)
        
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


class FFT_correlation_scoring(Scoring):
    '''
    This class has been once written by Niko but these functions are outdated
    and newer functions will replace the old ones step by step.
    The new functions have a name like this 'this_is_a_new_function', whereas
    old functions have a name like this 'thisWouldBeAnOldFunction'.

    It is important to set the protein interior to 1 and the solvent to 0, if
    one is interested in the geometric complementarity.
    '''
    def __init__(self, box_array, fft_module = None):
        # select FFT module
        if fft_module:
            if fft_module == 'anfft':
                if not __anfft__:
                    raise ImportError('Module anfft unavailable')
                else:
                    self.fft_module = 'anfft'
            elif fft_module in ('numpy', 'np', 'npfft'):
                self.fft_module = 'npfft'
            else:
                raise ImportError('Module {} not supported'.format(fft_module))
        else:
            if __anfft__:
                self.fft_module = 'anfft'
            else:
                self.fft_module = 'npfft'
        
        # initialize box
        if isinstance(box_array, np.ndarray):
            Scoring.__init__(self, box_array)
        elif isinstance(box_array, DXBox):
            Scoring.__init__(self, box_array)
        else:
            print("Error in FFT_correlation_scoring, this function requires "
                  "an array or a DXBox!!!")
            sys.exit(1)
        
        # compute FFT
        self.box = self.do_fft(self.box)
        self.box = self.box.conjugate()
        

    def do_fft(self, signal_matrix):
        """
        Maybe zero-padding is necessary, it could be done with:
        fftn(signal_matrix, [next pow of 2, next pow of 2, next pow of 2])
        """
        if self.fft_module == 'anfft':
            return anfft.fftn(signal_matrix)
        else:
            return np.fft.fftn(signal_matrix)

    def do_ifft(self, signal_matrix):
        if self.fft_module == 'anfft':
            return anfft.ifftn(signal_matrix)
        else:
            return np.fft.ifftn(signal_matrix)

    def do_fftshift(self, signal_matrix):
        """
        Shift the zero-frequency component to the center of the spectrum.
        Fallback method, usually :meth:`shift_fft` does the job faster.
        """
        return np.fft.fftshift(signal_matrix)

    def do_ifftshift(self, signal_matrix):
        """
        Shift the zero-frequency component to the left of the spectrum.
        Fallback method, usually :meth:`shift_fft` does the job faster.
        """
        return np.fft.ifftshift(signal_matrix)

    def shift_fft(self, signal_matrix, center = None):
        """
        Shift the zero-frequency component to the center of the spectrum,
        or to the left of the spectrum if applied twice.
        Faster than np.fft.fftshift() and np.fft.ifftshift().
        """
        if center is None:
            x_center = int(round(signal_matrix.shape[0] / 2.))
            y_center = int(round(signal_matrix.shape[1] / 2.))
            z_center = int(round(signal_matrix.shape[2] / 2.))
        else:
            x_center = center[0]
            y_center = center[1]
            z_center = center[2]

        return FFTCorrelationScoringCython.fft_shift(signal_matrix, x_center,
                                                     y_center, z_center)

    def score(self, signal_matrix, position = None, maxRadius = None):
        """
        Computes the FFT-form of the comparison box and scores it against the
        reference box, which is in a conjugated-complex FFT-form.
        """
        if isinstance(signal_matrix, DXBox):
            signal_clone = signal_matrix.box.copy()
        else:
            signal_clone = signal_matrix.copy()

        if signal_clone.shape != self.box.shape:
            raise NameError("ERROR in FFT_correlation_scoring.score: Grid dime"
                            "nsions differ. Cannot complete surface scoring!")

        #=======================================================================
        # Fourier Correlation Scoring
        #=======================================================================
        signal_clone = self.do_fft(signal_clone)
        signal_clone = signal_clone * self.box
        signal_clone = self.do_ifft(signal_clone)
        if position is None or maxRadius is None :
            return self.determine_best_correlation(self.do_ifftshift(signal_clone).real)
        else:
            return self.determineBestCorrelationWithinRadius(self.do_ifftshift(signal_clone).real, position, maxRadius)

    def determine_best_correlation(self, signal_matrix):
        best_correlation_score = np.amax(signal_matrix)
        best_position = np.nonzero(signal_matrix == best_correlation_score)
        correlation_vector = np.zeros(3)
        correlation_vector[0] = best_position[0][0]
        correlation_vector[1] = best_position[1][0]
        correlation_vector[2] = best_position[2][0]
        # calculate the translation vector
        correlation_vector = self.calculate_fft_translation(correlation_vector)

        result = Result()
        #print('***best score: {0}'.format(best_correlation_score))
        result.setScore(best_correlation_score)
        '''
        Vector points to optimal translation of fixed structure relative to
        rotated. Thus the optimal superposition of the rotated structure to the
        fixed is the inverted vector (-v)
        '''
        result.setGridTranslationX(-correlation_vector[0])
        result.setGridTranslationY(-correlation_vector[1])
        result.setGridTranslationZ(-correlation_vector[2])
        return result

    def determineBestCorrelationWithinRadius(self, signal_matrix, position, maxRadius):
        """
        This will be very slow ... I think the intention of this function is
        to find good correlations near the current position, so that it would
        for example be possible to search around a binding pocket and not get
        distracted by another (maybe stronger) binding pocket.
        """
        bestCorrelationScore = signal_matrix.real[0][0][0]
        correlationVector = np.zeros(3)
        # width of array
        w = signal_matrix.shape[0]
        # h height of array
        h = signal_matrix.shape[1]
        # depth of array
        d = signal_matrix.shape[2]
        for x in range(w):
            for y in range(h):
                for z in range(d):
                    currentCorrelationScore = signal_matrix.real[x][y][z]
                    if currentCorrelationScore > bestCorrelationScore:
                        v = np.array(x, y, z) - position
                        v = np.sqrt(np.dot(v, v))
                        if v < maxRadius:
                            bestCorrelationScore = currentCorrelationScore
                            correlationVector = [x, y, z]

        result = Result()
        result.setScore(bestCorrelationScore)
        correlationVector = self.calculate_fft_translation(correlationVector)
        '''
        Vector points to optimal translation of fixed structure relative ti
        rotated. Thus the optimal superposition of the rotated structure to the
        fixed is the inverted vector (-v)
        '''
        result.setGridTranslationX(-correlationVector[0])
        result.setGridTranslationY(-correlationVector[1])
        result.setGridTranslationZ(-correlationVector[2])
        return result
    def clone(self):
        return FFT_correlation_scoring(self.box.copy())

    def get_correlation(self, signal_matrix):
        """
        This method calculates the correlation with fft's.
        In this notation matrixA is fixed, whereas matrixB is the rotated one.
        The method returns only the real values, as they are the interesting
        part and it also shifts the zero frequencies to the middle.
        """
        signal_copy = self.do_fft(signal_matrix)
        signal_copy = self.box * signal_copy
        signal_copy = self.do_ifft(signal_copy)
        return signal_copy.real

    def get_imag_correlation(self, signal_matrix):
        """
        This method calculates the correlation with fft's.
        In this notation matrixA is fixed, whereas matrixB is the rotated one.
        The method returns only the real values, as they are the interesting
        part and it also shifts the zero frequencies to the middle.
        """
        signal_copy = self.do_fft(signal_matrix)
        signal_copy = self.box * signal_copy
        signal_copy = self.do_ifft(signal_copy)
        return signal_copy


