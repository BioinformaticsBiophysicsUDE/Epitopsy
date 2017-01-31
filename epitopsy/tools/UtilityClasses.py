__author__     = "Christoph Wilms, Jean-Noel Grad"
__copyright__  = "Copyright 2011, Epitopsy"
__date__       = "Jun 29,  2011"
__credits__    = ["Christoph Wilms", "Jean-Noel Grad"]
__license__    = "LGPL v3"
__doc__        = """Various tools."""

import time
import sys
import datetime
import numpy as np


class Progress_bar:
    '''
    This class builds a nice progress bar:

    0%                                    100%
    #                                        #
    #********                                #
    #********************                    #


    There are 50 characters between the #.
    '''
    def __init__(self, number_of_item, refresh = 10):
        self.number_of_item = number_of_item
        self.finished = 0.0
        self.counter = 0
        self.start_time = time.time()
        self.progress_bar = ''
        self.refresh = 15
        print('0%                                             100%\n'
              '#                                                 #')
    def add(self):
        self.finished = self.finished + 1
        if self.finished > self.number_of_item:
            return
        elif self.finished == self.number_of_item:
            time_spent = round((time.time() - self.start_time) / 60., 2)
            print('#*****************    DONE    ********************#')
            print('Execution time: {0} min.'.format(time_spent))
            return
        elif int(self.finished / self.number_of_item * 100. / 2.) > self.counter * 50. / self.refresh:
            self.counter = self.counter + 1
            self.progress_bar = '*' * int(self.finished / self.number_of_item * 100. / 2.)
            self.progress_bar = self.progress_bar + ' ' * int(49. - len(self.progress_bar)) + '#'
            print('{0}{1}'.format('#', self.progress_bar))

class Progress_bar_countdown(object):
    """
    Display a progress bar, 70 characters wide, with countdown if requested.
    The countdown uses a linear combination of the 4 last stored timings to
    compute the remaining time.
    
    :param iterations: total number of iterations
    :type  iterations: int
    :param refresh: number of iteration to complete before refreshing the
        progress bar (optional), default is to update every 10 %, namely for
        **iterations** = 150, the progress bar is updated every 15 increments;
        can be given as a percentage (0.10) or as a number of increments (15)
    :type refresh: int or float
    :param countdown: display a countdown (optional), default is ``False``
    :type  countdown: bool
    
    .. attribute:: iterations
    
        total number of iterations
    
    .. attribute:: counter
    
        current state of the iteration
    
    .. attribute:: finished
    
        ``True`` if iteration has finished, prevents the progress bar from
        overflowing in the case where **iterations** was underestimated
    
    .. attribute:: checkpoints
    
        all checkpoints to pass to get an update of the progress bar
    
    .. attribute:: timer
    
        timing of each iteration
    
    
    The progress bar takes 70 characters with countdown,  50 otherwise.
    
    Example::
    
        >>> from epitopsy.progress_bar import Progress_bar
        >>> progress_bar = Progress_bar_countdown(len(array), refresh = 1, countdown = True)
        >>> for i in range(len(array)):
        ...     function(array[i])
        ...     progress_bar.iterate()
        0%                                            100%   time left
        #+++++++++++++++++++++++++++++++++++++++         #   0:01:36        

    """
    
    def __init__(self, iterations, refresh = None, countdown = False):
        self.iterations = iterations
        self.counter = 0
        self.finished = False
        self.start_time = time.time()
        self.timer = []
        
        # rangecheck: refresh
        if not refresh or refresh > iterations/4 or refresh < 0: # if no value or illegal value
            refresh = max(1, int(iterations/10.))                # set to 10% of **iterations** (minimum value of 1)
        elif type(refresh) == float:
            if 0.0 < refresh < 1.0:                              # if percentage, multiply by **iterations**
                refresh = int(refresh * iterations)
            else:                                                # if not, consider as an int
                refresh = int(refresh)
        elif type(refresh) == int:                               # if already int, use it as is
            pass
        
        # use refresh to create checkpoints
        self.checkpoints = np.array(range(refresh, iterations, refresh) + [iterations])
        
        # print first line
        if countdown:
            self.timer.append([0,0])
            print("0%                                            100%   time left")
        else:
            print("0%                                            100%")
    
    
    def increment(self, amount = 1):
        """
        Increment the counter of **amount** unit(s).
        
        :param amount: value of the increment (optional), default is unity
        :type  amount: int
        
        :returns: 0 if no error, 1 if :attr:`finished` was ``True``
        """
        
        if self.finished:
            return 1
        
        self.counter += amount
        
        # managing countdown
        if self.timer:
            self.timer.append([self.counter, time.time() - self.start_time])
        remains = ''
        
        if np.any(self.counter >= self.checkpoints):
            if self.counter >= self.iterations:
                self.terminate()
            else:
                if len(self.timer) >= 4:
                    # the countdown is a linear combination of the last 4
                    # durations with weights [1, 0.5, 0.33, 0.25]
                    weights = 1 / np.arange(1, 5, dtype = float)[::-1]
                    val = np.array(self.timer[-5:], dtype = float)
                    spent = (val[1:,1] - val[:-1,1]) / (val[1:,0] - val[:-1,0])
                    spent = np.sum(spent * weights) / np.sum(weights)
                    remains = int(spent * (self.iterations - self.counter))
                    remains = "   {0:<15}".format(datetime.timedelta(seconds = remains))
                bar = "\r#{0:<48}#".format(int(48 * self.counter / float(self.iterations)) * '+')
                sys.stdout.write(bar + remains)
                sys.stdout.flush()
                self.checkpoints = np.array([x for x in self.checkpoints if x > self.counter])
        
        return 0
        
    def add(self):
        """
        Alias of :meth:`Progress_bar.increment` for retro-compatibility.
        Add 1 to the counter.
        """
        
        self.increment(amount = 1)
    
    def terminate(self):
        """
        Terminate the timer and print "DONE" to screen.
        
        :returns: 0 if no error, 1 if :attr:`finished` was already ``True``
        """
        
        if self.finished:
            return 1
        
        self.finished = True
        
        if self.timer:
            print "\r#{0:+^48}#   {1:<15}".format("   DONE   ", datetime.timedelta(seconds = 0))
        else:
            print "\r#{0:+^48}#".format("   DONE   ")
        
        return 0

class Fix_xxmer:
    '''
    This class fixes xxmers, where xx copies of one protein are in one pdb.
    '''
    def __init__(self, pdb_filename):
        self.filename = pdb_filename
        self.fix_name = '_fixed'

    def hack_linenumbering(self):
        reader = open(self.filename, 'r')
        writer = open(self.filename[:-4] + self.fix_name + '.pdb', 'w')
        atom_counter = 'not set'
        residue_counter = 'not set'
        current_res = 0
        for line in reader:
            if atom_counter == 'not set' and residue_counter == 'not_set':
                atom_counter = int(line[6:11])
                residue_counter = int(line[22:26])
                current_res = residue_counter
                continue
            # not working yet

    def hack_chain_ids(self):
        '''
        Make double occuring chain ids unique!
        '''
        unique_id = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N',
                     'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', '0', '1',
                     '2', '3', '4', '5', '6', '7', '8', '9', ' ']
        reader = open(self.filename, 'r')
        writer = open(self.filename[:-4] + self.fix_name + '.pdb', 'w')
        chainID = None
        residue_counter = None
        counter = 0
        for line in reader:
            if line[0:4] == 'ATOM':
                if chainID is None and residue_counter is None:
                    chainID = line[21:22]
                    residue_counter = int(line[22:26])
                elif line[21:22] != chainID or int(line[22:26]) < residue_counter:
                    counter = counter + 1
                    chainID = line[21:22]
                if counter > len(unique_id):
                    print("ERROR: Too many chains, chain hack does not work!")
                    reader.close()
                    writer.close()
                    sys.exit(1)
                newline = line[0:21] + unique_id[counter] + line[22:]
                residue_counter = int(line[22:26])
                writer.write(newline)
        reader.close()
        writer.close()

class ESP_Profile_Manager:
    '''
    This class writes ESP-profiles to file.
    '''

    def __init__(self):
        pass

    def write_profile(self, box, indices, id, filename, append = False):
        '''
        Input:
        box: array with 3 dimensions
        indices: tupel of arrays, which you get for example with nonzero(a==1)
        id; identifier for one ESP-profile in case you want to append it to an
            existing file
        append: True or False

        The formatting looks like this (seperated by \t):
        X        [coord_1    coord_2_... coord_n]
        Y        [coord_1    coord_2_... coord_n]
        Z        [coord_1    coord_2_... coord_n]
        [ID_1]   [phi_1        phi_2 ...   phi_n]
        '''
        try:
            # write x line
            xline = 'X'
            for item in indices[0][:]:
                xline = "{0}\t{1}".format(xline, item)
            # write y line
            yline = 'Y'
            for item in indices[1][:]:
                yline = "{0}\t{1}".format(yline, item)
            # write z line
            zline = 'Z'
            for item in indices[2][:]:
                zline = "{0}\t{1}".format(zline, item)
            # write phi values, 1 because it is the fist entry
            philine = str(id)
            phi_value = box[indices]
            for item in phi_value:
                philine = "{0}\t{1:e}".format(philine, item)
        except:
            print("ERROR: Writing file {0} failed!".format(filename))
            sys.exit(1)
        # write
        if append == False:
            writer = open(filename, 'w')
        elif append == True:
            writer = open(filename, 'a')
        writer.write(xline + '\n' + yline + '\n' + zline + '\n' +
                     philine + '\n\n')
        writer.close()

    def read_hull_coordinates(self, hull_coordinates_file, phi_values = False, id = 1):
        """
        Reads hull coordinates from a file and returns them as a list, so
        they can directly access the values in the grid.

        There should be only one set of hull_coordinates per file!

        The method can read the corresponding phi values as well, then it only
        needs an id, if the id differs from 1.

        The expected format is like the one shown in write_hull_coordinates.
        """
        infile = open(hull_coordinates_file, 'r')

        # 3 lists for x, y and z
        xlist = []
        ylist = []
        zlist = []
        phi_list = []

        for line in infile:
            if line[0] == 'X':
                line_content = (line.rstrip()).split('\t')[1:]
                for item in line_content:
                    xlist.append(int(item))

            if line[0] == 'Y':
                line_content = (line.rstrip()).split('\t')[1:]
                for item in line_content:
                    ylist.append(int(item))

            if line[0] == 'Z':
                line_content = (line.rstrip()).split('\t')[1:]
                for item in line_content:
                    zlist.append(int(item))

            if phi_values == True and line[0] == str(id):
                line_content = (line.rstrip()).split('\t')[1:]
                for item in line_content:
                    phi_list.append(float(item))
        if phi_values == True:
            return ([np.array(xlist), np.array(ylist), np.array(zlist)],
                                                      np.array(phi_list))
        else:
            return [np.array(xlist), np.array(ylist), np.array(zlist)]

    def write_hull_coordinates(self, hull_coordinates, filename = 'hull_coordinates.txt'):
        """
        Write the hull coordinates like this (seperated by '\t'):

        X        [coord_1    coord_2_... coord_n]
        Y        [coord_1    coord_2_... coord_n]
        Z        [coord_1    coord_2_... coord_n]
        """
        #if os.path.exists(filename):
        #    print("File {0} already exists, hull coordinates will not be written!".format(
        #        filename))
        #    return
        # write x line
        xline = 'X'
        for item in hull_coordinates[0][:]:
            xline = xline + '\t' + str(item)
        # write y line
        yline = 'Y'
        for item in hull_coordinates[1][:]:
            yline = yline + '\t' + str(item)
        # write z line
        zline = 'Z'
        for item in hull_coordinates[2][:]:
            zline = zline + '\t' + str(item)
        infile = open(filename, 'w')
        infile.write(xline + '\n' + yline + '\n' + zline + '\n')

