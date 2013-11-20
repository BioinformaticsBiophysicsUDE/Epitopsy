'''
Created on Jun 29,  2011

@author: Christoph Wilms
'''

import time
import sys

from numpy import ceil, mod, array

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
        print('0%                                             100%\n#                                                 #')
    def add(self):
        self.finished = self.finished + 1
        if self.finished > self.number_of_item:
            return
        elif self.finished == self.number_of_item:
            print('#*****************    DONE    ********************#')
            print('Calculations finished after {0} min.'.format(round((time.time() - self.start_time) / 60., 2)))
            return
        elif int(self.finished / self.number_of_item * 100. / 2.) > self.counter * 50. / self.refresh:
            self.counter = self.counter + 1
            self.progress_bar = '*' * int(self.finished / self.number_of_item * 100. / 2.)
            self.progress_bar = self.progress_bar + ' ' * int(49. - len(self.progress_bar)) + '#'
            print('{0}{1}'.format('#', self.progress_bar))

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
            return [array(xlist), array(ylist), array(zlist)], array(phi_list)
        else:
            return [array(xlist), array(ylist), array(zlist)]

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

