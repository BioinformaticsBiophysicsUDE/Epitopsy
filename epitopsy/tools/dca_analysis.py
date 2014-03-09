"""
.. moduleauthor:: Ludwig Ohl <Ludwig.Ohl@uni-due.de>
.. module:: dca_analysis
    :platform: Unix, (Windows)
    :synopsis: A module for analysing dca data.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as color
import math


class DcaAnalysis(object):

    '''
    This class performs multiple operations,
    visualisations and analysis tools on the dca data.
    '''

    def __init__(self, prot_id, filename, alignment_file=None):
        '''Initialises the class by loading all information from the
        specified file.

        :param filename: path and name of the file.
        :type filename: str.
        :returns: None.

        '''
        self.info, self.mi_data, self.di_data = self.__read_in(
            prot_id + '/dca/' + filename)
        self.info['name'] = prot_id
        self.bindings = None
        if alignment_file:
            self.di_data_reduce = None
            self.__get_aligned_seq(prot_id + '/alignment/' + alignment_file)
            self.SeqMatReduce()

    def __read_in(self, filename):
        '''Reads in the parameters (N, M, Meff, q, pCount, x) and the data

        :param filename: path and name of the file.
        :type  filename: str.
        :returns: tuple -- (info, Mutual information, Direct information)

        '''
        with open(filename) as filedata:
            lines = filedata.readlines()
        begin = lines[0].find('N=')  # parameters begin here
        datainfo = lines[0][begin:].split()
        info = {'seqlen': int(datainfo[0].split('=')[1]),
                'M': int(datainfo[1].split('=')[1]),
                'Meff': round(float(datainfo[2].split('=')[1])),
                'q': int(datainfo[3].split('=')[1]),
                'pCount': float(datainfo[4].split('=')[1]),
                'x': float(datainfo[5].split('=')[1]),
                'method': lines[0][5:begin - 1]
                }
        # Create square matrices with seqlen dimension
        mi_data = np.zeros([info['seqlen'], info['seqlen']])
        di_data = np.zeros([info['seqlen'], info['seqlen']])
        for line in lines[2:]:
            inf = line.split()
            mi_data[int(inf[0]) - 1, int(inf[1]) - 1] = float(inf[2])
            di_data[int(inf[0]) - 1, int(inf[1]) - 1] = float(inf[3])
        return info, mi_data, di_data

    def heatmap(self, name=None):
        '''Generates a heat map of the direct information data.

        :param name: (optional) Name is displayed as title of plot
        :type  name: str

        Returns:
            None.
        '''
        fig = plt.figure(num=None, figsize=(10, 8), dpi=120,
                         facecolor='w', edgecolor='k')
        fig.add_subplot(111)
        # Encoding of the values to colors in the heat map
        cdict = {
            'red': ((0.0, 1.0, 1.0), (0.2, 0.0, 0.0), (1.0, 0.5, 0.5)),
            'green': ((0.0, 1.0, 1.0), (0.2, 0.5, 0.5), (1.0, 0.0, 0.0)),
            'blue': ((0.0, 1.0, 1.0), (0.2, 1.0, 1.0), (1.0, 0.0, 0.0))
        }
        # colorbar object
        colmap = color.LinearSegmentedColormap('my_colormap', cdict, 1024)
        if name:
            plt.title('DCA pair map of ' + name, size=20)
        else:
            plt.title('DCA pair map', size=20)
        plt.imshow(self.di_data, cmap=colmap)
        # Parameters are placed inside the plot
        plt.annotate("N = " +
                     str(self.info['seqlen']), xy=(self.info['seqlen'] /
                                                   4, 3 *
                                                   self.info['seqlen'] /
                                                   6), xytext=(self.info['seqlen'] /
                                                               4, 3 *
                                                               self.info['seqlen'] /
                                                               6), size='20', color='r')
        plt.annotate("M = " +
                     str(self.info['M']), xy=(self.info['seqlen'] /
                                              4, 4 *
                                              self.info['seqlen'] /
                                              6), xytext=(self.info['seqlen'] /
                                                          4, 4 *
                                                          self.info['seqlen'] /
                                                          6), size='20', color='r')
        plt.annotate("Meff = " +
                     str(self.info['Meff']), xy=(self.info['seqlen'] /
                                                 4, 5 *
                                                 self.info['seqlen'] /
                                                 6), xytext=(self.info['seqlen'] /
                                                             4, 5 *
                                                             self.info['seqlen'] /
                                                             6), size='20', color='r')
        plt.colorbar()

    def freq_hist(self, method='di', cutoff=None):
        '''Returns a histogram that shows the
        frequency of the values of the Direct information.

        :param method: (optional) Name of method ('di' or 'mi').
        :type  method: str.
        :returns: None.

        '''
        if method == 'di':
            data = self.di_data.reshape(-1)
        elif method == 'mi':
            data = self.mi_data.reshape(-1)
        else:
            raise Exception('No valid method specified.')

        if cutoff:
            hist, bin_edges = np.histogram(data[data > cutoff], bins=100)
        else:
            hist, bin_edges = np.histogram(data[np.nonzero(data)], bins=100)
        width = 0.7 * (bin_edges[1] - bin_edges[0])
        center = (bin_edges[:-1] + bin_edges[1:]) / 2
        plt.bar(center, hist, align='center', width=width)

    def binding_candidates(self, method='di'):
        '''The best candidates of binding sites are displayed in
        order from best to worse. The lower bound is the cutoff, if used.

        :param method: (optional) Name of method ('di' or 'mi').
        :type  method: str.

        :param cutoff: lower bound for best candidates.
        :type cutoff: float.

        returns: list of in order of strongest candidates.
        '''
        # if self.bindings.size():
        #    return self.bindings
        if method == 'di':
            data = self.di_data
        elif method == 'mi':
            data = self.mi_data
        # sorts the bindings according to strength,
        # beginning with thw weakest.
        best_di = data.reshape(-1).argsort()
        self.bindings = np.array(
            [(int(math.floor(i / data.shape[0])), i % data.shape[0]) for i in best_di[::-1]])

    def show_bindings(self, num=20):
        '''Plots the binding


        '''
        fig = plt.figure(num=None, figsize=(12, 10), dpi=120,
                         facecolor='w', edgecolor='k')
        fig.add_subplot(111)
        plt.scatter(self.bindings[0:num].T[0], self.bindings[0:num].T[1],
                    c=[self.di_data[x[0], x[1]] for x in self.bindings[0:num]], s=15)
        plt.colorbar()

    def SeqMatReduce(self):
        '''This lovely function uses the gaps of the
        aligned sequence and filteres out the gapped
        positions in the di matrix.

        input:
            matrix: dca matrix
            alignedSeq: aligned sequence of the protein
         output:
             filtered matrix
        '''
        sequence = np.array([alpha for alpha in self.info['aligned_seq']])
        # numpy magic: creates a bool array for filtering
        sieve = (sequence != '-')
        self.di_data_reduce = self.di_data[sieve][:, sieve]

    def __get_aligned_seq(self, alignment_file):
        # read alignment file
        with open(alignment_file) as data:
            alignments = data.read()
        self.info['aligned_seq'] = alignments.split(
            '>target sequence\n')[1].split('>')[0].replace('\n', '')
