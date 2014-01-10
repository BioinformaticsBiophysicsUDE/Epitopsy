import numpy as np
import matplotlib.pyplot as plt

class DcaAnalysis(object) :
    '''
    This class performs multiple operations, 
    visualisations and analysis tools on the dca data.
    '''
    def __init__(self,filename):
        self.info, self.miData, self.diData = self.__read_in(filename)

    def __read_in(self, filename):
        '''
        Reads in the parameters (N, M, Meff, q, pCount, x) and the dca data.
        
        Args:
            filename -> path and name of the file
        
        Returns:
            A tuple: info, Mutual information, Direct information
        '''
        with open(filename) as file:
            lines = file.readlines()
        DataInfo = lines[0].split()
        info = {'seqlen': int(DataInfo[3].split('=')[1]),
	    			'M': int(DataInfo[4].split('=')[1]),
	    			'Meff': round(float(DataInfo[5].split('=')[1])),
	    			'q': int(DataInfo[6].split('=')[1]),
	    			'pCount': float(DataInfo[7].split('=')[1]),
	    			'x': float(DataInfo[8].split('=')[1])
                }
        miData = np.zeros([info['seqlen'],info['seqlen']])
        diData = np.zeros([info['seqlen'],info['seqlen']])
        for line in lines[2:]:
            a = line.split()
            miData[int(a[0])-1,int(a[1])-1] = float(a[2])
            diData[int(a[0])-1,int(a[1])-1] = float(a[3])
        return info, miData, diData
    
    def pltHeatMap(self, name=None):
        '''
        Generates a heat map of the direct information data.
        
        Args:
            name -> (optional) Name is displayed as title of plot
            
        Returns:
            None.
        '''
        fig = plt.figure(num=None, figsize=(10, 8), dpi=120, facecolor='w', edgecolor='k')
        fig.add_subplot(111)
        # Encoding of the values to colors in the heat map
        cdict = {
          'red'  :  ( (0.0, 1.0, 1.0), (0.2, 0.0, 0.0), (1.0, 0.5, .5)),
          'green':  ( (0.0, 1.0, 1.0), (0.2, 0.5, 0.5), (1.0, 0.0, 0.0)),
          'blue' :  ( (0.0, 1.0, 1.0), (0.2, 1.0, 1.0), (1.0, 0.0, 0.0))
        }
        # colorbar object
        cm = matplotlib.colors.LinearSegmentedColormap('my_colormap', cdict, 1024)
        if name:
            plt.title('DCA pair map of '+name, size=20)
        else:
            plt.title('DCA pair map', size=20)
        heatmap = imshow(self.diData, cmap=cm)
        # Parameters are placed inside the plot
        plt.annotate("N = "+str(self.info['seqlen']), xy=(self.info['seqlen']/4, 3*self.info['seqlen']/6), xytext=(self.info['seqlen']/4, 3*self.info['seqlen']/6), size='20', color='r')
        plt.annotate("M = "+str(self.info['M']), xy=(self.info['seqlen']/4, 4*self.info['seqlen']/6), xytext=(self.info['seqlen']/4, 4*self.info['seqlen']/6), size='20', color='r')
        plt.annotate("Meff = "+str(self.info['Meff']), xy=(self.info['seqlen']/4, 5*self.info['seqlen']/6), xytext=(self.info['seqlen']/4, 5*self.info['seqlen']/6), size='20', color='r')
        plt.colorbar()
    
    def freqHist(self, method='di'):
        '''Returns a histogram that shows the 
        frequency of the values of the Direct information.
        '''
        if method == 'di':
            a = self.diData.reshape(-1)
        elif method == 'mi':
            a = self.miData.reshape(-1)
        else:
            raise Exception('No valid method specified.')
        hist, bin_edges = np.histogram(a[nonzero(a)], bins=100)
        width = 0.7 * (bin_edges[1] - bin_edges[0])
        center = (bin_edges[:-1] + bin_edges[1:]) / 2
        plt.bar(center, hist, align='center', width=width)