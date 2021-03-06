import os
import re
import numpy as np

from epitopsy.Structure import PDBFile, PQRFile
from epitopsy.DXFile import read_dxfile, DXBox, VDWBox
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.metrics import silhouette_samples


def energy2occupancy(energy_path, microstates_path, occupancy_path=None,
                     normalize=False):
    '''
    Convert an energy grid to an occupancy grid using the formula
    :math:`occ = e^{-E}`.
    
    :param energy_path: path to the input energy grid
    :type  energy_path: str
    :param microstates_path: path to the input microstates grid
    :type  microstates_path: str
    :param occupancy_path: path to the output occupancy grid
    :type  occupancy_path: str
    :param normalize: normalize the DXBox
    :type  normalize: bool
    
    Example::
    
        >>> from epitopsy import EnergyGrid
        >>> EnergyGrid.energy2occupancy('protein-centered_epi.dx',
        ...                             'protein-centered_mic.dx.gz',
        ...                             'protein-centered_occ.dx')
    
    '''
    if occupancy_path is None:
        if '_epi.dx' in energy_path:
            occupancy_path = energy_path.split('_epi.dx')[0] + '_occ.dx'
        else:
            raise ValueError('Please provide an output path for the occupancy.')
    epi = read_dxfile(energy_path, 'esp')
    mic = read_dxfile(microstates_path, 'vdw')
    occ = epi
    # zero out values inside the protein (avoid overflows with exp of large values)
    epi.box *= np.array(mic.box >= 1, dtype=int)
    # compute exponential
    occ.box = np.exp(-epi.box)
    # zero out again inside the protein
    occ.box *= np.array(mic.box >= 1, dtype=int)
    # count how many values overflowed
    NaN = np.sum(np.isnan(occ.box))
    if NaN:
        occ.box[np.nonzero(np.isnan(occ.box))] = 1e12
        print('Warning: {} grid points had values too large and were '
              'replaced by {:.2g}'.format(NaN, np.max(occ.box)))
        if normalize:
            print('Warning: cannot normalize')
    elif normalize:
        occ.box /= np.sum(occ.box)
    occ.write(occupancy_path)


def isosurface_ramp(dxbox, output=None, name=None, percents=None):
    '''
    Find the volume of space in which there is a X% chance of finding a ligand
    (Highest Density Region). Write levels to a PyMOL script.
    
    :param dxbox: DXBox
    :type  dxbox: :class:`DXBox.DXBox` or np.ndarray or str
    :param output: output file for the PyMOL script
    :type  output: str
    :param name: custom name for the occupancy map (optional), default is DXBox
       filename
    :type  name: str
    :param percents: percentage level(s) (optional), default is 5-70%
    :type  percents: tuple(float) or float
    
    Example::
    
        >>> from epitopsy import EnergyGrid as EG
        >>> EG.isosurface_ramp('protein_occ.dx', percents=[0.01, 0.05, 0.10])
    
    '''
    if percents is None:
        percents = np.arange(0.05, 0.75, 0.05)
    if not output or not name:
        if isinstance(dxbox, DXBox):
            filename = dxbox.filename
        elif isinstance(dxbox, str):
            filename = dxbox
        else:
            raise TypeError('Cannot guess output name for ramp.')
        if not output:
            output = filename.split('.dx')[0] + '.pml'
        if not name:
            name = os.path.basename(filename).split('.dx')[0]
    
    levels = occupancy_HDR(dxbox, percent=percents)
    with open(output, 'w') as f:
        for level, percent in zip(levels, percents):
            f.write('# {:.0f}% probability of presence\n'.format(100 * percent))
            f.write('isosurface {0}_iso_{1:.0f}p, {0}, {2:.4g}\n'.format(
                    name, 100 * percent, level))


def occupancy_HDR(dxbox, percent=0.25):
    '''
    Find the volume of space in which there is a X% chance of finding a ligand
    (Highest Density Region).
    
    The DXBox is transformed into a linear vector and sorted in descending
    order. The HDR isosurface threshold is found at position *i* in the sorted
    array when the sum of the first *i* elements in the array account for X%
    of the sum of all grid points.
    
    :param dxbox: DXBox
    :type  dxbox: :class:`DXBox.DXBox` or np.ndarray or str
    :param percent: percentage level(s) (optional), default is 25%
    :type  percent: tuple(float) or float
    :resturns: Isosurface threshold values
    :rettype: tuple(float) or float
    
    Example::
    
        >>> from epitopsy import EnergyGrid
        >>> EnergyGrid.occupancy_HDR('protein_occ.dx', 0.25)
        0.08388
    
    '''
    if isinstance(percent, float):
        percent = [percent]
    for p in percent:
        if p <= 0 or p >= 1:
            raise ValueError('Percentages must lie between 0 and 1')
    
    if isinstance(dxbox, DXBox):
        box = dxbox.box
    elif isinstance(dxbox, np.ndarray):
        box = dxbox
    elif isinstance(dxbox, str):
        box = read_dxfile(dxbox, 'esp').box
    else:
        raise TypeError('Cannot parse type {}'.format(type(dxbox)))
    
    linear_vector = np.sort(box.flatten())[::-1]
    total = np.sum(linear_vector)
    steps = 10 ** np.arange(0, np.log10(linear_vector.shape[0]), 1)
    steps = np.array(steps, dtype=int)[::-1]
    thresholds = []
    for p in percent:
        target = p * total
        start = 0
        for step in steps: # simulated annealing
            for i in range(start, linear_vector.shape[0], step):
                if np.sum(linear_vector[0:i]) > target:
                    start = i - step # break one step before target value
                    break
        thresholds.append(linear_vector[i])
    
    if len(thresholds) == 1:
        return thresholds[0]
    else:
        return tuple(thresholds)



def microstates2LEV(microstates_path, output_path=None):
    '''
    Convert a microstate grid into a Ligand Excluded Volume grid, by listing
    all grid points where no rotation was allowed. The LEV grid is 1's inside
    and 0's outside. Can be plotted as an isosurface in PyMOL to represent the
    Ligand Accessible Surface (LAS).
    
    :param microstates_path: path to the input microstates grid
    :type  microstates_path: str
    :param output_path: path to the output LEV grid
    :type  output_path: str
    
    Example::
    
        >>> from epitopsy import EnergyGrid
        >>> EnergyGrid.microstates2LEV('protein-centered_mic.dx.gz',
        ...                            'protein-centered_lev.dx')
    
    '''
    if output_path is None:
        if '_mic.dx' in microstates_path:
            output_path = microstates_path.split('_mic.dx')[0] + '_lev.dx'
        else:
            raise ValueError('Please provide an output path for the LEV.')
    dxb = read_dxfile(microstates_path, 'vdw')
    eps = np.sum(np.abs(dxb.box - dxb.box.astype(int)))
    eps /= np.prod(dxb.box.shape)
    if eps > 1e-5:
        raise ValueError('The DXBox contains floating point values.\nPlease '
                         'check the argument, it should be an integer DXBox.')
    dxb.box = np.array(dxb.box < 1e-10, dtype=float)
    out = VDWBox(dxb.box, dxb.box_mesh_size, dxb.box_offset)
    out.setCommentHeader([
        'Output of function Epitopsy.EnergyGrid.tools.microstates2LEV(',
        '  "{}"'.format(microstates_path),
        '  "{}")'.format(output_path),
        '1\'s inside the ligand excluded volume and 0\'s outside'])
    out.write(output_path)


def ligand_grids(ligand_path, box_mesh_size, box_dim=None, box_offset=None,
                 box_type=('vdw', 'crg')):
    '''
    Create grids for the ligand charges and ligand Van der Waals surface.
    Help visualize how the ligand is interpreted by Epitopsy.
    
    :param ligand_path: paths to the ligand PQR file
    :type  ligand_path: str
    :param box_mesh_size: mesh size
    :type  box_mesh_size: tuple(float,float,float)
    :param box_dim: box dimesions (optional), by default use the ligand size
    :type  box_dim: tuple(int,int,int)
    :param box_offset: box origin (optional), default is the ligand center
    :type  box_offset: tuple(float,float,float)
    :param box_type: grids to compute: charges ('crg') or radii ('vdw')
    :type  box_type: tuple(str) or str
    
    Example::
    
        >>> from epitopsy import EnergyGrid
        >>> EnergyGrid.ligand_grids('sulfate.pqr', 3*[0.8,])
    
    '''
    if isinstance(box_type, basestring):
        box_type = (box_type,)
    
    # read ligand
    ligand = PQRFile(ligand_path)
    if box_dim is None:
        box_dim = ligand.get_dxbox_dim(box_mesh_size, extend=5,
          APBS_grid_size=False)
    if box_offset is None:
        box_offset = ligand.get_dxbox_offset(box_mesh_size, box_dim,
          ligand.determine_center_of_extremes())
    
    # compute grids
    for btype in box_type:
        output_file = '{1}_{0}.dx'.format(btype,
           os.path.basename(ligand_path)[::-1].split('.', 1)[1][::-1])
        if btype == 'vdw':
            box = ligand.snap_vdw_to_box(box_mesh_size, box_dim, box_offset)
            vdwbox = DXBox(box=box, meshsize=box_mesh_size, offset=box_offset)
            vdwbox.setCommentHeader([
               'Output of function Epitopsy.EnergyGrid.tools.ligand_grids()',
               'Van der Waals grid of ligand {}'.format(ligand_path)])
            vdwbox.write(output_file)
        elif btype == 'crg':
            box = ligand.snap_crg_to_box(box_mesh_size, box_dim, box_offset)
            crgbox = DXBox(box=box, meshsize=box_mesh_size, offset=box_offset)
            crgbox.setCommentHeader([
               'Output of function Epitopsy.EnergyGrid.tools.ligand_grids()',
               'Partial charges grid of ligand {}'.format(ligand_path)])
            crgbox.write(output_file)
        else:
            raise ValueError('Unknown box type "{}"'.format(btype))


def remove_trailing_digits(input_path, output_path, decimals):
    '''
    Truncate values in OpenDX grids to fewer decimal places. The default is to
    store 6 digits after the decimal point. Such a high accuracy is sometimes
    unnecessary and creates a lot of entropy that prevents efficient
    data compression.
    
    Use this tool with extreme caution. It is strongly advised not to remove
    trailing digits in APBS potential grid if they might be re-used in the
    future. Also, don't remove trailing digits in Epitopsy energy grids when
    they might be used as input for multiconformational averaging.
    
    :param input_path: path to the DXBox to truncate
    :type  input_path: str
    :param output_path: path to the output DXBox
    :type  output_path: str
    :param decimals: number of digits after the decimal point, between 2 and 5
    :type  decimals: int
    
    Example::
    
        >>> from epitopsy import EnergyGrid
        >>> EnergyGrid.remove_trailing_digits('protein_epi.dx',
        ...                                   'protein_epi_truncated.dx', 3)
    
    '''
    if decimals < 2:
        raise ValueError('Cannot truncate to less than 2 decimal places')
    if decimals >= 6:
        raise ValueError('Please choose a value between 2 and 5')
    pattern = '([0-9]\.[0-9]{{{0}}})[0-9]{{{1}}}e'.format(decimals, 6-decimals)
    replacement = '\g<1>{0}e'.format((6-decimals) * '0')
    content = open(input_path, 'r').read()
    m = re.search('\nobject 3 class array type double rank 0 items \d+ data'
                  ' follows', content)
    if not m:
        raise ValueError('Cannot parse DXBox header')
    header = content[:m.span()[1]]
    lines = content[m.span()[1]:]
    lines = re.sub(pattern, replacement, lines)
    open(output_path, 'w').write(header + lines)


def table_minima(minima):
    '''
    Count occurences of minima.
    
    :param minima: 3D position of local minima (as integers)
    :type  minima: np.array
    :returns: Minima counts
    :returntype: list(tuple(tuple(int),int))
    '''
    table = {}
    for v in minima:
        vector = tuple(v)
        table[vector] = table.get(vector, 0) + 1
    
    table = [(v, freq) for v, freq in table.items()]
    table = sorted(table, key=lambda x: x[1], reverse=True)
    return table


def cluster_kmeans(minima, k_max=12, discard_size=1, fun=np.median, seed=None):
    '''
    Cluster local minima of an energy grid detected by simulated annealing.
    Use silhouette coefficients to select the best *k*.
    
    To reduce the impact of outliers, rarely sampled local minima are ignored.
    
    :param minima: 3D position of local minima
    :type  minima: np.array
    :param cluster_max: maximum value for *k*
    :type  cluster_max: int
    :param discard_size: keep only local minima sampled more frequently than
       **discard_size**
    :type  discard_size: int
    :param fun: function to recompute clusters center, defautl is the median
    :type  fun: function
    :param seed: random generator seed (positive integer)
    :type  seed: int
    :returns: Clusters center and size
    :rettype: tuple(np.array,np.array)
    '''
    if seed is None:
        seed = np.random.randint(500)
    
    # discard minima with low occurence
    table = table_minima(minima)
    table = [(v, freq) for v, freq in table if freq > discard_size]
    minima = np.array([v for v, freq in table for _ in range(freq)])
    k_max = min(k_max, len(table) - 1)
    if k_max <= 2:
        raise ValueError('There are only {} unique data points, '
                         'clustering is meaningless'.format(k_max))
    
    # find optimal number of clusters
    best_n = None
    best_labels = None
    best_sil = -1
    for n in range(2, k_max):
        kmeans = KMeans(n_clusters=n, random_state=seed)
        kmeans.fit(minima)
        s = silhouette_samples(minima, kmeans.labels_, metric='euclidean')
        sil = np.mean(s)
        if sil > best_sil:
            best_sil = sil
            best_n = n
            best_labels = kmeans.labels_
    
    # compute clusters median
    centers = []
    for n in range(best_n):
        centers.append([(fun(minima[best_labels == n, 0]),
                         fun(minima[best_labels == n, 1]),
                         fun(minima[best_labels == n, 2])),
                        np.sum(best_labels==n)])
    centers = sorted(centers, key=lambda x: x[1], reverse=True)
    centers, freq = zip(*centers)
    return np.array(centers), np.array(freq)


def pseudo_atoms(group_name, name_fmt='P{}', **kwargs):
    '''
    Print PyMOL syntax for pseudo-atoms.
    
    :param group_name: name for the pseudo-atom object
    :type  group_name: str
    :param name_fmt: name formatting string, unless names are already in
       \*\*kwargs
    :type  name_fmt: str
    :param \*\*kwargs: optional key/value pairs, values can be
       str/int/float/np.array
    :type  \*\*kwargs: dict
    '''
    args_common = []
    args_iterable = []
    for key,val in kwargs.items():
        if isinstance(val, list) or isinstance(val, tuple) or isinstance(val, np.ndarray):
            args_iterable.append(key)
        else:
            args_common.append(key)
    if 'name' not in kwargs:
        kwargs['name'] = [name_fmt.format(i) for i in range(len(kwargs['pos']))]
        args_iterable.append('name')
    if len(kwargs['name'][-1]) > 4:
        raise ValueError('Cannot use name "{}": more than 4 characters for an '
                         'atom will crash PyMOL'.format(kwargs['name'][-1]))
    if isinstance(kwargs['pos'], basestring):
        kwargs['pos'] = [kwargs['pos']]
    lines = []
    for i in range(len(kwargs['pos'])):
        listing = ['pseudoatom {}'.format(group_name)]
        for key in args_iterable:
            val = kwargs[key][i]
            if isinstance(val, list) or isinstance(val, tuple) or isinstance(val, np.ndarray):
                val = list(val)
            listing.append('{}={}'.format(key, val))
        for key in args_common:
            val = kwargs[key]
            if isinstance(val, list) or isinstance(val, tuple) or isinstance(val, np.ndarray):
                val = list(val)
            listing.append('{}={}'.format(key, val))
        lines.append(', '.join(listing))
    return lines


