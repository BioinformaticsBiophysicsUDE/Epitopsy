# @date Created on 30.11.2011
# @author: chris

import numpy as np
import matplotlib.pyplot as plt


def calculate_rotation_matrix(phi, theta, psi):
    '''
    http://mathworld.wolfram.com/EulerAngles.html

    phi = 0 - 360 about the z-axis
    
    theta = 0 - 180 about the new x'-axis
    
    psi = 0 - 360 about the new z'-axis

    Same formalism as in biopython Bio.PDB.Superimposer
    -> angles can be derived from Superimposer.rotran!
    '''
    if 0 <= phi <= 360:
        phi = np.deg2rad(phi)
    else:
        # print('{0} should be in the range [0, 360] and is {1}'.format('phi', phi))
        phi = np.deg2rad(phi)
    if 0 <= theta <= 180:
        theta = np.deg2rad(theta)
    else:
        # print('{0} should be in the range [0, 180] and is {1}'.format('theta', theta))
        theta = np.deg2rad(theta)
    if 0 <= psi <= 360:
        psi = np.deg2rad(psi)
    else:
        # print('{0} should be in the range [0, 360] and is {1}'.format('psi', psi))
        psi = np.deg2rad(psi)

    euler_rotation = np.zeros((3, 3))
    euler_rotation[0][0] = np.cos(psi) * np.cos(phi) - np.cos(theta) * np.sin(phi) * np.sin(psi)
    euler_rotation[0][1] = np.cos(psi) * np.sin(phi) + np.cos(theta) * np.cos(phi) * np.sin(psi)
    euler_rotation[0][2] = np.sin(psi) * np.sin(theta)
    euler_rotation[1][0] = -np.sin(psi) * np.cos(phi) - np.cos(theta) * np.sin(phi) * np.cos(psi)
    euler_rotation[1][1] = -np.sin(psi) * np.sin(phi) + np.cos(theta) * np.cos(phi) * np.cos(psi)
    euler_rotation[1][2] = np.cos(psi) * np.sin(theta)
    euler_rotation[2][0] = np.sin(theta) * np.sin(phi)
    euler_rotation[2][1] = -np.sin(theta) * np.cos(phi)
    euler_rotation[2][2] = np.cos(theta)
    return euler_rotation


def get_euler_angles_from_rotation_matrix(euler_rotation):
    '''
    This method returns the angles from a rotation matrix, derived from
    euler angles. Because only angles in the range [0,360] (and [0,180])
    are used, combinations due to the 2 pi periodicity are possible, but
    are of no interest.

    To get the angles of a combined rotation by two sets of
    [phi_1, theta_1, psi_1] and [phi_2, theta_2, psi_2]
    one has to multiply the rotationmatrices with 'np.dot(...)'!
    '''
    theta = np.rad2deg(np.arccos(euler_rotation[2, 2]))
    psi = np.rad2deg(np.arctan2(euler_rotation[0, 2], euler_rotation[1, 2]))
    phi = np.rad2deg(np.arctan2(euler_rotation[2, 0], -euler_rotation[2, 1]))

    return [phi, theta, psi]


def rotation_matrix_to_axis_angle(R, epsilon=1e-12):
    cos_angle = (np.trace(R) - 1) / 2.
    if cos_angle > 1 - epsilon:
        # no rotation (identity matrix)
        return np.array([1,0,0]), 0
    elif cos_angle < -1 + epsilon:
        # 180 degrees rotation
        for x in (-1,0,1):
            for y in (-1,0,1):
                for z in (-1,0,1):
                    R2 = Euler_Rodrigues_rotation_matrix([x,y,z], np.pi)
                    if np.sum(np.abs(R2 - R)) < epsilon:
                        axis = np.array([x,y,z])
                        angle = 180
                        return axis, angle
    else:
        angle = np.arccos(cos_angle)
        x = R[2,1] - R[1,2]
        y = R[0,2] - R[2,0]
        z = R[1,0] - R[0,1]
        axis = np.array([x,y,z])
        axis /= np.linalg.norm(axis)
        return axis, angle
    
    raise RuntimeError('Cannot find the axis')


def get_euler_angles_for_equal_distributed_rotation(N, fun=None):
    '''
    Compute unit vectors uniformly distributed on a sphere.
    
    :param N: number of vectors
    :type  N: int
    :param fun: function use to compute the coordinates of a spherical lattice
       (optional), default is :func:`spherical_lattice_default`
    :type  fun: function
    :return: Rotation angles
    :rtype: list(list(float,float,float))
    '''
    if fun is None:
        fun = spherical_lattice_default
    
    lattice = fun(N) # uniform lattice on a sphere
    probe = [1, 0, 0] # reference vector
    euler_angles = []
    for point in lattice:
        euler_angles.append(get_euler_angles(fixed_point=probe,
                                             rotate_point=point))
    
    return euler_angles


def spherical_lattice_default(N, sequence_offset=0, complementary=False):
    '''
    Calculate Cartesian coordinates of dots uniformly distributed on a sphere
    using a Fibonacci lattice generative spiral.
    
    :param N: number of dots in the lattice
    :type  N: int
    :param complementary: use the complementary spiral if ``True``, eg. golden
       angle :math:`360^{\\circ}\\phi^{-1} \\simeq 222.5^{\\circ}` instead of
       the default :math:`360^{\\circ}(1 - \\phi^{-1}) \\simeq 137.5^{\\circ}`
    :type  complementary: bool
    :param sequence_offset: Fibonacci sequence offset (optional)
    :type  sequence_offset: int
    :return: Coordinates in Cartesian space
    :rtype: list(tuple(float,float,float))
    
    Implementation adapted from the following sources:
    
    * Patrick Boucher, `Points on a sphere
      <http://www.softimageblog.com/archives/115>`_, *Softimage Blog* **2006**
    * Fnord, reply to `Evenly distributing n points on a sphere
      <https://stackoverflow.com/a/26127012>`_, *Stack Overflow* **2014**
    * Richard Swinbank and James Purser, Fibonacci grids: A novel approach to
      global modelling, *Quarterly Journal of the Royal Meteorological Society*
      **2006**, *132*: 1769--1793.
      DOI: `10.1256/qj.05.227 <http://doi.org/10.1256/qj.05.227>`_
    * Alvaro Gonzalez, Measurement of Areas on a Sphere Using Fibonacci and
      Latitude--Longitude Lattices, *Mathematical Geosciences*, **2010**,
      *42*: 49. DOI: `10.1007/s11004-009-9257-x
      <https://doi.org/10.1007/s11004-009-9257-x>`_
    
    Example::
    
        >>> from epitopsy.tools import MathTools
        >>> coord = MathTools.spherical_lattice_default(150)
        >>> pdb_fmt = ('ATOM   {:>4}  C{:0>2} ARG A{:>4}    {:>8.3f}{:>8.3f}'
        ...            '{:>8.3f}  1.00 00.00           C \\n')
        >>> with open('sphere.pdb', 'w') as f:
        ...     for i,(x,y,z) in enumerate(coord):
        ...         f.write(pdb_fmt.format(i+1, i%99+1, i//99+1, x, y, z))
        ...     f.write('END')
    
    '''
    golden_ratio = (1 + np.sqrt(5)) / 2
    offset = 2. / N
    polar_offset = 0.5 # uniform dot distribution at the poles (Swinbank 2006)
    
    # clockwise or counter-clockwise spiral (Gonzalez 2010)
    if complementary:
        golden_angle = 2 * np.pi * (1 / golden_ratio) # 222.5 degrees
    else:
        golden_angle = 2 * np.pi * (1 - 1 / golden_ratio) # 137.5 degrees
    
    # numpy implementation of the original spiral algorithm (Boucher 2006)
    i = np.arange(N)
    y = offset * (i + polar_offset) - 1
    r = np.sqrt(1 - y**2)
    angle = golden_angle * (i + sequence_offset) # sequence offset (Fnord 2006)
    x = r * np.cos(angle)
    z = r * np.sin(angle)
    coordinates = list(zip(x, y, z))
    
    return coordinates


def spherical_lattice_Boucher(N):
    '''
    Calculate Cartesian coordinates of dots uniformly distributed on a sphere
    using a Fibonacci lattice generative spiral.
    
    :param N: number of dots in the lattice
    :type  N: int
    :return: Coordinates in Cartesian space
    :rtype: list(tuple(float,float,float))
    
    The algorithm was published in:
    
    * Patrick Boucher, `Points on a sphere
      <http://www.softimageblog.com/archives/115>`_, *Softimage Blog* **2006**
    
    '''
    points = []
    
    inc = np.pi * (3 - np.sqrt(5)) # golden_ratio
    off = 2. / N
    for k in range(0, N):
        y = k * off - 1 + (off / 2)
        r = np.sqrt(1 - y * y)
        phi = k * inc
        points.append((np.cos(phi) * r, y, np.sin(phi) * r))
    
    return points


def spherical_lattice_Gonzalez(N):
    '''
    Calculate Cartesian coordinates of dots uniformly distributed on a sphere
    using a Fibonacci lattice generative spiral.
    
    :param N: number of dots in the lattice, must be even
    :type  N: int
    :return: *N* or *N* + 1 coordinates in Cartesian space
    :rtype: list(tuple(float,float,float))
    
    The algorithm was published in:
    
    * Alvaro Gonzalez, Measurement of Areas on a Sphere Using Fibonacci and
      Latitude--Longitude Lattices, *Mathematical Geosciences*, **2010**,
      *42*: 49. DOI: `10.1007/s11004-009-9257-x
      <https://doi.org/10.1007/s11004-009-9257-x>`_
    
    In the original algorithm, parameter *N* would return :math:`2N + 1`
    points. Due to the way this function is called in :mod:`EnergyGrid`, for a
    uneven argument *N* the function returns *N* coordinates, but for even *N*
    the function returns :math:`N + 1` coordinates.
    
    '''
    if N % 2 == 1:
        N = N - 1
    N = N // 2
    phi = (1 + np.sqrt(5)) / 2 # golden_ratio
    i = np.arange(-N, N + 1)
    latitude = np.arcsin(i * 2. / (2 * N + 1))
    longitude = np.mod(i, phi) * 2 * np.pi / phi
    longitude = np.mod(longitude, 2 * np.pi)
    x = np.cos(latitude) * np.cos(longitude)
    y = np.cos(latitude) * np.sin(longitude)
    z = np.sin(latitude)
    coordinates = list(zip(x, y, z))
    
    return coordinates


def spherical_lattice_Saff_Kuijlaars(N):
    '''
    Calculate Cartesian coordinates of dots uniformly distributed on a sphere
    using the Saff-Kuijlaars method.
    
    :param N: number of dots in the lattice
    :type  N: int
    :return: Coordinates in Cartesian space
    :rtype: list(tuple(float,float,float))
    
    The algorithm was published in:
    
    * Edward Saff and Arno Kuijlaars, Distributing many points on a
      sphere, *The Mathematical Intelligencer* **1997**, *19(1)*: 5--11.
      DOI: `10.1007/BF03024331 <https://doi.org/10.1007/BF03024331>`_
    
    Formula:
    
    :math:`h_k = \\frac{2(k - 1)}{N - 1} - 1, 1 \\leq k \\leq N`
    
    :math:`\\phi_k = \\phi_{k-1} + \\frac{3.6}{\\sqrt{N}}
    \\frac{1}{\\sqrt{1 - h_k^2}}, 2 \\leq k \\leq N - 1, \\phi_1 = \\phi_N = 0`
    
    '''
    s = 3.6 / np.sqrt(N)
    phi = 0
    coordinates = [(0, -1, 0)]
    for k in range(1, N - 1):
        y = 2 * k / float(N - 1) - 1 # (k-1) becomes k in a zero-based array
        r = np.sqrt(1 - y**2)
        phi = phi + s / r
        coordinates.append((r * np.cos(phi), y, r * np.sin(phi)))
    coordinates.append((0, 1, 0))
    
    return coordinates


def get_euler_angles(fixed_point, rotate_point):
    '''
    This method returns a rotation matrix which rotates the 'rotate_point'
    onto the fixed_point.
    '''
    # avoid catastrophic cancellation when the third component is zero
    if np.abs(rotate_point[2]) < 1e-6:
        rotate_point = (rotate_point[0], rotate_point[1], 1e-6)
    # add [0,0,0] so that the calculation works
    coords = [rotate_point, [0, 0, 0]]
    ref_coords = [fixed_point, [0, 0, 0]]

    correlation_matrix = np.dot(np.transpose(coords), ref_coords)
    u, d, vt = np.linalg.svd(correlation_matrix)
    rot = np.transpose(np.dot(np.transpose(vt), np.transpose(u)))
    # check for self reflection
    if np.linalg.det(rot) < 0:
        vt[2] = -vt[2]
        rot = np.transpose(np.dot(np.transpose(vt), np.transpose(u)))

    angles = get_euler_angles_from_rotation_matrix(rot)

    return angles


def get_neighbor_angle_set(n_rotations, max_dist):
    '''This method returns euler angles for rotations, that rotate the object
    to a oritentation, where the distance displacement is less than max_dist.
    Notice, that the returned angle set does not have the specified number
    of rotations!
    '''
    # get equal distributed points at first
    equal_points = spherical_lattice(n_rotations)

    # this is a reference probe
    probe = np.array([1., 0., 0.])
    angle_set = []
    for point in equal_points:
        point = np.array(point)
        if np.linalg.norm(point - probe) <= max_dist:
            angle_set.append(get_euler_angles(fixed_point=probe,
                                              rotate_point=point))

    return angle_set


def PCA_base(coord):
    '''
    Compute the unit vectors of a basis that maximizes the spread of atomic
    coordinates on the x-axis, then the y-axis and the z-axis.
    
    The PCA of a collection of atomic coordinates returns a base with
    determinant :math:`\pm 1`. Molecules with stereocenters are orientable,
    therefore a base with determinant :math:`+1` preserves the molecule
    topology while a base with determinant :math:`-1` introduces a reflection.
    The reflection component is removed to avoid inverting all stereocenters.
    
    :param coord: Atomic coordinates, usually the subset of C-alpha atoms
    :type  coord: :class:`numpy.ndarray`\[:,3\]
    :returns: Set of basis vectors, without reflection
    :returntype: :class:`numpy.ndarray`\[3,3\]
    '''
    covariance = np.cov(coord.T)
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    eigenvectors = eigenvectors[:, np.argsort(eigenvalues)[::-1]]
    if np.linalg.det(eigenvectors) < 0:
        eigenvectors *= -1
    return eigenvectors


def PCA_projection(eigenvectors, coord):
    '''
    Project atomic positions in a new basis set.
    
    :param coord: Atomic coordinates
    :type  coord: :class:`numpy.ndarray`\[:,3\]
    :param eigenvectors: Set of basis vectors
    :type  eigenvectors: :class:`numpy.ndarray`\[3,3\]
    :returns: Updated coordinates
    :rettype: :class:`numpy.ndarray`\[:,3\]
    '''
    return np.dot(eigenvectors.T, coord.T).T


def Lattman_distance(R):
    '''
    Compute the Lattman distance between rotation angles.
    
    :param R: rotation matrices
    :type  R: list(:class:`numpy.ndarray`\[3,3\])
    :returns: Lattmann distance matrix
    :returntype: :class:`numpy.ndarray`\[:,:\]
    
    The algorithm was published in:
    
    * Eaton Lattman, Optimal sampling of the rotation function,
      *Acta Crystallographica Section B* **1972**, *B28*: 1065--1068. DOI:
      `10.1107/S0567740872003723 <https://doi.org/10.1107/S0567740872003723>`_
    '''
    err_msg = 'Value outside the arc cosine definition domain: {:.17f}'
    chi_matrix = np.zeros([len(R), len(R)])
    for i in range(len(R) - 1):
        for j in range(i + 1, len(R)):
            a = (np.trace(np.dot(R[i], R[j].T)) - 1) / 2.
            if a < -1:
                if np.isclose(a, -1, rtol=0, atol=1e-12):
                    a = -1
                else:
                    raise RuntimeError(err_msg.format(a))
            elif a > 1:
                if np.isclose(a, 1, rtol=0, atol=1e-12):
                    a = 1
                else:
                    raise RuntimeError(err_msg.format(a))
            chi_matrix[i,j] = chi_matrix[j,i] = np.arccos(a)
    
    return chi_matrix


class RotationProtocol(object):
    '''
    Generic base class for a rotation protocol. A protocol initializes
    one or more rotation functions and applies them sequentially to a set of
    coordinates. This class is not supposed to be instanciated. When deriving
    a subclass, simply change :meth:`__init__`.
    
    Examples::
    
        >>> from epitopsy.tools import MathTools
        >>> import matplotlib.pyplot as plt
        >>> from matplotlib.backends.backend_pdf import PdfPages
        >>> 
        >>> # get one example from each protocol
        >>> p1 = MathTools.RotationProtocolFibonacci(400)
        >>> p2 = MathTools.RotationProtocolFibonacciSpin(100, 4)
        >>> p3 = MathTools.RotationProtocolLattman(4)
        >>> p4 = MathTools.RotationProtocolEulerUniform(7)
        >>> 
        >>> # plot
        >>> plt.subplot(121)
        >>> plot = p4.plot(plt,      markersize=5, label=p4.desc_short)
        >>> plot = p3.plot_add(plot, markersize=5, label=p3.desc_short)
        >>> plot = p2.plot_add(plot, markersize=5, label=p2.desc_short)
        >>> plot = p1.plot_add(plot, markersize=5, label=p1.desc_short)
        >>> plot.title('Sampling quality')
        >>> plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        >>> plt.margins(0.05, 0.05)
        >>> plt.tight_layout()
        >>> plt.show()
        >>> 
        >>> # Fibonacci spiral: constant N.S
        >>> f1 = MathTools.RotationProtocolFibonacciSpin(600, 1)
        >>> f2 = MathTools.RotationProtocolFibonacciSpin(300, 2)
        >>> f3 = MathTools.RotationProtocolFibonacciSpin(200, 3)
        >>> f4 = MathTools.RotationProtocolFibonacciSpin(150, 4)
        >>> f5 = MathTools.RotationProtocolFibonacciSpin(120, 5)
        >>> f6 = MathTools.RotationProtocolFibonacciSpin(100, 6)
        >>> # Fibonacci spiral: constant N
        >>> g1 = MathTools.RotationProtocolFibonacciSpin(150, 1)
        >>> g2 = MathTools.RotationProtocolFibonacciSpin(150, 2)
        >>> g3 = MathTools.RotationProtocolFibonacciSpin(150, 3)
        >>> g4 = MathTools.RotationProtocolFibonacciSpin(150, 4)
        >>> g5 = MathTools.RotationProtocolFibonacciSpin(150, 5)
        >>> g6 = MathTools.RotationProtocolFibonacciSpin(150, 6)
        >>> # plot constant N.S
        >>> plt.subplot(231)
        >>> f6.plot(plt,     metric='min', linemarker='-', label=f6.desc_short)
        >>> f5.plot_add(plt, metric='min', linemarker='-', label=f5.desc_short)
        >>> f4.plot_add(plt, metric='min', linemarker='-', label=f4.desc_short)
        >>> f3.plot_add(plt, metric='min', linemarker='-', label=f3.desc_short)
        >>> f2.plot_add(plt, metric='min', linemarker='-', label=f2.desc_short)
        >>> f1.plot_add(plt, metric='min', linemarker='-', label=f1.desc_short)
        >>> plt.title('')
        >>> plt.ylabel('Distance to closest\\nneighbor (degrees)')
        >>> plt.margins(0.05, 0.05)
        >>> plt.ylim(ymin=0)
        >>> plt.subplot(232)
        >>> f6.plot(plt,     metric='max', linemarker='-', label=f6.desc_short)
        >>> f5.plot_add(plt, metric='max', linemarker='-', label=f5.desc_short)
        >>> f4.plot_add(plt, metric='max', linemarker='-', label=f4.desc_short)
        >>> f3.plot_add(plt, metric='max', linemarker='-', label=f3.desc_short)
        >>> f2.plot_add(plt, metric='max', linemarker='-', label=f2.desc_short)
        >>> f1.plot_add(plt, metric='max', linemarker='-', label=f1.desc_short)
        >>> plt.title('')
        >>> plt.ylabel('Distance to farthest\\nneighbor (degrees)')
        >>> plt.margins(0.05, 0.05)
        >>> plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        >>> # plot constant N
        >>> plt.subplot(234)
        >>> g6.plot(plt,     metric='min', linemarker='-', label=g6.desc_short)
        >>> g5.plot_add(plt, metric='min', linemarker='-', label=g5.desc_short)
        >>> g4.plot_add(plt, metric='min', linemarker='-', label=g4.desc_short)
        >>> g3.plot_add(plt, metric='min', linemarker='-', label=g3.desc_short)
        >>> g2.plot_add(plt, metric='min', linemarker='-', label=g2.desc_short)
        >>> g1.plot_add(plt, metric='min', linemarker='-', label=g1.desc_short)
        >>> plt.title('')
        >>> plt.ylabel('Distance to closest\\nneighbor (degrees)')
        >>> plt.margins(0.05, 0.05)
        >>> plt.ylim(ymin=0)
        >>> plt.subplot(235)
        >>> g6.plot(plt,     metric='max', linemarker='-', label=g6.desc_short)
        >>> g5.plot_add(plt, metric='max', linemarker='-', label=g5.desc_short)
        >>> g4.plot_add(plt, metric='max', linemarker='-', label=g4.desc_short)
        >>> g3.plot_add(plt, metric='max', linemarker='-', label=g3.desc_short)
        >>> g2.plot_add(plt, metric='max', linemarker='-', label=g2.desc_short)
        >>> g1.plot_add(plt, metric='max', linemarker='-', label=g1.desc_short)
        >>> plt.title('')
        >>> plt.ylabel('Distance to farthest\\nneighbor (degrees)')
        >>> plt.margins(0.05, 0.05)
        >>> plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        >>> plt.tight_layout()
        >>> plt.show()
    
    '''
    suitable_for_FFT = False
    
    def __init__(self):
        '''
        Initialize rotation matrices.
        '''
        raise RuntimeError('Base class cannot be instanciated')
    
    def _reset(self):
        '''
        Initialize the rotation matrix generator.
        '''
        self.R = []
        self.R_hash = None
        self.R_chi_matrix = None
        self.desc = type(self).__name__
        self.desc_short = type(self).__name__.replace('RotationProtocol', '')
    
    def __len__(self):
        return len(self.R)
    
    def apply(self, coords, idx):
        '''
        Apply a rotation onto a set of coordinates.
        
        :param coords: atomic coordinates
        :type  coords: np.array[:,3]
        :param idx: rotation matrix number
        :type  idx: int
        :returns: Rotated coordinates.
        :returntype: np.array[:,3]
        '''
        R = self.R[idx]
        return np.array([np.dot(R, coords[i]) for i in range(len(coords))])
    
    def Lattman_distance(self):
        '''
        Compute the Lattman distance between rotation angles.
        '''
        # check if Lattman distance matrix already available, if not or if
        # rotation matrices have changed, (re-)compute the distane matrix
        R_copy = np.copy(self.R)
        R_copy.flags.writeable = False
        R_hash = hash(R_copy.data)
        if self.R_hash is None or self.R_hash != R_hash or self.R_chi_matrix is None:
            self.R_hash = R_hash # update hash
            self.R_chi_matrix = Lattman_distance(self.R)
        return self.R_chi_matrix
    
    def plot(self, plot=None, metric='min', linemarker='o-', **kwargs):
        '''
        Plot the Lattman distance as a function of the rotation matrices,
        ordered by increasing distance. Call method show() on the resulting
        object to display the plot in a new window.
        
        :param plot: plotting device, default is matplotlib.pyplot
        :type  plot: :class:`matplotlib.pyplot`
        :param \*\*kwargs: optional arguments for :meth:`matplotlib.pyplot.plot`
        :type  \*\*kwargs: dict
        :returns: Plotting device.
        :returntype: :class:`matplotlib.pyplot`
        '''
        if plot is None:
            plot = plt
        plot = self.plot_add(plot, metric=metric, linemarker=linemarker, **kwargs)
        plot.xlabel('Rotation matrix')
        if metric == 'min':
            plot.ylabel('Distance to closest neighbor (degrees)')
            # plot.ylim(ymin=0)
        elif metric == 'max':
            plot.ylabel('Distance to farthest neighbor (degrees)')
            # plot.ylim(ymax=180)
        plot.title(self.desc)
        plot.grid(True)
        return plot
    
    def plot_add(self, plot, metric='min', linemarker='o-', **kwargs):
        '''
        Plot the Lattman distance. Call method show() on the resulting object
        to display the plot in a new window.
        
        :param plot: plotting device, default is matplotlib.pyplot
        :type  plot: :class:`matplotlib.pyplot`
        :param \*\*kwargs: optional arguments for :meth:`matplotlib.pyplot.plot`
        :type  \*\*kwargs: dict
        :returns: Plotting device.
        :returntype: :class:`matplotlib.pyplot`
        '''
        self.Lattman_distance()
        if metric == 'min':
            dist = np.amin(self.R_chi_matrix + 2 * np.pi * np.eye(len(self)), 1)
        elif metric == 'max':
            dist = np.amax(self.R_chi_matrix, 1)
        else:
            raise RuntimeError('Unknown parameter metric="{}"'.format(metric))
        plot.plot(range(len(dist)), np.sort(dist) * 180 / np.pi, linemarker, **kwargs)
        return plot


class RotationProtocolFibonacciSpin(RotationProtocol):
    '''
    Compute unit vectors uniformly distributed on a sphere using a Fibonacci
    generative spiral, with an additional spin.
    
    :param N: number of vectors (range 150-300)
    :type  N: int
    :param spins: number of spins around theta (range 1-5)
    :type  spins: int
    :param fun: function use to compute the coordinates of a spherical lattice
       (optional), default is :func:`spherical_lattice_default`
    :type  fun: function
    :param R0: custom pre-rotation
    :type  R0: np.array[3,3]
    '''
    suitable_for_FFT = True
    
    def __init__(self, N, spins, fun=None, R0=None):
        if fun is None:
            fun = spherical_lattice_default
        
        self.desc = '{}: N={}, S={}, F={}'.format(
           type(self).__name__, N, spins, fun.__name__)
        self.desc_short = '{}: N={}, S={}'.format(
           type(self).__name__.replace('RotationProtocol', ''), N, spins)
        
        self.R = []
        if R0 is None:
            R0 = np.eye(3)
        for point in fun(N):
            psi, theta, phi = get_euler_angles(fixed_point=[1, 0, 0],
                                               rotate_point=point)
            R2 = calculate_rotation_matrix(psi, theta, phi)
            for spin in range(0, 360, 360 // spins):
                R1 = calculate_rotation_matrix(0, spin, 0)
                R = np.dot(R2, np.dot(R1, R0))
                self.R.append(R)
        
        self.R_hash = None
        self.R_chi_matrix = None


class RotationProtocolFibonacci(RotationProtocolFibonacciSpin):
    '''
    Compute unit vectors uniformly distributed on a sphere using a Fibonacci
    generative spiral.
    
    :param N: number of vectors
    :type  N: int
    :param fun: function use to compute the coordinates of a spherical lattice
       (optional), default is :func:`spherical_lattice_default`
    :type  fun: function
    :param R0: custom pre-rotation
    :type  R0: np.array[3,3]
    '''
    suitable_for_FFT = True
    
    def __init__(self, N, fun=None, R0=None):
        super(RotationProtocolFibonacci, self).__init__(N, 1, fun=fun, R0=R0)
        self.desc = ','.join(x for x in self.desc.split(',') if not x.startswith(' S='))
        self.desc_short = ','.join(x for x in self.desc_short.split(',') if not x.startswith(' S='))


class RotationProtocolEulerUniform(RotationProtocol):
    '''
    Compute rotation matrices based on uniformly sampled Euler angles
    :math:`(\\theta_1, \\theta_2, \\theta_3)` with
    :math:`\\theta_1 \\in [-\\pi, \\pi[`,
    :math:`\\theta_2 \\in [-\\pi, \\pi[`,
    :math:`\\theta_3 \\in [0, \\pi]`.
    
    The resulting rotations are not suitable for FFT correlation.
    
    :param level: number of samples for theta, in the range [3, 10]
    :type  level: int
    :param R0: custom pre-rotation
    :type  R0: np.array[3,3]
    '''
    suitable_for_FFT = 'uniformly sampled Euler angles are not uniform in SO(3)'
    
    def __init__(self, level, R0=None):
        self.desc = '{}: L={}'.format(type(self).__name__, level)
        self.desc_short = self.desc.replace('RotationProtocol', '')
        
        self.R = []
        if R0 is None:
            R0 = np.eye(3)
        
        for psi in np.arange(-np.pi, np.pi, 2 * np.pi / level):
            R1 = Euler_Rodrigues_rotation_matrix([0,0,1], psi)
            for theta in np.arange(-np.pi, np.pi, 2 * np.pi / level):
                R2 = Euler_Rodrigues_rotation_matrix([0,1,0], theta)
                for phi in np.arange(0, np.pi + np.pi / level / 2., np.pi / level):
                    R3 = Euler_Rodrigues_rotation_matrix([1,0,0], phi)
                    R = np.dot(R3, np.dot(R2, np.dot(R1, R0)))
                    self.R.append(R)
        
        self.R_hash = None
        self.R_chi_matrix = None


class RotationProtocolEulerNonRedundant(RotationProtocol):
    '''
    Compute rotation matrices based on uniformly sampled Euler angles
    :math:`(\\theta_1, \\theta_2, \\theta_3)` with
    :math:`\\theta_1 \\in [-\\pi, \\pi[`,
    :math:`\\theta_2 \\in [-\\pi, \\pi[`,
    :math:`\\theta_3 \\in [0, \\pi]`. Degenerate rotations are filtered out.
    
    :param step: angle increment (degrees), in the range [10, 20]
    :type  step: float
    :param degeneracy_threshold: angle threshold (degrees) below which two
       rotations are considered degenerate, in the range [0, 2]
    :type  degeneracy_threshold: float
    :param R0: custom pre-rotation
    :type  R0: np.array[3,3]
    
    Implementation adapted from the following source:
    
    * Henry Gabb, Richard Jackson and Michael Sternberg, Modelling protein
      docking using shape complementarity, electrostatics and biochemical
      information, *Journal of Molecular Biology* **1997**, *272(1)*: 106--120.
      DOI: `10.1006/jmbi.1997.1203 <https://doi.org/10.1006/jmbi.1997.1203>`_
    '''
    suitable_for_FFT = True
    
    def __init__(self, step, degeneracy_threshold=None, R0=None, verbose=False):
        self.desc = '{}: A={}'.format(type(self).__name__, step)
        self.desc_short = self.desc.replace('RotationProtocol', '')
        
        self.R = []
        if R0 is None:
            R0 = np.eye(3)
        
        for psi in np.arange(0, 360, step):
            R1 = Euler_Rodrigues_rotation_matrix([0,0,1], psi*np.pi/180)
            for theta in np.arange(0, 360, step):
                R2 = Euler_Rodrigues_rotation_matrix([0,1,0], theta*np.pi/180)
                for phi in np.arange(0, 180, step):
                    R3 = Euler_Rodrigues_rotation_matrix([1,0,0], phi*np.pi/180)
                    R = np.dot(R3, np.dot(R2, np.dot(R1, R0)))
                    self.R.append(R)
        
        self.R_hash = None
        self.R_chi_matrix = None
        
        # extremely slow redundancy check, need to find a better algorithm
        if degeneracy_threshold is not None and degeneracy_threshold > 0:
            if verbose:
                print('Computing Lattman distance...')
            self.Lattman_distance()
            while True:
                smallest = np.min(np.amin(self.R_chi_matrix + 2 * np.pi * np.eye(len(self)), 1))
                if smallest * 180/np.pi > degeneracy_threshold:
                    break
                if verbose:
                    print(smallest * 180/np.pi, len(self))
                # find the most degenerate state
                idx = [i[0] for i in sorted(enumerate(np.sort(self.R_chi_matrix, axis=1).tolist()), key=lambda x:x[1])][0]
                del self.R[idx] # delete rotation
                self.R_chi_matrix = np.delete(self.R_chi_matrix, (idx), axis=0) # delete row
                self.R_chi_matrix = np.delete(self.R_chi_matrix, (idx), axis=1) # delete column
        
            # update hash
            R_copy = np.copy(self.R)
            R_copy.flags.writeable = False
            self.R_hash = hash(R_copy.data)


class RotationProtocolRandomUniform(RotationProtocol):
    '''
    Compute rotation matrices based on a random uniform sampling of Euler
    angles :math:`(\\theta_1, \\theta_2, \\theta_3)` with
    :math:`\\theta_1 \\in [-\\pi, \\pi[`,
    :math:`\\theta_2 \\in [-\\pi, \\pi[`,
    :math:`\\theta_3 \\in [0, \\pi]`.
    
    The resulting rotations are not suitable for FFT correlation.
    
    :param N: number of points
    :type  N: int
    :param R0: custom pre-rotation
    :type  R0: np.array[3,3]
    '''
    suitable_for_FFT = 'random Euler angles are not uniform in SO(3)'
    
    def __init__(self, N, seed=0, R0=None):
        self.desc = '{}: N={}'.format(type(self).__name__, N)
        self.desc_short = self.desc.replace('RotationProtocol', '')
        
        self.R = []
        if R0 is None:
            R0 = np.eye(3)
        
        np.random.seed(seed)
        psi_list   = np.pi * (2 * np.random.rand(N) - 1)
        theta_list = np.pi * (2 * np.random.rand(N) - 1)
        phi_list   = np.pi * np.random.rand(N)
        for psi, theta, phi in zip(psi_list, theta_list, phi_list):
            R1 = Euler_Rodrigues_rotation_matrix([0,0,1], psi)
            R2 = Euler_Rodrigues_rotation_matrix([0,1,0], theta)
            R3 = Euler_Rodrigues_rotation_matrix([1,0,0], phi)
            R = np.dot(R3, np.dot(R2, np.dot(R1, R0)))
            self.R.append(R)
        
        self.R_hash = None
        self.R_chi_matrix = None


class RotationProtocolLattman(RotationProtocol):
    '''
    Compute rotation matrices using the orthogonal Lattman angles.
    
    The resulting rotations are not suitable for FFT correlation.
    
    :param level: number of samples for theta, in the range [3, 10]
    :type  level: int
    :param R0: custom pre-rotation
    :type  R0: np.array[3,3]
    
    Implementation adapted from the following source:
    
    * Eaton Lattman, Optimal sampling of the rotation function,
      *Acta Crystallographica Section B* **1972**, *B28*: 1065--1068. DOI:
      `10.1107/S0567740872003723 <https://doi.org/10.1107/S0567740872003723>`_
    '''
    suitable_for_FFT = 'uniformly sampled Lattman angles are not uniform in SO(3)'
    
    def __init__(self, level, R0=None):
        self.desc = '{}: L={}'.format(type(self).__name__, level)
        self.desc_short = self.desc.replace('RotationProtocol', '')
        
        self.R = []
        if R0 is None:
            R0 = np.eye(3)
        
        for thetaPlus, theta, thetaMinus in Lattman(level):
            R1 = Euler_Rodrigues_rotation_matrix([0,0,1], thetaMinus)
            R2 = Euler_Rodrigues_rotation_matrix([0,1,0], theta)
            R3 = Euler_Rodrigues_rotation_matrix([1,0,0], thetaPlus)
            R = np.dot(R3, np.dot(R2, np.dot(R1, R0)))
            self.R.append(R)
        
        self.R_hash = None
        self.R_chi_matrix = None


def Lattman(N):
    result = []
    step = np.pi / N
    
    def sampleTheta(step):
        return np.arange(0, np.pi*1.001, step)
    
    def sampleThetaPlus(theta, step):
        if np.isclose(np.sin(theta/2 + np.pi/2), 0, rtol=0, atol=1e-12):
            return [0] # degenerate, any value will do
        else:
            return np.arange(0, 4 * np.pi, step/np.sin(theta/2 + np.pi/2))
    
    def sampleThetaMinus(theta, step):
        if np.isclose(np.sin(np.pi/2 - theta/2), 0, rtol=0, atol=1e-12):
            return [0] # degenerate, any value will do
        else:
            return np.arange(0, 2 * np.pi, step/np.sin(np.pi/2 - theta/2))
    
    for theta in sampleTheta(step):
        for thetaPlus in sampleThetaPlus(theta, step):
            for thetaMinus in sampleThetaMinus(theta, step):
                result.append((thetaPlus, theta, thetaMinus))
    return result


class RotationProtocolExact(RotationProtocol):
    '''
    Uniform sampling in SO(3) based on symmetry groups.
    
    The resulting rotations are not suitable for FFT correlation (N=60 max).
    
    :param symmetry: symmetry group (Ih=dodecahedron=icosahedron has N=60,
       Oh=octahedron=cube has N=24, Td=tetrahedron has N=12)
    :type  symmetry: str
    '''
    suitable_for_FFT = 'number of rotations is too small'
    
    axis_angle = [
       np.array([1.000000, 0.000000, 0.000000, 0.000000,
                 1.000000, 0.000000, 0.000000, 3.141593,
                 0.000000, 0.000000,-1.000000, 3.141593,
                 0.000000,-1.000000, 0.000000, 3.141593,
                -0.577350, 0.577350,-0.577350, 2.094395,
                 0.577350,-0.577350, 0.577350, 2.094395,
                -0.577350,-0.577350, 0.577350, 2.094395,
                 0.577350, 0.577350,-0.577350, 2.094395,
                 0.577350,-0.577350,-0.577350, 2.094395,
                -0.577350, 0.577350, 0.577350, 2.094395,
                 0.577350, 0.577350, 0.577350, 2.094395,
                -0.577350,-0.577350,-0.577350, 2.094395]).reshape([12,4]),
       np.array([1.000000, 0.000000, 0.000000, 0.000000,
                 0.707107, 0.000000, 0.707107, 3.141593,
                 0.707107, 0.000000,-0.707107, 3.141593,
                 0.000000,-1.000000,-1.000000, 3.141593,
                 0.707107, 0.707107, 0.000000, 3.141593,
                 0.000000,-1.000000, 1.000000, 3.141593,
                 0.707107,-0.707107, 0.000000, 3.141593,
                 0.000000,-1.000000, 0.000000, 3.141593,
                -0.577350,-0.577350, 0.577350, 2.094395,
                -0.577350, 0.577350, 0.577350, 2.094395,
                 0.577350,-0.577350,-0.577350, 2.094395,
                 0.577350, 0.577350,-0.577350, 2.094395,
                -0.577350, 0.577350,-0.577350, 2.094395,
                 0.577350,-0.577350, 0.577350, 2.094395,
                 0.577350, 0.577350, 0.577350, 2.094395,
                -0.577350,-0.577350,-0.577350, 2.094395,
                 1.000000, 0.000000, 0.000000, 3.141593,
                 0.000000, 0.000000,-1.000000, 3.141593,
                 1.000000, 0.000000, 0.000000, 1.570796,
                 0.000000, 0.000000,-1.000000, 1.570796,
                -1.000000, 0.000000, 0.000000, 1.570796,
                 0.000000, 0.000000, 1.000000, 1.570796,
                 0.000000, 1.000000, 0.000000, 1.570796,
                 0.000000,-1.000000, 0.000000, 1.570796]).reshape([24,4]),
       np.array([1.000000, 0.000000, 0.000000, 0.000000,
                 0.309017, 0.809017, 0.500000, 3.141593,
                 0.500000,-0.309017,-0.809017, 3.141593,
                 1.000000, 0.000000, 0.000000, 3.141593,
                 0.500000,-0.309017, 0.809017, 3.141593,
                 0.309017, 0.809017,-0.500000, 3.141593,
                 0.309017,-0.809017,-0.500000, 3.141593,
                 0.309017,-0.809017, 0.500000, 3.141593,
                 0.809017,-0.500000,-0.309017, 3.141593,
                 0.000000,-1.000000, 0.000000, 3.141593,
                 0.809017,-0.500000, 0.309017, 3.141593,
                 0.500000, 0.309017, 0.809017, 3.141593,
                 0.000000, 0.000000,-1.000000, 3.141593,
                 0.809017, 0.500000, 0.309017, 3.141593,
                 0.809017, 0.500000,-0.309017, 3.141593,
                 0.500000, 0.309017,-0.809017, 3.141593,
                 0.577350,-0.577350, 0.577350, 2.094395,
                 0.000000, 0.525731,-0.850651, 2.513274,
                 0.850651,-0.000000,-0.525731, 2.513274,
                -0.934172, 0.356822, 0.000000, 2.094395,
                -0.000000,-0.525731, 0.850651, 1.256637,
                -0.850651, 0.000000, 0.525731, 2.513274,
                 0.000000,-0.525731, 0.850651, 2.513274,
                -0.850651, 0.000000, 0.525731, 1.256637,
                 0.850651,-0.000000,-0.525731, 1.256637,
                 0.934172,-0.356822, 0.000000, 2.094395,
                 0.000000, 0.525731,-0.850651, 1.256637,
                -0.577350, 0.577350,-0.577350, 2.094395,
                -0.000000,-0.934172, 0.356822, 2.094395,
                 0.525731, 0.850651, 0.000000, 2.513274,
                 0.850651,-0.000000, 0.525731, 2.513274,
                -0.850651, 0.000000,-0.525731, 1.256637,
                -0.525731,-0.850651, 0.000000, 1.256637,
                -0.850651, 0.000000,-0.525731, 2.513274,
                -0.525731,-0.850651, 0.000000, 2.513274,
                 0.000000, 0.934172,-0.356822, 2.094395,
                 0.525731, 0.850651, 0.000000, 1.256637,
                 0.850651, 0.000000, 0.525731, 1.256637,
                -0.000000,-0.934172,-0.356822, 2.094395,
                 0.000000, 0.525731, 0.850651, 2.513274,
                 0.000000,-0.525731,-0.850651, 2.513274,
                -0.000000,-0.525731,-0.850651, 1.256637,
                -0.000000, 0.525731, 0.850651, 1.256637,
                -0.000000, 0.934172, 0.356822, 2.094395,
                 0.577350,-0.577350,-0.577350, 2.094395,
                -0.577350, 0.577350, 0.577350, 2.094395,
                -0.934172,-0.356822, 0.000000, 2.094395,
                -0.577350,-0.577350, 0.577350, 2.094395,
                 0.577350, 0.577350,-0.577350, 2.094395,
                 0.934172, 0.356822, 0.000000, 2.094395,
                 0.577350, 0.577350, 0.577350, 2.094395,
                -0.577350,-0.577350,-0.577350, 2.094395,
                -0.356822, 0.000000,-0.934172, 2.094395,
                 0.356822, 0.000000, 0.934172, 2.094395,
                -0.356822, 0.000000, 0.934172, 2.094395,
                 0.356822,-0.000000,-0.934172, 2.094395,
                 0.525731,-0.850651,-0.000000, 1.256637,
                -0.525731, 0.850651,-0.000000, 1.256637,
                -0.525731, 0.850651, 0.000000, 2.513274,
                 0.525731,-0.850651,-0.000000, 2.513274]).reshape([60,4])
    ]
    
    def __init__(self, symmetry='dodecahedron'):
        '''
        Compute a rotation matrix.
        
        :param symmetry: symmetry group (Ih=dodecahedron=icosahedron has N=60,
           Oh=octahedron=cube has N=24, Td=tetrahedron has N=12)
        :type  symmetry: str
        '''
        if symmetry in ('Ih', 'dodecahedron', 'dodecahedral',
                        'icosahedron', 'icosahedral'):
            group = 2
        elif symmetry in ('Oh', 'octahedron', 'octahedral', 'cube', 'cubic'):
            group = 1
        elif symmetry in ('Td', 'tetrahedron', 'tetrahedral'):
            group = 0
        else:
            raise ValueError('Molecular symmetry group "{}" unknown'.format(
                                                                     symmetry))
        
        self.desc = '{}, N={}'.format(symmetry, len(self.axis_angle[group]))
        self.desc_short = self.desc
        
        self.R = []
        for x, y, z, theta in self.axis_angle[group]:
            self.R.append(Euler_Rodrigues_rotation_matrix([x, y, z], theta))
        
        self.R_hash = None
        self.R_chi_matrix = None


def Euler_Rodrigues_rotation_matrix(axis, theta):
    '''
    Compute a rotation matrix.
    
    :param axis: rotation axis
    :type  axis: list
    :param theta: rotation angle in radians
    :type  theta: float
    :returns: Rotation matrix.
    :returntype: np.array[3,3]
    '''
    # adapted from Gaurav Dhama <https://stackoverflow.com/a/38792397>
    axis = np.array(axis) / np.linalg.norm(axis)  # need a unit vector
    a = np.cos(theta / 2)
    b, c, d = axis * np.sin(theta / 2)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2 * (bc+ad), 2 * (bd-ac)],
                     [2 * (bc-ad), aa+cc-bb-dd, 2 * (cd+ab)],
                     [2 * (bd+ac), 2 * (cd-ab), aa+dd-bb-cc]])


def rotate_around_axis(coordinates, axis, angle):
    '''
    Perform a rotation around an axis.

    :param coordinates: coordinate(s)
    :type  coordinates: np.array
    :param axis: rotation axis
    :type  axis: list
    :param theta: rotation angle in radians
    :type  theta: float
    :returns: Rotated coordinates.
    :returntype: np.array
    '''
    shape = np.copy(coordinates).shape
    if shape == (3,):
        coord = [coordinates]
    elif shape[1] == 3:
        coord = np.copy(coordinates)
    else:
        raise ValueError('Wrong array shape: ({})'.format(shape))
    R = Euler_Rodrigues_rotation_matrix(axis, angle)
    coord = np.array([np.dot(R, v) for v in coord])
    if shape == (3,):
        coord = coord[0]
    return coord


