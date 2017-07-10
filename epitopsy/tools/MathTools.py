# @date Created on 30.11.2011
# @author: chris

import numpy as np

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
#        print('{0} should be in the range [0, 360] and is {1}'.format('phi', phi))
        phi = np.deg2rad(phi)
    if 0 <= theta <= 180:
        theta = np.deg2rad(theta)
    else:
#        print('{0} should be in the range [0, 180] and is {1}'.format('theta', theta))
        theta = np.deg2rad(theta)
    if 0 <= psi <= 360:
        psi = np.deg2rad(psi)
    else:
#        print('{0} should be in the range [0, 360] and is {1}'.format('psi', psi))
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

    #if psi < 0 :
    #    # add 2*np.pi
    #    psi = psi + 360
    #
    #if phi < 0:
    #    # add 2*np.pi
    #    phi = phi + 360

    return [phi, theta, psi]


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
      DOI: `10.1256/qj.05.227 <http://dx.doi.org/10.1256/qj.05.227>`_
    * Alvaro Gonzalez, Measurement of Areas on a Sphere Using Fibonacci and
      Latitude--Longitude Lattices, *Mathematical Geosciences*, **2010**,
      *42*: 49. DOI: `10.1007/s11004-009-9257-x
      <http://dx.doi.org/10.1007/s11004-009-9257-x>`_
    
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
      <http://dx.doi.org/10.1007/s11004-009-9257-x>`_
    
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
    
    * Edward B. Saff and A. B. J. Kuijlaars, Distributing many points on a
      sphere, *The Mathematical Intelligencer* **1997**, *19(1)*: 5--11.
      DOI: `10.1007/BF03024331
      <http://dx.doi.org/10.1007/BF03024331>`_
    
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
            angle_set.append(get_euler_angles(fixed_point = probe,
                    rotate_point = point))

    return angle_set


def fix_grid_size(box_dim, nlev=4):
    '''
    Due to a multilevel approach APBS requires the grid to be of certain sizes
    (see the APBS manual for more information). Valid dimensions will be as
    large or larger than **box_dim**.

    :param box_dim: proposed APBS grid dimensions
    :type  box_dim: tuple(int,int,int)
    :param nlev: depth of the multilevel hierarchy
    :type  nlev: int
    :returns: Valid APBS grid dimensions
    :returntype: :class:`numpy.ndarray[3]`
    '''
    def calculate_valid_dimension(c, nlev=4):
        '''
        Test a grid dimension c.
        
        :param c: test grid dimension
        :type  c: int
        :param nlev: depth of the multilevel hierarchy
        :type  nlev: int
        :returns: APBS grid dimension
        :returntype: int
        '''
        return int((c * (np.power(2.0, nlev + 1)) + 1))

    # validate dimensions
    fixed_dimension = []
    for i in range(len(box_dim)):
        c = 0
        while(box_dim[i] > calculate_valid_dimension(c, nlev)):
            c = c + 1
        fixed_dimension.append(calculate_valid_dimension(c, nlev))

    return np.array(fixed_dimension)

