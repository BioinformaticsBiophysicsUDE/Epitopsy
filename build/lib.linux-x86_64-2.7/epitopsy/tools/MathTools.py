'''
Created on 30.11.2011

@author: chris
'''

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


def get_euler_angles_for_equal_distributed_rotation(number_of_rotations):
    '''
    This method returns an angle distribution that is pretty uniformley
    distributed.
    '''
    # get equal distributed points at first
    equal_points = points_on_sphere(number_of_rotations)

    # this is a reference probe
    probe = [1, 0, 0]
    euler_angles = []
    for point in equal_points:
        euler_angles.append(get_euler_angles(fixed_point = probe,
                                             rotate_point = point))

    return euler_angles


def points_on_sphere(number_of_points):
    '''
    This code is based on the golden section (by Patrick Boucher 2006).
    '''
    points = []

    inc = np.pi * (3 - np.sqrt(5))
    off = 2 / float(number_of_points)
    for k in range(0, number_of_points):
        y = k * off - 1 + (off / 2)
        r = np.sqrt(1 - y * y)
        phi = k * inc
        points.append([np.cos(phi) * r, y, np.sin(phi) * r])

    return points

def get_euler_angles(fixed_point, rotate_point):
    '''
    This method returns a rotation matrix which rotates the 'rotate_point'
    onto the fixed_point.
    '''
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
    equal_points = points_on_sphere(n_rotations)

    # this is a reference probe
    probe = np.array([1., 0., 0.])
    angle_set = []
    for point in equal_points:
        point = np.array(point)
        if np.linalg.norm(point - probe) <= max_dist:
            angle_set.append(get_euler_angles(fixed_point = probe,
                    rotate_point = point))

    return angle_set


def fix_grid_size(box_dim):
    '''
    Fix the grid size.

    Args:
        box_dim -> list with dimensions of an APBS box

    Returns:
        A numpy array with fixed grid dimensions.
    '''
    def calculate_valid_dimension(c, nlev=4):
        '''
        Due to a multilevel approach APBS requires the grid to be of certain sizes.
        (See APBS manual for more information)

        self method ensures,  that chosen grid dimensions meet these requirements.
        Current grid dimensions will be enlarged accordingly.

        Args:
            c -> Test grid dimension.
            nlev -> ?

        Returns:
            Integer number that has the correct dimension.
        '''
        return int((c * (np.power(2.0, nlev + 1)) + 1))

    # validate x dimension
    fixed_dimension = []
    for i in range(len(box_dim)):
        c = 0
        while(box_dim[i] > calculate_valid_dimension(c)):
            c = c + 1
        fixed_dimension.append(calculate_valid_dimension(c))

    return np.array(fixed_dimension)

