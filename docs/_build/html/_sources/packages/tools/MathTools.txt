
:mod:`MathTools` --- Mathematical tools
=======================================

.. module:: MathTools
   :synopsis: Mathematical tools.
.. moduleauthor:: Christoph Wilms <christoph.wilms@uni-due.de>
.. sectionauthor:: Jean-Noel Grad <jean-noel.grad@uni-due.de>


This module provides Mathematical tools.


.. _MathTools-syntax:

Module Syntax
-------------

Empty.

.. _contents-of-module-MathTools:

Module Contents
---------------

.. function:: calculate_rotation_matrix(phi, theta, psi)

    .. seealso::
            Euler Angles

            http://mathworld.wolfram.com/EulerAngles.html

            * phi = 0 - 360 about the z-axis
            * theta = 0 - 180 about the new x'-axis
            * psi = 0 - 360 about the new z'-axis

            Same formalism as in biopython Bio.PDB.Superimposer

                :math:`\rightarrow` angles can be derived from Superimposer.rotran!


.. function:: get_euler_angles_from_rotation_matrix(euler_rotation)

    This method returns the angles from a rotation matrix, derived from
    euler angles. Because only angles in the range [0,360] (and [0,180])
    are used, combinations due to the 2 pi periodicity are possible, but
    are of no interest.

    To get the angles of a combined rotation by two sets of
    :math:`[\phi_1, \theta_1, \psi_1]` and
    :math:`[\phi_2, \theta_2, \psi_2]`
    one has to multiply the rotationmatrices with 'np.dot(...)'!


.. function:: get_euler_angles_for_equal_distributed_rotation(number_of_rotations)

    This method returns an angle distribution that is pretty uniformley
    distributed.


.. function:: points_on_sphere(number_of_points)

    This code is based on the golden section (by Patrick Boucher 2006).


.. function:: get_euler_angles(fixed_point, rotate_point)

    This method returns a rotation matrix which rotates the 'rotate_point'
    onto the fixed_point.


.. function:: get_neighbor_angle_set(n_rotations, max_dist)

    This method returns euler angles for rotations, that rotate the object
    to a oritentation, where the distance displacement is less than max_dist.
    Notice, that the returned angle set does not have the specified number
    of rotations!


.. function:: fix_grid_size(box_dim)

    Fix the grid size.

    :param box_dim: list with dimensions of an APBS box

    :returns: a numpy array with fixed grid dimensions.


.. function:: calculate_valid_dimension(c, nlev=4)

    Due to a multilevel approach APBS requires the grid to be of certain sizes.
    (See APBS manual for more information.)

    Self method ensures, that chosen grid dimensions meet these requirements.
    Current grid dimensions will be enlarged accordingly.

    :param c: test grid dimension
    :param nlev: number of levels

    :returns: integer number which has the correct dimension

