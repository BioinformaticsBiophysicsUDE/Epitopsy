
:mod:`UtilityClasses` --- Various utility functions
===================================================

.. module:: UtilityClasses
   :synopsis: Various utility functions.
.. moduleauthor:: Christoph Wilms <christoph.wilms@uni-due.de>
.. sectionauthor:: Jean-Noel Grad <jean-noel.grad@uni-due.de>


This module provides various utility functions.


.. _UtilityClasses-syntax:

Module Syntax
-------------

Empty.

.. _contents-of-module-UtilityClasses:

Module Contents
---------------

.. class:: Progress_bar

    This class builds a nice progress bar::

            0%                                    100%
            #                                        #
            #********                                #
            #********************                    #
            
    There are 50 characters between the #.


    .. method:: add()

        Docstring missing.

.. class:: Fix_xxmer

    This class fixes xxmers, where xx copies of one protein are in one pdb.


    .. method:: hack_linenumbering()

        Docstring missing.

    .. method:: hack_chain_ids()

        Make double occuring chain ids unique!


.. class:: ESP_Profile_Manager

    This class writes ESP-profiles to file.


    .. method:: write_profile(box, indices, id, filename, append = False)

        :param box: array with 3 dimensions
        :param indices: tupel of arrays, which you get for example with :func:`Numpy.nonzero(a==1)`
        :param id: identifier for one ESP-profile in case you want to append it to an
                 existing file
        :param append: ``True`` or ``False``

        The formatting looks like this (seperated by \\t)::

            X        [coord_1    coord_2_... coord_n]
            Y        [coord_1    coord_2_... coord_n]
            Z        [coord_1    coord_2_... coord_n]
            [ID_1]   [phi_1        phi_2 ...   phi_n]


    .. method:: read_hull_coordinates(hull_coordinates_file, phi_values = False, id = 1)

        Docstring missing.

    .. method:: write_hull_coordinates(hull_coordinates, filename = 'hull_coordinates.txt')

        Docstring missing.

