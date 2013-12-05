
:mod:`calc_vol_and_surface` --- Calculation of Volume and Surface
=================================================================

.. module:: calc_vol_and_surface
   :synopsis: Calculation of Volume and Surface.
.. moduleauthor:: Christoph Wilms <christoph.wilms@uni-due.de>
.. sectionauthor:: Jean-Noel Grad <jean-noel.grad@uni-due.de>


This module provides direct coupling analysis operations.


.. _calc_vol_and_surface-syntax:

Module Syntax
-------------

Empty.

.. _contents-of-module-calc_vol_and_surface:

Module Contents
---------------

.. function:: check_indices(indices, box, raise_error)

    Docstring missing.

.. function:: process_vol(energy_box, vdw_box, energy_cutoff, result_dict, raise_error, vdw_type)

    This function is not well designed.

    :param probability_box: DXBox contains the energies.
    :param vdw_box: flooded VDWBox.
    :param energy_cutoff: count grid points with an energy < limit
    :param result_dict: dictionary

    :returns: updated dictionary
    
.. function:: process_surface(energy_box, surface_index, energy_cutoff, result_dict, vdw_type)

    Docstring missing.

.. function:: calc_vol_and_surface(pdb_path, energy_cutoff, zipped=True, result_path="result_data.txt", energy_box=None, counter_box=None, raise_error=True, conc=1.)

    :param pdb_path: path to the pdb file
    :param energy_cutoff: positive float
    :param zipped: are the matrices gzipped or not? Default is True.
    :param result_path: path to which the results are to be written.

    :returns: a dictionary with following keys:
	
        * 'total_volume'
        * 'fav_volume'
        * 'unfav_volume'
        * 'neutral_volume'
        * 'total_surface'
        * 'fav_surface'
        * 'unfav_surface'
        * 'neutral_surface'
    
