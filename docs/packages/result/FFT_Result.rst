
:mod:`FFT_Result` --- YYYYYYY
======================================================

.. module:: FFT_Result
   :synopsis: YYYYYYYY.
.. moduleauthor:: Christoph Wilms <christoph.wilms@uni-due.de>
.. sectionauthor:: Jean-Noel Grad <jean-noel.grad@uni-due.de>


This module provides YYYYY.


.. _FFT_Result-syntax:

Module Syntax
-------------

Empty.

.. _contents-of-module-FFT_Result:

Module Contents
---------------

.. class:: FFT_Result(object)

    This class stores the results from FFT Correlation Docking.
    It stores them in a list, which contains a dictionary for every entry.
    
    The format looks like this:
    
    x_coord y_coord z_coord phi theta psi score.

    .. method:: _add_score(x, y, z, phi, theta, psi, score)

        This method adds a score to the list.

    .. method:: _sort()

        Sort the results by their score, with the highest beeing the first.

    .. method:: find_scores(box, phi, theta, psi, number_of_elements, dxbox = None)

        This method finds the highest scoring elements from the box array. It
        needs the dxbox to transform the coordinates to real space, if none is
        given, it will store the box indices.

    .. method:: _return_highest_score(mod_box)

        This method returns the highest score + the box coordinates of this 
        score. For further calls, it sets the highest score to the lowest, so 
        that the next highest score can be returned the next time.

    .. method:: write_to_disk(filename)

        Write results to disk. Before doing so, the list will be sorted by the
        score.

    .. method:: read_from_filename(filename)

        Read the results from disk.

    .. method:: make_pdb_results(pdb_path, storage_dir, new_name, best_x)

        This method takes the 'best_x' results and moves the pdb structure 
        to these positions. It stores the moved and rotated structures in the
        specified directory. The name of the pdb works like this:
        'new_name' + '_number.pdb', with number indexing the results from 
        0 as the best to the end.
        
        This method assumes that the given coordinates are in real space!!!

