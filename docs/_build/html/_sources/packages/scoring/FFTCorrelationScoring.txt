
:mod:`FFTCorrelationScoring` --- YYYYYYY
======================================================

.. module:: FFTCorrelationScoring
   :synopsis: YYYYYYYY.
.. moduleauthor:: Christoph Wilms <christoph.wilms@uni-due.de>
.. sectionauthor:: Jean-Noel Grad <jean-noel.grad@uni-due.de>


This module provides YYYYY.


.. _FFTCorrelationScoring-syntax:

Module Syntax
-------------

Empty.

.. _contents-of-module-FFTCorrelationScoring:

Module Contents
---------------

.. class:: FFT_correlation_scoring(Scoring)

    This class has been once written by Niko but these functions are outdated
    and newer functions will replace the old ones step by step.
    The new functions have a name like this 'this_is_a_new_function', whereas
    old functions have a name like this 'thisWouldBeAnOldFunction'.

    It is important to set the protein interior to 1 and the solvent to 0, if
    one is interested in the geometric complementarity.


    .. method:: do_fft(signal_matrix)

        Maybe zero-padding is necessary, it could be done with::

            >>> fftn(signal_matrix, [next pow of 2, next pow of 2, next pow of 2])

    .. method:: do_ifft(signal_matrix)

        Docstring missing.

    .. method:: do_fftshift(signal_matrix)

        At the moment this function is of no use! Returns the argument.

    .. method:: do_ifftshift(signal_matrix)

        At the moment this function is of no use! Returns the argument.

    .. method:: shift_fft(signal_matrix, center = None)

        Docstring missing.

    .. method:: score(signal_matrix, position = None, maxRadius = None)

        Computes the FFT-form of the comparison box and scores it against the
        reference box,  which is in a conjugated-complex FFT-form

    .. method:: determine_best_correlation(signal_matrix)

        Docstring missing.

    .. method:: determineBestCorrelationWithinRadius(signal_matrix, position, maxRadius)

        This will be very slow ... I think the intention of this function is
        to find good correlations near the current position, so that it would
        for example be possible to search around a binding pocket and not get
        distracted by another (maybe stronger) binding pocket.

    .. method:: clone()

        Docstring missing.

    .. method:: get_correlation(signal_matrix)

        This method calculates the correlation with fft's.
        In this notation matrixA is fixed, whereas matrixB is the rotated one.
        The method returns only the real values, as they are the interesting
        part and it also shifts the zero frequencies to the middle.

    .. method:: get_imag_correlation(signal_matrix)

        This method calculates the correlation with fft's.
        In this notation matrixA is fixed, whereas matrixB is the rotated one.
        The method returns only the real values, as they are the interesting
        part and it also shifts the zero frequencies to the middle.

