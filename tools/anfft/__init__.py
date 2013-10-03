#   ANFFT is an FFT package for Python, using FFTW.
#   Copyright 2010 Andrew Collette.
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
    ANFFT is an FFT package for Python, based on FFTW.

    All functions in this module take the following optional keywords:

    measure:    If True, use FFTW's planning system to determine the most
                efficient FFT strategy for this shape and type.  The first
                call will be delayed for planning, but subsequent FFTs
                of the same shape and type are much faster.  Default is False.
"""

from _guts import fftn, ifftn, fft, ifft, rfftn, irfftn, rfft, irfft


__all__ = ['fftn','ifftn','fft','ifft','rfftn','irfftn','rfft','irfft']
