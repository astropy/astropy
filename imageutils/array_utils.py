# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module includes helper functions for array operations.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import numpy as np

__all__ = ['extract_array_2d', 'add_array_2d', 'subpixel_indices',
           'fix_prf_nan']


def _get_slices(large_array_shape, small_array_shape, position):
    """
    Get slices for the overlapping part of a small and a large array.

    Given a certain position of the center of the small array, with
    respect to the large array, four slices are computed, which can be
    used to extract, add or subtract the small array at the given
    position. This function takes care of the correct behavior at the
    boundaries, where the small array is cut of appropriately.

    Parameters
    ----------
    large_array_shape : tuple
        Shape of the large array.
    small_array_shape : tuple
        Shape of the small array.
    position : tuple, (x, y)
        Position of the small array's center, with respect
        to the large array.

    Returns
    -------
    s_y : slice
        Slice in y direction for the large array.
    s_x : slice
        Slice in x direction for the large array.
    b_y : slice
        Slice in y direction for the small array.
    b_x : slice
        Slice in x direction for the small array.
    """
    # Get edge coordinates
    y_min = int(position[1] + 0.5) - small_array_shape[0] // 2
    x_min = int(position[0] + 0.5) - small_array_shape[1] // 2
    y_max = int(position[1] + 0.5) + small_array_shape[0] // 2 + 1
    x_max = int(position[0] + 0.5) + small_array_shape[1] // 2 + 1

    # Set up slices in x direction
    s_x = slice(max(0, x_min), min(large_array_shape[1], x_max))
    b_x = slice(max(0, -x_min), min(large_array_shape[1] - x_min,
                                    x_max - x_min))

    # Set up slices in y direction
    s_y = slice(max(0, y_min), min(large_array_shape[0], y_max))
    b_y = slice(max(0, -y_min), min(large_array_shape[0] - y_min,
                                    y_max - y_min))
    return s_y, s_x, b_y, b_x


def extract_array_2d(array_large, shape, position):
    """
    Extract smaller array of given shape and position out of a larger array.

    Parameters
    ----------
    array_large : ndarray
        Array to extract another array from.
    shape : tuple
        Shape of the extracted array.
    position : tuple, (x, y)
        Position of the small array's center, with respect
        to the large array.

    Examples
    --------
    We consider a large array of zeros with the shape 11x10 and a small
    array of ones with a shape of 3x4:

    >>> import numpy as np
    >>> from photutils.utils import extract_array_2d
    >>> large_array = np.zeros((11, 10))
    >>> large_array[4:9, 4:9] = np.ones((5, 5))
    >>> extract_array_2d(large_array, (3, 4), (7, 7))
    array([[ 1.,  1.,  1.,  1.,  0.],
           [ 1.,  1.,  1.,  1.,  0.],
           [ 1.,  1.,  1.,  1.,  0.]])
    """
    # Check if larger array is really larger
    if array_large.shape >= shape:
        s_y, s_x, _, _ = _get_slices(array_large.shape, shape, position)
        return array_large[s_y, s_x]
    else:
        raise Exception("Can't extract array. Shape too large.")


def add_array_2d(array_large, array_small, position):
    """
    Add a smaller 2D array at a given position in a larger 2D array.

    Parameters
    ----------
    array_large : ndarray
        Large array.
    array_small : ndarray
        Small array to add.
    position : tuple, (x, y)
        Position of the small array's center, with respect
        to the large array.

    Examples
    --------
    We consider a large array of zeros with the shape 5x5 and a small
    array of ones with a shape of 3x3:

    >>> import numpy as np
    >>> from photutils.utils import add_array_2d
    >>> large_array = np.zeros((5, 5))
    >>> small_array = np.ones((3, 3))
    >>> add_array_2d(large_array, small_array, (2, 1))
    array([[ 0.,  1.,  1.,  1.,  0.],
           [ 0.,  1.,  1.,  1.,  0.],
           [ 0.,  1.,  1.,  1.,  0.],
           [ 0.,  0.,  0.,  0.,  0.],
           [ 0.,  0.,  0.,  0.,  0.]])
    """
    # Check if larger array is really larger
    if array_large.shape >= array_small.shape:
        s_y, s_x, b_y, b_x = _get_slices(array_large.shape,
                                         array_small.shape, position)
        array_large[s_y, s_x] += array_small[b_y, b_x]
        return array_large
    else:
        raise Exception("Can't add array. Small array too large.")


def subpixel_indices(position, subsampling):
    """
    Convert decimal points to indices, given a subsampling factor.

    Parameters
    ----------
    position : ndarray [x, y]
        Positions in pixels.
    subsampling : int
        Subsampling factor per pixel.
    """
    # Get decimal points
    x_frac, y_frac = np.modf((position[0] + 0.5,
                             position[1] + 0.5))[0]

    # Convert to int
    x_sub = np.int(x_frac * subsampling)
    y_sub = np.int(y_frac * subsampling)
    return x_sub, y_sub


def fix_prf_nan(extracted_prf, prf_nan):
    """
    Fix NaN values in an extracted PRF image.

    The NaN values are fixed by replacing with the mirrored
    value with respect to the PRF's center

    Parameters
    ----------
    extracted_prf : array
        PRF array to be fixed.
    prf_nan : array, bool
        Mask indicating where NaN values are present
        ``extracted_prf``.
    """
    # Allow at most 3 NaN values to prevent the unlikely case,
    # that the mirrored values are also NaN.
    y_nan_coords, x_nan_coords = np.nonzero(prf_nan)
    for y_nan, x_nan in zip(y_nan_coords, x_nan_coords):
        if not np.isnan(extracted_prf[-y_nan - 1, -x_nan - 1]):
            extracted_prf[y_nan, x_nan] = \
                extracted_prf[-y_nan - 1, -x_nan - 1]
        elif not np.isnan(extracted_prf[y_nan, -x_nan]):
            extracted_prf[y_nan, x_nan] = \
                extracted_prf[y_nan, -x_nan - 1]
        else:
            extracted_prf[y_nan, x_nan] = \
                extracted_prf[-y_nan - 1, x_nan]
    return extracted_prf
