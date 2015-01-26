# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module includes helper functions for array operations.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np

__all__ = ['extract_array', 'add_array', 'subpixel_indices',
           'overlap_slices']


class NoOverlapError(ValueError):
    '''Raised when determining the overlap of non-overlapping arrays.'''
    pass


def overlap_slices(large_array_shape, small_array_shape, position):
    """
    Get slices for the overlapping part of a small and a large array.

    Given a certain position of the center of the small array, with
    respect to the large array, tuples of slices are returned which can be
    used to extract, add or subtract the small array at the given
    position. This function takes care of the correct behavior at the
    boundaries, where the small array is cut of appropriately.

    Parameters
    ----------
    large_array_shape : tuple
        Shape of the large array.
    small_array_shape : tuple
        Shape of the small array.
    position : tuple of integers
        Position of the small array's center, with respect to the large array.
        Coordinates should be in the same order as the array shape.
        When determining the center for a coordinate with an even number of
        elements, the position is rounded up. So the coordinates of the
        center of an array with shape=(2,3) will be (0,1).

    Returns
    -------
    slices_large : tuple of slices
        Slices in all directions for the large array, such that
        ``large_array[slices_large]`` extracts the region of the large array
        that overlaps with the small array.
    slices_small : slice
        Slices in all directions for the small array, such that
        ``small_array[slices_small]`` extracts the region that is inside the
        large array.
    """
    if len(small_array_shape) != len(large_array_shape):
        raise ValueError("Both arrays must have the same number of dimensions.")

    if len(small_array_shape) != len(position):
        raise ValueError("Position must have the same number of dimensions as array.")

    # Get edge coordinates
    edges_min = [int(np.floor(pos + 1 - small_shape / 2.))
                 for (pos, small_shape) in zip(position, small_array_shape)]
    edges_max = [int(np.floor(pos + 1 + small_shape / 2.))
                 for (pos, small_shape) in zip(position, small_array_shape)]

    for e_max in edges_max:
        if e_max <= 0:
            raise NoOverlapError('Arrays do not overlap.')
    for e_min, large_shape in zip(edges_min, large_array_shape):
        if e_min >= large_shape:
            raise NoOverlapError('Arrays do not overlap.')

    # Set up slices
    slices_large = tuple(slice(max(0, edge_min), min(large_shape, edge_max))
                         for (edge_min, edge_max, large_shape) in
                         zip(edges_min, edges_max, large_array_shape))
    slices_small = tuple(slice(max(0, -edge_min),
                               min(large_shape - edge_min, edge_max - edge_min))
                         for (edge_min, edge_max, large_shape) in
                         zip(edges_min, edges_max, large_array_shape))

    return slices_large, slices_small


def extract_array(array_large, shape, position):
    """
    Extract smaller array of given shape and position out of a larger array.

    Parameters
    ----------
    array_large : `~numpy.ndarray`
        Array to extract another array from.
    shape : tuple
        Shape of the extracted array.
    position : tuple
        Position of the small array's center, with respect to the large array.
        Coordinates should be in the same order as the array shape.

    Returns
    -------
    array_small : `~numpy.ma.ndarray`
        The extracted array. Values that do not overlap with ``array_large``
        are masked.

    Examples
    --------
    We consider a large array with the shape 11x10, from which we extract
    a small array of shape 3x5:

    >>> import numpy as np
    >>> from astropy.nddata.utils import extract_array
    >>> large_array = np.arange(110).reshape((11, 10))
    >>> large_array[4:9, 4:9] = np.ones((5, 5))
    >>> extract_array(large_array, (3, 5), (7, 7))
    masked_array(data =
     [[1 1 1 1 69]
     [1 1 1 1 79]
     [1 1 1 1 89]],
                 mask =
     [[False False False False False]
     [False False False False False]
     [False False False False False]],
           fill_value = 999999)

    """
    large_slices, small_slices = overlap_slices(array_large.shape,
                                                shape, position)
    extracted_array = np.ma.zeros(shape, dtype=array_large.dtype)
    extracted_array.mask = np.ma.masked
    extracted_array[small_slices] = array_large[large_slices]
    return extracted_array


def add_array(array_large, array_small, position):
    """
    Add a smaller array at a given position in a larger array.

    Parameters
    ----------
    array_large : `~numpy.ndarray`
        Large array.
    array_small : `~numpy.ndarray`
        Small array to add.
    position : tuple
        Position of the small array's center, with respect to the large array.
        Coordinates should be in the same order as the array shape.

    Returns
    -------
    new_array : `~numpy.ndarray`
        The new array formed from the sum of ``array_large`` and
        ``array_small``.

    Notes
    -----
    The addition is done in-place.

    Examples
    --------
    We consider a large array of zeros with the shape 5x5 and a small
    array of ones with a shape of 3x3:

    >>> import numpy as np
    >>> from astropy.nddata.utils import add_array
    >>> large_array = np.zeros((5, 5))
    >>> small_array = np.ones((3, 3))
    >>> add_array(large_array, small_array, (1, 2))
    array([[ 0.,  1.,  1.,  1.,  0.],
           [ 0.,  1.,  1.,  1.,  0.],
           [ 0.,  1.,  1.,  1.,  0.],
           [ 0.,  0.,  0.,  0.,  0.],
           [ 0.,  0.,  0.,  0.,  0.]])
    """
    # Check if large array is really larger
    if all(large_shape > small_shape for (large_shape, small_shape)
           in zip(array_large.shape, array_small.shape)):
        large_slices, small_slices = overlap_slices(array_large.shape,
                                                    array_small.shape, position)
        array_large[large_slices] += array_small[small_slices]
        return array_large
    else:
        raise ValueError("Can't add array. Small array too large.")


def subpixel_indices(position, subsampling):
    """
    Convert decimal points to indices, given a subsampling factor.

    This discards the integer part of the position and uses only the decimal
    place, and converts this to a subpixel position depending on the
    subsampling specified. The center of a pixel corresponds to an integer
    position.

    Parameters
    ----------
    position : `~numpy.ndarray` or array-like
        Positions in pixels.
    subsampling : int
        Subsampling factor per pixel.

    Returns
    -------
    indices : `~numpy.ndarray`
        The integer subpixel indices corresponding to the input positions.

    Examples
    --------

    If no subsampling is used, then the subpixel indices returned are always 0:

    >>> from astropy.nddata.utils import subpixel_indices
    >>> subpixel_indices([1.2, 3.4, 5.6],1)
    array([ 0.,  0.,  0.])

    If instead we use a subsampling of 2, we see that for the two first values
    (1.1 and 3.4) the subpixel position is 1, while for 5.6 it is 0. This is
    because the values of 1, 3, and 6 lie in the center of pixels, and 1.1 and
    3.4 lie in the left part of the pixels and 5.6 lies in the right part.

    >>> subpixel_indices([1.2, 3.4, 5.5],2)
    array([ 1.,  1.,  0.])
    """
    # Get decimal points
    fractions = np.modf(np.asanyarray(position) + 0.5)[0]
    return np.floor(fractions * subsampling)
