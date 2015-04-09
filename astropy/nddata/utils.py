# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module includes helper functions for array operations.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np

__all__ = ['extract_array', 'add_array', 'subpixel_indices',
           'overlap_slices', 'NoOverlapError', 'PartialOverlapError']


class NoOverlapError(ValueError):
    '''Raised when determining the overlap of non-overlapping arrays.'''
    pass


class PartialOverlapError(ValueError):
    '''Raised when arrays only partially overlap.'''
    pass


def _round(a):
    '''Always round up.

    ``np.round`` cannot be used here, because it round .5 to the nearest
    even number.
    '''
    return int(np.floor(a + 0.5))


def _offset(a):
    '''Offset by 0.5 for an even array.

    For an array with an odd number of elements, the center is
    symmetric, e.g. for 3 elements, it's center +/-1 elements, but for
    four elements it's center -2 / +1
    This function introduces that offset.
    '''
    if np.mod(a, 2) == 0:
        return -0.5
    else:
        return 0.


def overlap_slices(large_array_shape, small_array_shape, position, mode='partial'):
    """
    Get slices for the overlapping part of a small and a large array.

    Given a certain position of the center of the small array, with
    respect to the large array, tuples of slices are returned which can be
    used to extract, add or subtract the small array at the given
    position. This function takes care of the correct behavior at the
    boundaries, where the small array is cut of appropriately.
    Integer positions are at the pixel centers.

    Parameters
    ----------
    large_array_shape : tuple or number
        Shape of the large array (for 1D arrays, this can be an int or float).
    small_array_shape : tuple or number
        Shape of the small array (for 1D arrays, this can be an int or float).
    position : tuple of integers or int
        Position of the small array's center, with respect to the large array.
        Coordinates should be in the same order as the array shape.
        For a coordinate with an even number of elements, the position is
        rounded up, e.g. extracting two elements with a center of ``1`` will
        give positions ``[0, 1]``.
    mode : ['partial', 'strict']
        In "partial" mode, a partial overlap of the small and the large
        array is sufficient. In the "strict" mode, the small array has to be
        fully contained in the large array, otherwise an
        `~astropy.nddata.utils.PartialOverlapError` is raised. In both modes,
        non-overlapping arrays will raise a `~astropy.nddata.utils.NoOverlapError`.

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
    if mode not in ['partial', 'strict']:
        raise ValueError('Mode can only be "partial" or "strict".')
    if np.isscalar(small_array_shape):
        small_array_shape = (small_array_shape, )
    if np.isscalar(large_array_shape):
        large_array_shape = (large_array_shape, )
    if np.isscalar(position):
        position = (position, )

    if len(small_array_shape) != len(large_array_shape):
        raise ValueError("Both arrays must have the same number of dimensions.")

    if len(small_array_shape) != len(position):
        raise ValueError("Position must have the same number of dimensions as array.")

    # Get edge coordinates
    edges_min = [_round(pos + 0.5 - small_shape / 2. + _offset(small_shape))
                 for (pos, small_shape) in zip(position, small_array_shape)]
    edges_max = [_round(pos + 0.5 + small_shape / 2. + _offset(small_shape))
                 for (pos, small_shape) in zip(position, small_array_shape)]

    for e_max in edges_max:
        if e_max <= 0:
            raise NoOverlapError('Arrays do not overlap.')
    for e_min, large_shape in zip(edges_min, large_array_shape):
        if e_min >= large_shape:
            raise NoOverlapError('Arrays do not overlap.')

    if mode == 'strict':
        for e_min in edges_min:
            if e_min < 0:
                raise PartialOverlapError('Arrays overlap only partially.')
        for e_max, large_shape in zip(edges_max, large_array_shape):
            if e_max >= large_shape:
                raise PartialOverlapError('Arrays overlap only partially.')

    # Set up slices
    slices_large = tuple(slice(max(0, edge_min), min(large_shape, edge_max))
                         for (edge_min, edge_max, large_shape) in
                         zip(edges_min, edges_max, large_array_shape))
    slices_small = tuple(slice(max(0, -edge_min),
                               min(large_shape - edge_min, edge_max - edge_min))
                         for (edge_min, edge_max, large_shape) in
                         zip(edges_min, edges_max, large_array_shape))

    return slices_large, slices_small


def extract_array(array_large, shape, position, mode='partial', fill_value=np.nan,
                  return_position=False):
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
    mode : ['partial', 'trim', 'strict']
        In "partial" and "trim" mode, a partial overlap of the small and the large
        array is sufficient. In the "strict" mode, the small array has to be
        fully contained in the large array, otherwise an
        `~astropy.nddata.utils.PartialOverlapError` is raised. In all modes,
        non-overlapping arrays will raise a `~astropy.nddata.utils.NoOverlapError`.
        In "partial" mode, positions in the extracted array, that do not overlap
        with the original array, will be filled with ``fill_value``. In "trim" mode
        only the overlapping elements are returned, thus the resulting array may
        be smaller than requested.

    fill_value : object of type array_large.dtype
        In "partial" mode ``fill_value`` set the values in the extracted array that
        do not overlap with ``large_array``.

    return_position : boolean
        If true, return the coordinates of ``position`` in the coordinate system
        of the returned array.

    Returns
    -------
    array_small : `~numpy.ndarray`
        The extracted array.

    new_position : tuple
        If ``return_position`` is true, this tuple will contain the coordinates
        of the input ``position`` in the coordinate system of ``array_small``. Note 
        that for partially overlapping arrays, ``new_position`` might actually be 
        outside of the ``array_small``; ``array_small[new_position]`` might give 
        wrong results if any element in ``new_position`` is negative.

    Examples
    --------
    We consider a large array with the shape 11x10, from which we extract
    a small array of shape 3x5:

    >>> import numpy as np
    >>> from astropy.nddata.utils import extract_array
    >>> large_array = np.arange(110).reshape((11, 10))
    >>> extract_array(large_array, (3, 5), (7, 7))
    array([[65, 66, 67, 68, 69],
           [75, 76, 77, 78, 79],
           [85, 86, 87, 88, 89]])
    """
    if np.isscalar(shape):
        shape = (shape, )
    if np.isscalar(position):
        position = (position, )

    if mode in ['partial', 'trim']:
        slicemode = 'partial'
    elif mode == 'strict':
        slicemode = mode
    else:
        raise ValueError("Valid modes are 'partial', 'trim', and 'strict'.")
    large_slices, small_slices = overlap_slices(array_large.shape,
                                                shape, position, mode=slicemode)
    extracted_array = array_large[large_slices]
    if return_position:
        new_position = [i - s.start for i, s in zip(position, large_slices)]
    # Extracting on the edges is presumably a rare case, so treat special here.
    if (extracted_array.shape != shape) and (mode == 'partial'):
        extracted_array = np.zeros(shape, dtype=array_large.dtype)
        extracted_array[:] = fill_value
        extracted_array[small_slices] = array_large[large_slices]
        if return_position:
            new_position = [i + s.start for i, s in zip(new_position, small_slices)]
    if return_position:
        return extracted_array, tuple(new_position)
    else:
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
