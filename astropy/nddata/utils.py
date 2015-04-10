# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module includes helper functions for array operations.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import numpy as np
from .decorators import support_nddata


__all__ = ['extract_array', 'add_array', 'subpixel_indices',
           'overlap_slices', 'block_reduce', 'block_replicate']


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
    position : tuple
        Position of the small array's center, with respect to the large array.
        Coordinates should be in the same order as the array shape.

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
    # Get edge coordinates
    edges_min = [int(pos + 0.5 - small_shape / 2.) for (pos, small_shape) in
                 zip(position, small_array_shape)]
    edges_max = [int(pos + 0.5 + small_shape / 2.) for (pos, small_shape) in
                 zip(position, small_array_shape)]

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
    array_small : `~numpy.ndarray`
        The extracted array

    Examples
    --------
    We consider a large array with the shape 11x10, from which we extract
    a small array of shape 3x5:

    >>> import numpy as np
    >>> from astropy.nddata.utils import extract_array
    >>> large_array = np.arange(110).reshape((11, 10))
    >>> large_array[4:9, 4:9] = np.ones((5, 5))
    >>> extract_array(large_array, (3, 5), (7, 7))
    array([[ 1,  1,  1,  1, 69],
           [ 1,  1,  1,  1, 79],
           [ 1,  1,  1,  1, 89]])
    """
    # Check if larger array is really larger
    if all(large_shape > small_shape for (large_shape, small_shape)
           in zip(array_large.shape, shape)):
        large_slices, _ = overlap_slices(array_large.shape, shape, position)
        return array_large[large_slices]
    else:
        raise ValueError("Can't extract array. Shape too large.")


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


@support_nddata
def block_reduce(data, block_size, func=np.sum):
    """
    Downsample a data array by applying a function to local blocks.

    If ``data`` is not perfectly divisible by ``block_size`` along a
    given axis then the data will be trimmed (from the end) along that
    axis.

    Parameters
    ----------
    data : array_like
        The data to be resampled.

    block_size : int or array_like (int)
        The integer block size along each axis.  If ``block_size`` is a
        scalar and ``data`` has more than one dimension, then
        ``block_size`` will be used for for every axis.

    func : callable
        The method to use to downsample the data.  Must be a callable
        that takes in a `~numpy.ndarray` along with an ``axis`` keyword,
        which defines the axis along which the function is applied.  The
        default is `~numpy.sum`, which provides block summation (and
        conserves the data sum).

    Returns
    -------
    output : array-like
        The resampled data.

    Examples
    --------
    >>> import numpy as np
    >>> from astropy.nddata.utils import block_reduce
    >>> data = np.arange(16).reshape(4, 4)
    >>> block_reduce(data, 2)
    array([[10, 18],
           [42, 50]])

    >>> block_reduce(data, 2, func=np.mean)
    array([[  2.5,   4.5],
           [ 10.5,  12.5]])
    """

    from skimage.measure import block_reduce

    data = np.asanyarray(data)

    block_size = np.atleast_1d(block_size)
    if data.ndim > 1 and len(block_size) == 1:
        block_size = np.repeat(block_size, data.ndim)

    if len(block_size) != data.ndim:
        raise ValueError('`block_size` must be a scalar or have the same '
                         'length as `data.shape`')

    block_size = np.array([int(i) for i in block_size])
    size_resampled = np.array(data.shape) // block_size
    size_init = size_resampled * block_size

    # trim data if necessary
    for i in range(data.ndim):
        if data.shape[i] != size_init[i]:
            data = data.swapaxes(0, i)
            data = data[:size_init[i]]
            data = data.swapaxes(0, i)

    return block_reduce(data, tuple(block_size), func=func)


@support_nddata
def block_replicate(data, block_size, conserve_sum=True):
    """
    Upsample a data array by block replication.

    Parameters
    ----------
    data : array_like
        The data to be block replicated.

    block_size : int or array_like (int)
        The integer block size along each axis.  If ``block_size`` is a
        scalar and ``data`` has more than one dimension, then
        ``block_size`` will be used for for every axis.

    conserve_sum : bool
        If `True` (the default) then the sum of the output
        block-replicated data will equal the sum of the input ``data``.

    Returns
    -------
    output : array_like
        The block-replicated data.

    Examples
    --------
    >>> import numpy as np
    >>> from astropy.nddata.utils import block_replicate
    >>> data = np.array([[0., 1.], [2., 3.]])
    >>> block_replicate(data, 2)
    array([[ 0.  ,  0.  ,  0.25,  0.25],
           [ 0.  ,  0.  ,  0.25,  0.25],
           [ 0.5 ,  0.5 ,  0.75,  0.75],
           [ 0.5 ,  0.5 ,  0.75,  0.75]])

    >>> block_replicate(data, 2, conserve_sum=False)
    array([[ 0.,  0.,  1.,  1.],
           [ 0.,  0.,  1.,  1.],
           [ 2.,  2.,  3.,  3.],
           [ 2.,  2.,  3.,  3.]])
    """

    data = np.asanyarray(data)

    block_size = np.atleast_1d(block_size)
    if data.ndim > 1 and len(block_size) == 1:
        block_size = np.repeat(block_size, data.ndim)

    if len(block_size) != data.ndim:
        raise ValueError('`block_size` must be a scalar or have the same '
                         'length as `data.shape`')

    for i in range(data.ndim):
        data = np.repeat(data, block_size[i], axis=i)

    if conserve_sum:
        data = data / float(np.prod(block_size))

    return data
