# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module includes helper functions for array operations.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
from numpy.lib.index_tricks import index_exp

__all__ = ['extract_array', 'add_array', 'subpixel_indices',
           'mask_to_mirrored_num', 'overlap_slices']


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
    array_large : ndarray
        Array to extract another array from.
    shape : tuple
        Shape of the extracted array.
    position : tuple
        Position of the small array's center, with respect to the large array.
        Coordinates should be in the same order as the array shape.

    Examples
    --------
    We consider a large array with the shape 11x10, from which we extract
    a small array of shape 3x5:

    >>> import numpy as np
    >>> from astropy.image.array_utils import extract_array
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
    array_large : ndarray
        Large array.
    array_small : ndarray
        Small array to add.
    position : tuple
        Position of the small array's center, with respect to the large array.
        Coordinates should be in the same order as the array shape.

    Notes
    -----
    The addition is done in-place.

    Examples
    --------
    We consider a large array of zeros with the shape 5x5 and a small
    array of ones with a shape of 3x3:

    >>> import numpy as np
    >>> from astropy.image.array_utils import add_array
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

    Parameters
    ----------
    position : ndarray or array-like
        Positions in pixels.
    subsampling : int
        Subsampling factor per pixel.
    """
    # Get decimal points
    fractions = np.modf(np.asanyarray(position) + 0.5)[0]
    return np.floor(fractions * subsampling)


def mask_to_mirrored_num(image, mask_image, center_position, bbox=None,
                         copy=True):
    """
    Replace masked pixels with the value of the pixel mirrored across a
    given ``center_position``.  If the mirror pixel is unavailable (i.e.
    itself masked or outside of the image), then the masked pixel value
    is set to zero.

    Parameters
    ----------
    image : `numpy.ndarray`, 2D
        The 2D array of the image.

    mask_image : array-like, bool
        A boolean mask with the same shape as ``image``, where a `True`
        value indicates the corresponding element of ``image`` is
        considered bad.

    center_position : 2-tuple
        (x, y) center coordinates around which masked pixels will be
        mirrored.

    bbox : list, tuple, `numpy.ndarray`, optional
        The bounding box (x_min, x_max, y_min, y_max) over which to
        replace masked pixels.

    copy : bool, optional
        Whether to return a new copy of the array, fixed (default), or whether
        to fix the original array in-place.

    Returns
    -------
    result : `numpy.ndarray`, 2D
        A 2D array with replaced masked pixels.

    Examples
    --------
    >>> import numpy as np
    >>> from astropy.image import mask_to_mirrored_num
    >>> image = np.arange(16).reshape(4, 4)
    >>> mask = np.zeros_like(image, dtype=bool)
    >>> mask[0, 0] = True
    >>> mask[1, 1] = True
    >>> mask_to_mirrored_num(image, mask, (1.5, 1.5))
    array([[15,  1,  2,  3],
           [ 4, 10,  6,  7],
           [ 8,  9, 10, 11],
           [12, 13, 14, 15]])
    """

    if bbox is None:
        ny, nx = image.shape
        bbox = [0, nx, 0, ny]

    region = index_exp[bbox[2]:bbox[3]+1, bbox[0]:bbox[1]+1]

    if copy:
        subdata = image[region].copy()
    else:
        subdata = image[region]

    submask = mask_image[region]
    y_masked, x_masked = np.nonzero(submask)
    x_mirror = (2 * (center_position[0] - bbox[0])
                - x_masked + 0.5).astype('int32')
    y_mirror = (2 * (center_position[1] - bbox[2])
                - y_masked + 0.5).astype('int32')

    # Reset mirrored pixels that go out of the image.
    outofimage = ((x_mirror < 0) | (y_mirror < 0) |
                  (x_mirror >= subdata.shape[1]) |
                  (y_mirror >= subdata.shape[0]))
    if outofimage.any():
        x_mirror[outofimage] = x_masked[outofimage].astype('int32')
        y_mirror[outofimage] = y_masked[outofimage].astype('int32')

    subdata[y_masked, x_masked] = subdata[y_mirror, x_mirror]

    # Set pixels that mirrored to another masked pixel to zero.
    # This will also set to zero any pixels that mirrored out of
    # the image.
    mirror_is_masked = submask[y_mirror, x_mirror]
    x_bad = x_masked[mirror_is_masked]
    y_bad = y_masked[mirror_is_masked]
    subdata[y_bad, x_bad] = 0.0

    if copy:
        outimage = image.copy()
    else:
        outimage = image

    outimage[region] = subdata
    return outimage
