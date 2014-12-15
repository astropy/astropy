# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module includes helper functions for array operations.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import numpy as np
import copy

__all__ = ['extract_array_2d', 'add_array_2d', 'subpixel_indices',
           'mask_to_mirrored_num', 'overlap_slices']


def overlap_slices(large_array_shape, small_array_shape, position):
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
    y_min = int(position[1] + 0.5 - small_array_shape[0] / 2.)
    x_min = int(position[0] + 0.5 - small_array_shape[1] / 2.)
    y_max = int(position[1] + 0.5 + small_array_shape[0] / 2.)
    x_max = int(position[0] + 0.5 + small_array_shape[1] / 2.)

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
    We consider a large array with the shape 11x10, from which we extract
    a small array of shape 3x5:

    >>> import numpy as np
    >>> from astropy.image.array_utils import extract_array_2d
    >>> large_array = np.arange(110).reshape((11, 10))
    >>> large_array[4:9, 4:9] = np.ones((5, 5))
    >>> extract_array_2d(large_array, (3, 5), (7, 7))
    array([[ 1,  1,  1,  1, 69],
           [ 1,  1,  1,  1, 79],
           [ 1,  1,  1,  1, 89]])
    """
    # Check if larger array is really larger
    if array_large.shape >= shape:
        s_y, s_x, _, _ = overlap_slices(array_large.shape, shape, position)
        return array_large[s_y, s_x]
    else:
        raise ValueError("Can't extract array. Shape too large.")


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
    >>> from astropy.image.array_utils import add_array_2d
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
        s_y, s_x, b_y, b_x = overlap_slices(array_large.shape,
                                         array_small.shape, position)
        array_large[s_y, s_x] += array_small[b_y, b_x]
        return array_large
    else:
        raise ValueError("Can't add array. Small array too large.")


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


def mask_to_mirrored_num(image, mask_image, center_position, bbox=None):
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
    subdata = copy.deepcopy(image[bbox[2]:bbox[3]+1, bbox[0]:bbox[1]+1])
    submask = mask_image[bbox[2]:bbox[3]+1, bbox[0]:bbox[1]+1]
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

    outimage = copy.deepcopy(image)
    outimage[bbox[2]:bbox[3]+1, bbox[0]:bbox[1]+1] = subdata
    return outimage
