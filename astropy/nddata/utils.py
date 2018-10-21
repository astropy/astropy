# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module includes helper functions for array operations.
"""
import os
from copy import deepcopy
from functools import partial
import math
import errno

import numpy as np

from .decorators import support_nddata
from .. import units as u
from ..coordinates import SkyCoord
from ..utils import lazyproperty
from ..wcs import Sip, WCS, NoConvergence
from ..wcs.utils import skycoord_to_pixel, proj_plane_pixel_scales
from ..io import fits
from ..io.fits import PrimaryHDU, ImageHDU, CompImageHDU
from ..table import Table, QTable
from .nddata import NDData
from .. import log


__all__ = ['extract_array', 'add_array', 'subpixel_indices',
           'overlap_slices', 'block_reduce', 'block_replicate',
           'NoOverlapError', 'PartialOverlapError', 'Cutout2D']


class NoOverlapError(ValueError):
    '''Raised when determining the overlap of non-overlapping arrays.'''
    pass


class PartialOverlapError(ValueError):
    '''Raised when arrays only partially overlap.'''
    pass


def _round(a):
    '''Always round up.

    ``np.round`` cannot be used here, because it rounds .5 to the nearest
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


def overlap_slices(large_array_shape, small_array_shape, position,
                   mode='partial'):
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
    large_array_shape : tuple of int or int
        The shape of the large array (for 1D arrays, this can be an
        `int`).
    small_array_shape : tuple of int or int
        The shape of the small array (for 1D arrays, this can be an
        `int`).  See the ``mode`` keyword for additional details.
    position : tuple of numbers or number
        The position of the small array's center with respect to the
        large array.  The pixel coordinates should be in the same order
        as the array shape.  Integer positions are at the pixel centers.
        For any axis where ``small_array_shape`` is even, the position
        is rounded up, e.g. extracting two elements with a center of
        ``1`` will define the extracted region as ``[0, 1]``.
    mode : {'partial', 'trim', 'strict'}, optional
        In ``'partial'`` mode, a partial overlap of the small and the
        large array is sufficient.  The ``'trim'`` mode is similar to
        the ``'partial'`` mode, but ``slices_small`` will be adjusted to
        return only the overlapping elements.  In the ``'strict'`` mode,
        the small array has to be fully contained in the large array,
        otherwise an `~astropy.nddata.utils.PartialOverlapError` is
        raised.  In all modes, non-overlapping arrays will raise a
        `~astropy.nddata.utils.NoOverlapError`.

    Returns
    -------
    slices_large : tuple of slices
        A tuple of slice objects for each axis of the large array, such
        that ``large_array[slices_large]`` extracts the region of the
        large array that overlaps with the small array.
    slices_small : tuple of slices
        A tuple of slice objects for each axis of the small array, such
        that ``small_array[slices_small]`` extracts the region that is
        inside the large array.
    """

    if mode not in ['partial', 'trim', 'strict']:
        raise ValueError('Mode can be only "partial", "trim", or "strict".')
    if np.isscalar(small_array_shape):
        small_array_shape = (small_array_shape, )
    if np.isscalar(large_array_shape):
        large_array_shape = (large_array_shape, )
    if np.isscalar(position):
        position = (position, )

    if len(small_array_shape) != len(large_array_shape):
        raise ValueError('"large_array_shape" and "small_array_shape" must '
                         'have the same number of dimensions.')

    if len(small_array_shape) != len(position):
        raise ValueError('"position" must have the same number of dimensions '
                         'as "small_array_shape".')
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
    if mode == 'trim':
        slices_small = tuple(slice(0, slc.stop - slc.start)
                             for slc in slices_large)
    else:
        slices_small = tuple(slice(max(0, -edge_min),
                                   min(large_shape - edge_min,
                                       edge_max - edge_min))
                             for (edge_min, edge_max, large_shape) in
                             zip(edges_min, edges_max, large_array_shape))

    return slices_large, slices_small


def extract_array(array_large, shape, position, mode='partial',
                  fill_value=np.nan, return_position=False):
    """
    Extract a smaller array of the given shape and position from a
    larger array.

    Parameters
    ----------
    array_large : `~numpy.ndarray`
        The array from which to extract the small array.
    shape : tuple or int
        The shape of the extracted array (for 1D arrays, this can be an
        `int`).  See the ``mode`` keyword for additional details.
    position : tuple of numbers or number
        The position of the small array's center with respect to the
        large array.  The pixel coordinates should be in the same order
        as the array shape.  Integer positions are at the pixel centers
        (for 1D arrays, this can be a number).
    mode : {'partial', 'trim', 'strict'}, optional
        The mode used for extracting the small array.  For the
        ``'partial'`` and ``'trim'`` modes, a partial overlap of the
        small array and the large array is sufficient.  For the
        ``'strict'`` mode, the small array has to be fully contained
        within the large array, otherwise an
        `~astropy.nddata.utils.PartialOverlapError` is raised.   In all
        modes, non-overlapping arrays will raise a
        `~astropy.nddata.utils.NoOverlapError`.  In ``'partial'`` mode,
        positions in the small array that do not overlap with the large
        array will be filled with ``fill_value``.  In ``'trim'`` mode
        only the overlapping elements are returned, thus the resulting
        small array may be smaller than the requested ``shape``.
    fill_value : number, optional
        If ``mode='partial'``, the value to fill pixels in the extracted
        small array that do not overlap with the input ``array_large``.
        ``fill_value`` must have the same ``dtype`` as the
        ``array_large`` array.
    return_position : boolean, optional
        If `True`, return the coordinates of ``position`` in the
        coordinate system of the returned array.

    Returns
    -------
    array_small : `~numpy.ndarray`
        The extracted array.
    new_position : tuple
        If ``return_position`` is true, this tuple will contain the
        coordinates of the input ``position`` in the coordinate system
        of ``array_small``. Note that for partially overlapping arrays,
        ``new_position`` might actually be outside of the
        ``array_small``; ``array_small[new_position]`` might give wrong
        results if any element in ``new_position`` is negative.

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

    if mode not in ['partial', 'trim', 'strict']:
        raise ValueError("Valid modes are 'partial', 'trim', and 'strict'.")
    large_slices, small_slices = overlap_slices(array_large.shape,
                                                shape, position, mode=mode)
    extracted_array = array_large[large_slices]
    if return_position:
        new_position = [i - s.start for i, s in zip(position, large_slices)]
    # Extracting on the edges is presumably a rare case, so treat special here
    if (extracted_array.shape != shape) and (mode == 'partial'):
        extracted_array = np.zeros(shape, dtype=array_large.dtype)
        extracted_array[:] = fill_value
        extracted_array[small_slices] = array_large[large_slices]
        if return_position:
            new_position = [i + s.start for i, s in zip(new_position,
                                                        small_slices)]
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
    >>> add_array(large_array, small_array, (1, 2))  # doctest: +FLOAT_CMP
    array([[0., 1., 1., 1., 0.],
           [0., 1., 1., 1., 0.],
           [0., 1., 1., 1., 0.],
           [0., 0., 0., 0., 0.],
           [0., 0., 0., 0., 0.]])
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
    >>> subpixel_indices([1.2, 3.4, 5.6], 1)  # doctest: +FLOAT_CMP
    array([0., 0., 0.])

    If instead we use a subsampling of 2, we see that for the two first values
    (1.1 and 3.4) the subpixel position is 1, while for 5.6 it is 0. This is
    because the values of 1, 3, and 6 lie in the center of pixels, and 1.1 and
    3.4 lie in the left part of the pixels and 5.6 lies in the right part.

    >>> subpixel_indices([1.2, 3.4, 5.5], 2)  # doctest: +FLOAT_CMP
    array([1., 1., 0.])
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

    func : callable, optional
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
    >>> block_reduce(data, 2)    # doctest: +SKIP
    array([[10, 18],
           [42, 50]])

    >>> block_reduce(data, 2, func=np.mean)    # doctest: +SKIP
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

    conserve_sum : bool, optional
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
    >>> block_replicate(data, 2)  # doctest: +FLOAT_CMP
    array([[0.  , 0.  , 0.25, 0.25],
           [0.  , 0.  , 0.25, 0.25],
           [0.5 , 0.5 , 0.75, 0.75],
           [0.5 , 0.5 , 0.75, 0.75]])

    >>> block_replicate(data, 2, conserve_sum=False)  # doctest: +FLOAT_CMP
    array([[0., 0., 1., 1.],
           [0., 0., 1., 1.],
           [2., 2., 3., 3.],
           [2., 2., 3., 3.]])
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


class Cutout2D:
    """
    Create a cutout object from a 2D array.

    The returned object will contain a 2D cutout array.  If
    ``copy=False`` (default), the cutout array is a view into the
    original ``data`` array, otherwise the cutout array will contain a
    copy of the original data.

    If a `~astropy.wcs.WCS` object is input, then the returned object
    will also contain a copy of the original WCS, but updated for the
    cutout array.

    For example usage, see :ref:`cutout_images`.

    .. warning::

        The cutout WCS object does not currently handle cases where the
        input WCS object contains distortion lookup tables described in
        the `FITS WCS distortion paper
        <http://www.atnf.csiro.au/people/mcalabre/WCS/dcs_20040422.pdf>`__.

    Parameters
    ----------
    data : `~numpy.ndarray`
        The 2D data array from which to extract the cutout array.

    position : tuple or `~astropy.coordinates.SkyCoord`
        The position of the cutout array's center with respect to
        the ``data`` array.  The position can be specified either as
        a ``(x, y)`` tuple of pixel coordinates or a
        `~astropy.coordinates.SkyCoord`, in which case ``wcs`` is a
        required input.

    size : int, array-like, `~astropy.units.Quantity`
        The size of the cutout array along each axis.  If ``size``
        is a scalar number or a scalar `~astropy.units.Quantity`,
        then a square cutout of ``size`` will be created.  If
        ``size`` has two elements, they should be in ``(ny, nx)``
        order.  Scalar numbers in ``size`` are assumed to be in
        units of pixels.  ``size`` can also be a
        `~astropy.units.Quantity` object or contain
        `~astropy.units.Quantity` objects.  Such
        `~astropy.units.Quantity` objects must be in pixel or
        angular units.  For all cases, ``size`` will be converted to
        an integer number of pixels, rounding the the nearest
        integer.  See the ``mode`` keyword for additional details on
        the final cutout size.

        .. note::
            If ``size`` is in angular units, the cutout size is
            converted to pixels using the pixel scales along each
            axis of the image at the ``CRPIX`` location.  Projection
            and other non-linear distortions are not taken into
            account.

    wcs : `~astropy.wcs.WCS`, optional
        A WCS object associated with the input ``data`` array.  If
        ``wcs`` is not `None`, then the returned cutout object will
        contain a copy of the updated WCS for the cutout data array.

    mode : {'trim', 'partial', 'strict'}, optional
        The mode used for creating the cutout data array.  For the
        ``'partial'`` and ``'trim'`` modes, a partial overlap of the
        cutout array and the input ``data`` array is sufficient.
        For the ``'strict'`` mode, the cutout array has to be fully
        contained within the ``data`` array, otherwise an
        `~astropy.nddata.utils.PartialOverlapError` is raised.   In
        all modes, non-overlapping arrays will raise a
        `~astropy.nddata.utils.NoOverlapError`.  In ``'partial'``
        mode, positions in the cutout array that do not overlap with
        the ``data`` array will be filled with ``fill_value``.  In
        ``'trim'`` mode only the overlapping elements are returned,
        thus the resulting cutout array may be smaller than the
        requested ``shape``.

    fill_value : number, optional
        If ``mode='partial'``, the value to fill pixels in the
        cutout array that do not overlap with the input ``data``.
        ``fill_value`` must have the same ``dtype`` as the input
        ``data`` array.

    copy : bool, optional
        If `False` (default), then the cutout data will be a view
        into the original ``data`` array.  If `True`, then the
        cutout data will hold a copy of the original ``data`` array.

    Attributes
    ----------
    data : 2D `~numpy.ndarray`
        The 2D cutout array.

    shape : 2 tuple
        The ``(ny, nx)`` shape of the cutout array.

    shape_input : 2 tuple
        The ``(ny, nx)`` shape of the input (original) array.

    input_position_cutout : 2 tuple
        The (unrounded) ``(x, y)`` position with respect to the cutout
        array.

    input_position_original : 2 tuple
        The original (unrounded) ``(x, y)`` input position (with respect
        to the original array).

    slices_original : 2 tuple of slice objects
        A tuple of slice objects for the minimal bounding box of the
        cutout with respect to the original array.  For
        ``mode='partial'``, the slices are for the valid (non-filled)
        cutout values.

    slices_cutout : 2 tuple of slice objects
        A tuple of slice objects for the minimal bounding box of the
        cutout with respect to the cutout array.  For
        ``mode='partial'``, the slices are for the valid (non-filled)
        cutout values.

    xmin_original, ymin_original, xmax_original, ymax_original : float
        The minimum and maximum ``x`` and ``y`` indices of the minimal
        rectangular region of the cutout array with respect to the
        original array.  For ``mode='partial'``, the bounding box
        indices are for the valid (non-filled) cutout values.  These
        values are the same as those in `bbox_original`.

    xmin_cutout, ymin_cutout, xmax_cutout, ymax_cutout : float
        The minimum and maximum ``x`` and ``y`` indices of the minimal
        rectangular region of the cutout array with respect to the
        cutout array.  For ``mode='partial'``, the bounding box indices
        are for the valid (non-filled) cutout values.  These values are
        the same as those in `bbox_cutout`.

    wcs : `~astropy.wcs.WCS` or `None`
        A WCS object associated with the cutout array if a ``wcs``
        was input.

    Examples
    --------
    >>> import numpy as np
    >>> from astropy.nddata.utils import Cutout2D
    >>> from astropy import units as u
    >>> data = np.arange(20.).reshape(5, 4)
    >>> cutout1 = Cutout2D(data, (2, 2), (3, 3))
    >>> print(cutout1.data)  # doctest: +FLOAT_CMP
    [[ 5.  6.  7.]
     [ 9. 10. 11.]
     [13. 14. 15.]]

    >>> print(cutout1.center_original)
    (2.0, 2.0)
    >>> print(cutout1.center_cutout)
    (1.0, 1.0)
    >>> print(cutout1.origin_original)
    (1, 1)

    >>> cutout2 = Cutout2D(data, (2, 2), 3)
    >>> print(cutout2.data)  # doctest: +FLOAT_CMP
    [[ 5.  6.  7.]
     [ 9. 10. 11.]
     [13. 14. 15.]]

    >>> size = u.Quantity([3, 3], u.pixel)
    >>> cutout3 = Cutout2D(data, (0, 0), size)
    >>> print(cutout3.data)  # doctest: +FLOAT_CMP
    [[0. 1.]
     [4. 5.]]

    >>> cutout4 = Cutout2D(data, (0, 0), (3 * u.pixel, 3))
    >>> print(cutout4.data)  # doctest: +FLOAT_CMP
    [[0. 1.]
     [4. 5.]]

    >>> cutout5 = Cutout2D(data, (0, 0), (3, 3), mode='partial')
    >>> print(cutout5.data)  # doctest: +FLOAT_CMP
    [[nan nan nan]
     [nan  0.  1.]
     [nan  4.  5.]]
    """

    def __init__(self, data, position, size, wcs=None, mode='trim',
                 fill_value=np.nan, copy=False):
        if isinstance(position, SkyCoord):
            if wcs is None:
                raise ValueError('wcs must be input if position is a '
                                 'SkyCoord')
            position = skycoord_to_pixel(position, wcs, mode='all')  # (x, y)

        if np.isscalar(size):
            size = np.repeat(size, 2)

        # special handling for a scalar Quantity
        if isinstance(size, u.Quantity):
            size = np.atleast_1d(size)
            if len(size) == 1:
                size = np.repeat(size, 2)

        if len(size) > 2:
            raise ValueError('size must have at most two elements')

        shape = np.zeros(2).astype(int)
        pixel_scales = None
        # ``size`` can have a mixture of int and Quantity (and even units),
        # so evaluate each axis separately
        for axis, side in enumerate(size):
            if not isinstance(side, u.Quantity):
                shape[axis] = int(np.round(size[axis]))     # pixels
            else:
                if side.unit == u.pixel:
                    shape[axis] = int(np.round(side.value))
                elif side.unit.physical_type == 'angle':
                    if wcs is None:
                        raise ValueError('wcs must be input if any element '
                                         'of size has angular units')
                    if pixel_scales is None:
                        pixel_scales = u.Quantity(
                            proj_plane_pixel_scales(wcs), wcs.wcs.cunit[axis])
                    shape[axis] = int(np.round(
                        (side / pixel_scales[axis]).decompose()))
                else:
                    raise ValueError('shape can contain Quantities with only '
                                     'pixel or angular units')

        data = np.asanyarray(data)
        # reverse position because extract_array and overlap_slices
        # use (y, x), but keep the input position
        pos_yx = position[::-1]

        cutout_data, input_position_cutout = extract_array(
            data, tuple(shape), pos_yx, mode=mode, fill_value=fill_value,
            return_position=True)
        if copy:
            cutout_data = np.copy(cutout_data)
        self.data = cutout_data

        self.input_position_cutout = input_position_cutout[::-1]    # (x, y)
        slices_original, slices_cutout = overlap_slices(
            data.shape, shape, pos_yx, mode=mode)

        self.slices_original = slices_original
        self.slices_cutout = slices_cutout

        self.shape = self.data.shape
        self.input_position_original = position
        self.shape_input = shape

        ((self.ymin_original, self.ymax_original),
         (self.xmin_original, self.xmax_original)) = self.bbox_original

        ((self.ymin_cutout, self.ymax_cutout),
         (self.xmin_cutout, self.xmax_cutout)) = self.bbox_cutout

        # the true origin pixel of the cutout array, including any
        # filled cutout values
        self._origin_original_true = (
            self.origin_original[0] - self.slices_cutout[1].start,
            self.origin_original[1] - self.slices_cutout[0].start)

        if wcs is not None:
            self.wcs = deepcopy(wcs)
            self.wcs.wcs.crpix -= self._origin_original_true
            self.wcs._naxis = [self.data.shape[1], self.data.shape[0]]
            if wcs.sip is not None:
                self.wcs.sip = Sip(wcs.sip.a, wcs.sip.b,
                                   wcs.sip.ap, wcs.sip.bp,
                                   wcs.sip.crpix - self._origin_original_true)
        else:
            self.wcs = None

    def to_original_position(self, cutout_position):
        """
        Convert an ``(x, y)`` position in the cutout array to the original
        ``(x, y)`` position in the original large array.

        Parameters
        ----------
        cutout_position : tuple
            The ``(x, y)`` pixel position in the cutout array.

        Returns
        -------
        original_position : tuple
            The corresponding ``(x, y)`` pixel position in the original
            large array.
        """
        return tuple(cutout_position[i] + self.origin_original[i]
                     for i in [0, 1])

    def to_cutout_position(self, original_position):
        """
        Convert an ``(x, y)`` position in the original large array to
        the ``(x, y)`` position in the cutout array.

        Parameters
        ----------
        original_position : tuple
            The ``(x, y)`` pixel position in the original large array.

        Returns
        -------
        cutout_position : tuple
            The corresponding ``(x, y)`` pixel position in the cutout
            array.
        """
        return tuple(original_position[i] - self.origin_original[i]
                     for i in [0, 1])

    def plot_on_original(self, ax=None, fill=False, **kwargs):
        """
        Plot the cutout region on a matplotlib Axes instance.

        Parameters
        ----------
        ax : `matplotlib.axes.Axes` instance, optional
            If `None`, then the current `matplotlib.axes.Axes` instance
            is used.

        fill : bool, optional
            Set whether to fill the cutout patch.  The default is
            `False`.

        kwargs : optional
            Any keyword arguments accepted by `matplotlib.patches.Patch`.

        Returns
        -------
        ax : `matplotlib.axes.Axes` instance
            The matplotlib Axes instance constructed in the method if
            ``ax=None``.  Otherwise the output ``ax`` is the same as the
            input ``ax``.
        """

        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches

        kwargs['fill'] = fill

        if ax is None:
            ax = plt.gca()

        height, width = self.shape
        hw, hh = width / 2., height / 2.
        pos_xy = self.position_original - np.array([hw, hh])
        patch = mpatches.Rectangle(pos_xy, width, height, 0., **kwargs)
        ax.add_patch(patch)
        return ax

    @staticmethod
    def _calc_center(slices):
        """
        Calculate the center position.  The center position will be
        fractional for even-sized arrays.  For ``mode='partial'``, the
        central position is calculated for the valid (non-filled) cutout
        values.
        """
        return tuple(0.5 * (slices[i].start + slices[i].stop - 1)
                     for i in [1, 0])

    @staticmethod
    def _calc_bbox(slices):
        """
        Calculate a minimal bounding box in the form ``((ymin, ymax),
        (xmin, xmax))``.  Note these are pixel locations, not slice
        indices.  For ``mode='partial'``, the bounding box indices are
        for the valid (non-filled) cutout values.
        """
        # (stop - 1) to return the max pixel location, not the slice index
        return ((slices[0].start, slices[0].stop - 1),
                (slices[1].start, slices[1].stop - 1))

    @lazyproperty
    def origin_original(self):
        """
        The ``(x, y)`` index of the origin pixel of the cutout with
        respect to the original array.  For ``mode='partial'``, the
        origin pixel is calculated for the valid (non-filled) cutout
        values.
        """
        return (self.slices_original[1].start, self.slices_original[0].start)

    @lazyproperty
    def origin_cutout(self):
        """
        The ``(x, y)`` index of the origin pixel of the cutout with
        respect to the cutout array.  For ``mode='partial'``, the origin
        pixel is calculated for the valid (non-filled) cutout values.
        """
        return (self.slices_cutout[1].start, self.slices_cutout[0].start)

    @lazyproperty
    def position_original(self):
        """
        The ``(x, y)`` position index (rounded to the nearest pixel) in
        the original array.
        """
        return (_round(self.input_position_original[0]),
                _round(self.input_position_original[1]))

    @lazyproperty
    def position_cutout(self):
        """
        The ``(x, y)`` position index (rounded to the nearest pixel) in
        the cutout array.
        """
        return (_round(self.input_position_cutout[0]),
                _round(self.input_position_cutout[1]))

    @lazyproperty
    def center_original(self):
        """
        The central ``(x, y)`` position of the cutout array with respect
        to the original array.  For ``mode='partial'``, the central
        position is calculated for the valid (non-filled) cutout values.
        """
        return self._calc_center(self.slices_original)

    @lazyproperty
    def center_cutout(self):
        """
        The central ``(x, y)`` position of the cutout array with respect
        to the cutout array.  For ``mode='partial'``, the central
        position is calculated for the valid (non-filled) cutout values.
        """
        return self._calc_center(self.slices_cutout)

    @lazyproperty
    def bbox_original(self):
        """
        The bounding box ``((ymin, ymax), (xmin, xmax))`` of the minimal
        rectangular region of the cutout array with respect to the
        original array.  For ``mode='partial'``, the bounding box
        indices are for the valid (non-filled) cutout values.
        """
        return self._calc_bbox(self.slices_original)

    @lazyproperty
    def bbox_cutout(self):
        """
        The bounding box ``((ymin, ymax), (xmin, xmax))`` of the minimal
        rectangular region of the cutout array with respect to the
        cutout array.  For ``mode='partial'``, the bounding box indices
        are for the valid (non-filled) cutout values.
        """
        return self._calc_bbox(self.slices_cutout)


def make_cutouts(data, catalog, wcs=None, origin=0, verbose=True):
    """Make cutouts of catalog targets from a 2D image array.

    Notes
    -----
    The input Catalog must have the following columns, which must have
    `~astropy.unit.Unit`s where applicable:

        * ``'id'`` - ID string; no unit necessary.
        * ``'coords'`` - SkyCoord (Overrides ra, dec, x and y columns).
        * ``'ra'`` or ``'x'``- RA (angular units e.g., deg, H:M:S, arcsec etc..)
          or pixel x position (only in `~astropy.units.pix`).
        * ``'dec'`` or ``'y'`` - Dec (angular units e.g., deg, D:M:S, arcsec etc..)
          or pixel y position (only in `~astropy.units.pix`).
        * ``'cutout_width'`` - Cutout width (e.g., in arcsec, pix).
        * ``'cutout_height'`` - Cutout height (e.g., in arcsec, pix).

    Optional columns:
        * ``'cutout_pa'`` - Cutout angle (e.g., in deg, arcsec). This is only
          use if user chooses to rotate the cutouts. Positive value
          will result in a clockwise rotation.

    If saved to fits, cutouts are organized as follows:
        <output_dir>/
            <id>.fits

    Each cutout image is a simple single-extension FITS with updated WCS.
    Its header has the following special keywords:

        * ``OBJ_RA`` - RA of the cutout object in degrees.
        * ``OBJ_DEC`` - DEC of the cutout object in degrees.
        * ``OBJ_ROT`` - Rotation of cutout object in degrees.

    Parameters
    ----------
    data : 2D `~numpy.ndarray` or `~astropy.nddata.NDData`
        The 2D cutout array.
    catalog : `~astropy.table.table.Table`
        Catalog table defining the sources to cut out. Must contain
        unit information as the cutouttool does not assume default units.
    wcs : `~astropy.wcs.wcs.WCS`
        WCS if the input image is `~numpy.ndarray`.
    origin : int
        Whether SkyCoord.from_pixel should use 0 or 1-basedpixel coordinates.
    verbose : bool
        Print extra info. Default is `True`.

    Returns
    -------
    cutouts : list
        A list of NDData. If cutout failed for a target,
       `None` will be added as a place holder.
    """
    # Optional dependencies...
    try:
        from reproject.interpolation.high_level import reproject_interp
    except ImportError as e:
        raise ImportError("Optional requirement not met: " + e.msg)

    # Search for wcs:
    if isinstance(data, NDData):
            if wcs is not None:
                raise Exception("Ambiguous: WCS defined in NDData and parameters.")
            wcs = data.wcs
    elif not isinstance(data, np.ndarray):
        raise TypeError("Input image should be a 2D `~numpy.ndarray` "
                        "or `~astropy.nddata.NDData")
    elif wcs is None:
        raise Exception("WCS information was not provided.")

    # Calculate the pixel scale of input image:
    pixel_scales = proj_plane_pixel_scales(wcs)
    pixel_scale_width = pixel_scales[0] * u.Unit(wcs.wcs.cunit[0]) / u.pix
    pixel_scale_height = pixel_scales[1] * u.Unit(wcs.wcs.cunit[1]) / u.pix

    # Check if `SkyCoord`s are available:
    if 'coords' in catalog.colnames:
        coords = SkyCoord(catalog['coords'])
    elif 'ra' in catalog.colnames and 'dec' in catalog.colnames:
        if 'x' in catalog.colnames and 'y' in catalog.colnames:
            raise Exception("Ambiguous catalog: Both (ra, dec) and pixel positions provided.")
        if catalog['ra'].unit is None or catalog['dec'].unit is None:
            raise u.UnitsError("Units not specified for ra and/or dec columns.")
        coords = SkyCoord(catalog['ra'], catalog['dec'], unit=(catalog['ra'].unit,
                                                               catalog['dec'].unit))
    elif 'x' in catalog.colnames and 'y' in catalog.colnames:
        coords = SkyCoord.from_pixel(catalog['x'].astype(float), catalog['y'].astype(float), wcs, origin=origin)
    else:
        try:
            coords = SkyCoord.guess_from_table(catalog)
        except Exception as e:
            raise e

    # Figure out cutout size:
    if 'cutout_width' in catalog.colnames:
        if catalog['cutout_width'].unit is None:
            raise u.UnitsError("Units not specified for cutout_width.")
        if catalog['cutout_width'].unit == u.pix:
            width = catalog['cutout_width'].astype(float)  # pix
        else:
            width = (catalog['cutout_width'] / pixel_scale_width).decompose().value  # pix
    else:
        raise Exception("cutout_width column not found in catalog.")

    if 'cutout_height' in catalog.colnames:
        if catalog['cutout_height'].unit is None:
            raise u.UnitsError("Units not specified for cutout_height.")
        if catalog['cutout_height'].unit == u.pix:
            height = catalog['cutout_height'].astype(float)  # pix
        else:
            height = (catalog['cutout_height'] / pixel_scale_height).decompose().value  # pix
    else:
        raise Exception("cutout_height column not found in catalog.")

    # Do not rotate if column is missing.
    if 'cutout_pa' in catalog.colnames:
        if catalog['cutout_pa'].unit is None:
            raise u.UnitsError("Units not specified for cutout_pa.")
        apply_rotation = True
    else:
        apply_rotation = False

    cutcls = partial(Cutout2D, data.data, wcs=wcs, mode='partial')
    cutouts = []
    for position, x_pix, y_pix, row in zip(coords, width, height, catalog):
        if apply_rotation:
            pix_rot = row['cutout_pa'].to(u.degree).value

            # Construct new rotated WCS:
            cutout_wcs = WCS(naxis=2)
            cutout_wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']
            cutout_wcs.wcs.crval = [position.ra.deg, position.dec.deg]
            cutout_wcs.wcs.crpix = [(x_pix + 1) * 0.5, (y_pix + 1) * 0.5]

            try:
                cutout_wcs.wcs.cd = wcs.wcs.cd
                cutout_wcs.rotateCD(-pix_rot)
            except AttributeError:
                cutout_wcs.wcs.cdelt = wcs.wcs.cdelt
                cutout_wcs.wcs.crota = [0, -pix_rot]

            cutout_hdr = cutout_wcs.to_header()

            # Rotate the image using reproject
            try:
                cutout_arr = reproject_interp(
                    (data, wcs), cutout_hdr, shape_out=(math.floor(y_pix + math.copysign(0.5, y_pix)),
                                                        math.floor(x_pix + math.copysign(0.5, x_pix))), order=1)
            except Exception:
                if verbose:
                    log.info('reproject failed: '
                             'Skipping {0}'.format(row['id']))
                cutouts.append(None)
                continue

            cutout_arr = cutout_arr[0]  # Ignore footprint
            cutout_hdr['OBJ_ROT'] = (pix_rot, 'Cutout rotation in degrees')
        else:
            # Make cutout or handle exceptions by adding None to output list
            try:
                cutout = cutcls(position, size=(y_pix, x_pix))
            except NoConvergence:
                if verbose:
                    log.info('WCS solution did not converge: '
                             'Skipping {0}'.format(row['id']))
                cutouts.append(None)
                continue
            except NoOverlapError:
                if verbose:
                    log.info('Cutout is not on image: '
                             'Skipping {0}'.format(row['id']))
                cutouts.append(None)
                continue
            else:
                cutout_hdr = cutout.wcs.to_header()
                cutout_arr = cutout.data

        # If cutout result is empty, skip that target
        if np.array_equiv(cutout_arr, 0):
            if verbose:
                log.info('No data in cutout: Skipping {0}'.format(row['id']))
            cutouts.append(None)
            continue

        # Finish constructing header.
        cutout_hdr['OBJ_RA'] = (position.ra.deg, 'Cutout object RA in deg')
        cutout_hdr['OBJ_DEC'] = (position.dec.deg, 'Cutout object DEC in deg')

        cutouts.append(NDData(data=cutout_arr, wcs=WCS(cutout_hdr), meta=cutout_hdr))

    return cutouts


def cutouts_from_fits(image, catalog, image_ext=0, origin=0,
                      output_dir=None, overwrite=False, verbose=True):
    """Wrapper for the make_cutouts function. This function will take in a single
    fits image and return an array containing a list of cutouts as fits hdus. It
    will also save the cutouts to file if requested.

    Notes
    -----
    The input Catalog must have the following columns, which must have
    `~astropy.unit.Unit`s where applicable:

        * ``'id'`` - ID string; no unit necessary.
        * ``'coords'`` - SkyCoord (Overrides ra, dec, x and y columns).
        * ``'ra'`` or ``'x'``- RA (angular units e.g., deg, H:M:S, arcsec etc..)
          or pixel x position (only in `~astropy.units.pix`).
        * ``'dec'`` or ``'y'`` - Dec (angular units e.g., deg, D:M:S, arcsec etc..)
          or pixel y position (only in `~astropy.units.pix`).
        * ``'cutout_width'`` - Cutout width (e.g., in arcsec, pix).
        * ``'cutout_height'`` - Cutout height (e.g., in arcsec, pix).

    Optional columns:
        * ``'cutout_pa'`` - Cutout angle (e.g., in deg, arcsec). This is only
          use if user chooses to rotate the cutouts. Positive value
          will result in a clockwise rotation.

    If saved to fits, cutouts are organized as follows:
        <output_dir>/
            <id>.fits

    Each cutout image is a simple single-extension FITS with updated WCS.
    Its header has the following special keywords:

        * ``OBJ_RA`` - RA of the cutout object in degrees.
        * ``OBJ_DEC`` - DEC of the cutout object in degrees.
        * ``OBJ_ROT`` - Rotation of cutout object in degrees.

    Examples
    --------
    Given a list of Hubble Ultra Deep Field RA and Dec coords,
    you may use the tool as follows:
        >>> from astropy.table import Table
        >>> import astropy.units as u

        >>> ra = [53.18782913, 53.14794797, 53.15059559] * u.deg
        >>> dec = [-27.79405589, -27.77392421, -27.77158621] * u.deg
        >>> ids = ["Galax_0", 123, 53.15059559 * u.deg]
        >>> cutout_width = cutout_height = [3.0, 4.0, 3.0] * u.arcsec

        >>> catalog = Table(
        ...     data=[ids, ra, dec, cutout_width, cutout_height],
        ...     names=['id', 'ra', 'dec', 'cutout_width', 'cutout_height'])

        # To get a list of PrimaryHDU objects:
        >>> cutouts = cutouts_from_fits('h_udf_wfc_b_drz_img.fits', catalog)
        # To save to fits file provide an output dir:
        >>> cutouts = cutouts_from_fits('h_udf_wfc_b_drz_img.fits', catalog, output_dir='~/cutouts')

        # The input image can be read in before the function call:
        >>> image = fits.open('h_udf_wfc_b_drz_img.fits')
        >>> cutouts = cutouts_from_fits(image, catalog, image_ext=0)

        # If the above catalog table is saved in an ECSV file with the proper units information:
        >>> catalog.write('catalog.ecsv', format='ascii.ecsv')
        >>> cutouts = cutouts_from_fits('h_udf_wfc_b_drz_img.fits', 'catalog.ecsv')

    Parameters
    ----------
    image : filename or `HDUList` or `PrimaryHDU` or `ImageHDU` or `CompImageHDU`
        Image to cut from. If string is provided, it is assumed to be a
        fits file path.
    catalog : str or `~astropy.table.table.Table`
        Catalog table defining the sources to cut out. Must contain
        unit information as the cutouttool does not assume default units.
        Must be an astropy Table or a file name to an ECSV file containing sources.
    image_ext : int, optional
        If image is in an HDUList or read from file, use this image extension index
        to extract header and data from the primary image. Default is 0.
    origin : int
        Whether SkyCoord.from_pixel should use 0 or 1-basedpixel coordinates.
    output_dir : str
        Path to directory to save the cutouts in. If provided, each cutout will be
        saved to a separate file. The directory is created if it does not exist.
    overwrite: bool, optional
        Overwrite existing files. Default is `False`.
    verbose : bool, optional
        Print extra info. Default is `True`.

    Returns
    -------
    cutouts : list
        A list of NDData or fits PrimaryHDU. If cutout failed for a target,
       `None` will be added as a place holder.
    """

    save_to_file = output_dir is not None

    # read in the catalog file:
    if isinstance(catalog, str):
        catalog = QTable.read(catalog)
    elif not isinstance(catalog, Table):
        raise TypeError("Catalog should be an astropy.table.table.Table or"
                        " file name, got {0} instead".format(type(catalog)))

    # Load data and wcs:
    if isinstance(image, str):
        # Read data and WCS from file
        with fits.open(image) as pf:
            image_hdu = pf[image_ext]
            data = image_hdu.data
            wcs = WCS(image_hdu.header)
    else:
        # If image is HDUList or HDU:
        if isinstance(image, fits.hdu.hdulist.HDUList):
            image_hdu = image[image_ext]
        elif isinstance(image, (PrimaryHDU, ImageHDU, CompImageHDU)):
            image_hdu = image
        else:
            raise TypeError("Expected array, ImageHDU, HDUList, or file "
                            "name. Got {0} instead".format(type(image)))
        data = image_hdu.data
        wcs = WCS(image_hdu.header)

    # If image is empty, raise exception
    if np.array_equiv(data, 0):
        raise ValueError("No data in image.")

    # Sub-directory, relative to working directory.
    if save_to_file:
        path = output_dir
        if not os.path.exists(path):
            try:
                os.mkdir(path)
            except OSError as exc:
                if exc.errno != errno.EEXIST:
                    raise
                pass

    # Make the cutouts using the make_cutouts function:
    cutouts = make_cutouts(data=data, catalog=catalog, wcs=wcs,
                           origin=origin, verbose=verbose)

    # Convert cutouts to fits hdus.
    # Save the hdus to file if requested:
    fits_cutouts = []
    for idx, cutout in enumerate(cutouts):
        if cutout is None:
            fits_cutouts.append(None)
            continue

        row = catalog[idx]
        cutout_arr = cutout.data
        cutout_hdr = cutout.meta

        # Construct FITS HDU.
        hdu = fits.PrimaryHDU(cutout_arr)
        hdu.header.update(cutout_hdr)

        # Save to file if output directory is provided
        if save_to_file:
            fname = os.path.join(
                path, '{0}.fits'.format(row['id']))
            try:
                hdu.writeto(fname, overwrite=overwrite)
            except OSError as e:
                if not overwrite:
                    raise OSError(str(e)+" Try setting overwrite parameter to True.")
                else:
                    raise e
            if verbose:
                log.info('Wrote {0}'.format(fname))

        # Add cutout to output list
        fits_cutouts.append(hdu)

    return fits_cutouts
