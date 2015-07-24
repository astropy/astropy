# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module includes helper functions for array operations.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import numpy as np
from copy import deepcopy
from .decorators import support_nddata
from astropy.utils import lazyproperty
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import skycoord_to_pixel


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
    large_array_shape : tuple or int
        The shape of the large array (for 1D arrays, this can be an
        `int`).
    small_array_shape : tuple or int
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
    slices_small : slice
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
        large array.  The pixel cordinates should be in the same order
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


class Cutout2D(object):
    def __init__(self, data, position, shape, wcs=None, mode='trim',
                 fill_value=np.nan, copy=False):
        """
        Create a cutout object from a 2D array.

        The returned object will contain a 2D cutout array.  If
        ``copy=False`` (default), the cutout array is a view into the
        original ``data`` array, otherwise the cutout array will contain
        a copy of the original data.

        If a `~astropy.wcs.WCS` object is input, then the returned
        object will also contain a copy of the original WCS, but updated
        for the cutout array.

        For example usage, see :ref:`cutout_images`.

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

        shape : tuple
            The shape (``(ny, nx)``) of the cutout array in pixel
            coordinates (but see the ``mode`` keyword for additional
            details).

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

        Returns
        -------
        result : `~astropy.nddata.utils.Cutout2D`
            A cutout object containing the 2D cutout data array and the
            updated WCS, if ``wcs`` is input.

        Examples
        --------
        >>> import numpy as np
        >>> from astropy.nddata.utils import Cutout2D
        >>> data = np.arange(20.).reshape(5, 4)
        >>> c1 = Cutout2D(data, (2, 2), (3, 3))
        >>> print(c1.data)
        [[  5.   6.   7.]
         [  9.  10.  11.]
         [ 13.  14.  15.]]

        >>> print(c1.center_original)
        (2.0, 2.0)
        >>> print(c1.center_cutout)
        (1.0, 1.0)
        >>> print(c1.origin_original)
        (1, 1)

        >>> c2 = Cutout2D(data, (0, 0), (3, 3))
        >>> print(c2.data)
        [[ 0.  1.]
         [ 4.  5.]]

        >>> c3 = Cutout2D(data, (0, 0), (3, 3), mode='partial')
        >>> print(c3.data)
        [[ nan  nan  nan]
         [ nan   0.   1.]
         [ nan   4.   5.]]
        """

        if isinstance(position, SkyCoord):
            if wcs is None:
                raise ValueError('wcs must be input if position is a '
                                 'SkyCoord')
            position = skycoord_to_pixel(position, wcs, mode='all')  # (x, y)

        # extract_array and overlap_slices use (y, x) positions
        pos = position[::-1]

        data = np.asanyarray(data)
        cutout_data, input_position_cutout  = extract_array(
            data, shape, pos, mode=mode, fill_value=fill_value,
            return_position=True)
        if copy:
            cutout_data = np.copy(cutout_data)
        self.data = cutout_data

        self.input_position_cutout = input_position_cutout[::-1]    # (x, y)
        slices_original, slices_cutout = overlap_slices(data.shape, shape,
                                                    pos, mode=mode)
        self.slices_original = slices_original
        self.slices_cutout = slices_cutout

        self.shape = self.data.shape
        self.input_position_original = position
        self.shape_input = shape

        ((self.xmin_original, self.xmax_original),
         (self.ymin_original, self.ymax_original)) = self.bbox_original

        ((self.xmin_cutout, self.xmax_cutout),
         (self.ymin_cutout, self.ymax_cutout)) = self.bbox_cutout

        # the true origin pixel of the cutout array, including any
        # filled cutout values
        self._origin_original_true = (self.origin_original[0] -
                                   self.slices_cutout[1].start,
                                   self.origin_original[1] -
                                   self.slices_cutout[0].start)

        if wcs is not None:
            self.wcs = deepcopy(wcs)
            self.wcs.wcs.crpix -= self._origin_original_true
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
            `ax=None`.  Otherwise the output ``ax`` is the same as the
            input ``ax``.
        """

        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches

        kwargs['fill'] = fill

        if ax is None:
            ax = plt.gca()

        height, width = self.shape
        hw, hh = width / 2., height / 2.
        pos = self.position_original - np.array([hw, hh])
        patch = mpatches.Rectangle(pos, width, height, 0., **kwargs)
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
        Calculate a minimal bounding box in the form ``((xmin, xmax),
        (ymin, ymax))``.  Note these are pixel locations, not slice
        indices.  For ``mode='partial'``, the bounding box indices are
        for the valid (non-filled) cutout values.
        """
        # (stop - 1) to return the max pixel location, not the slice index
        return ((slices[1].start, slices[1].stop - 1),
                (slices[0].start, slices[0].stop - 1))

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
        The bounding box ``(ymin, xmin, ymax, xmax)`` of the minimal
        rectangular region of the cutout array with respect to the
        original array.  For ``mode='partial'``, the bounding box
        indices are for the valid (non-filled) cutout values.
        """
        return self._calc_bbox(self.slices_original)

    @lazyproperty
    def bbox_cutout(self):
        """
        The bounding box ``(ymin, xmin, ymax, xmax)`` of the minimal
        rectangular region of the cutout array with respect to the
        cutout array.  For ``mode='partial'``, the bounding box indices
        are for the valid (non-filled) cutout values.
        """
        return self._calc_bbox(self.slices_cutout)
