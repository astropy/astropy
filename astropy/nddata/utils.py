# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module includes helper functions for array operations.
"""

from copy import deepcopy

import numpy as np

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.utils import lazyproperty
from astropy.wcs.utils import skycoord_to_pixel, proj_plane_pixel_scales
from astropy.wcs import Sip

__all__ = ['extract_array', 'add_array', 'subpixel_indices',
           'overlap_slices', 'NoOverlapError', 'PartialOverlapError',
           'Cutout2D']


class NoOverlapError(ValueError):
    '''Raised when determining the overlap of non-overlapping arrays.'''
    pass


class PartialOverlapError(ValueError):
    '''Raised when arrays only partially overlap.'''
    pass


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
    small_array_shape : int or tuple thereof
        The shape of the small array (for 1D arrays, this can be an
        `int`).  See the ``mode`` keyword for additional details.
    position : number or tuple thereof
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
    slices_large : tuple of slice
        A tuple of slice objects for each axis of the large array, such
        that ``large_array[slices_large]`` extracts the region of the
        large array that overlaps with the small array.
    slices_small : tuple of slice
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

    if any(~np.isfinite(position)):
        raise ValueError('Input position contains invalid values (NaNs or '
                         'infs).')

    if len(small_array_shape) != len(large_array_shape):
        raise ValueError('"large_array_shape" and "small_array_shape" must '
                         'have the same number of dimensions.')

    if len(small_array_shape) != len(position):
        raise ValueError('"position" must have the same number of dimensions '
                         'as "small_array_shape".')

    # define the min/max pixel indices
    indices_min = [int(np.ceil(pos - (small_shape / 2.)))
                   for (pos, small_shape) in zip(position, small_array_shape)]
    indices_max = [int(np.ceil(pos + (small_shape / 2.)))
                   for (pos, small_shape) in zip(position, small_array_shape)]

    for e_max in indices_max:
        if e_max < 0:
            raise NoOverlapError('Arrays do not overlap.')
    for e_min, large_shape in zip(indices_min, large_array_shape):
        if e_min >= large_shape:
            raise NoOverlapError('Arrays do not overlap.')

    if mode == 'strict':
        for e_min in indices_min:
            if e_min < 0:
                raise PartialOverlapError('Arrays overlap only partially.')
        for e_max, large_shape in zip(indices_max, large_array_shape):
            if e_max > large_shape:
                raise PartialOverlapError('Arrays overlap only partially.')

    # Set up slices
    slices_large = tuple(slice(max(0, indices_min),
                               min(large_shape, indices_max))
                         for (indices_min, indices_max, large_shape) in
                         zip(indices_min, indices_max, large_array_shape))
    if mode == 'trim':
        slices_small = tuple(slice(0, slc.stop - slc.start)
                             for slc in slices_large)
    else:
        slices_small = tuple(slice(max(0, -indices_min),
                                   min(large_shape - indices_min,
                                       indices_max - indices_min))
                             for (indices_min, indices_max, large_shape) in
                             zip(indices_min, indices_max, large_array_shape))

    return slices_large, slices_small


def extract_array(array_large, shape, position, mode='partial',
                  fill_value=np.nan, return_position=False):
    """
    Extract a smaller array of the given shape and position from a
    larger array.

    Parameters
    ----------
    array_large : ndarray
        The array from which to extract the small array.
    shape : int or tuple thereof
        The shape of the extracted array (for 1D arrays, this can be an
        `int`).  See the ``mode`` keyword for additional details.
    position : number or tuple thereof
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
        ``fill_value`` will be changed to have the same ``dtype`` as the
        ``array_large`` array, with one exception. If ``array_large``
        has integer type and ``fill_value`` is ``np.nan``, then a
        `ValueError` will be raised.
    return_position : bool, optional
        If `True`, return the coordinates of ``position`` in the
        coordinate system of the returned array.

    Returns
    -------
    array_small : ndarray
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
        try:
            extracted_array[:] = fill_value
        except ValueError as exc:
            exc.args += ('fill_value is inconsistent with the data type of '
                         'the input array (e.g., fill_value cannot be set to '
                         'np.nan if the input array has integer type). Please '
                         'change either the input array dtype or the '
                         'fill_value.',)
            raise exc

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
    array_large : ndarray
        Large array.
    array_small : ndarray
        Small array to add. Can be equal to ``array_large`` in size in a given
        dimension, but not larger.
    position : tuple
        Position of the small array's center, with respect to the large array.
        Coordinates should be in the same order as the array shape.

    Returns
    -------
    new_array : ndarray
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
    # Check if large array is not smaller
    if all(large_shape >= small_shape for (large_shape, small_shape)
           in zip(array_large.shape, array_small.shape)):
        large_slices, small_slices = overlap_slices(array_large.shape,
                                                    array_small.shape,
                                                    position)
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
    position : ndarray or array-like
        Positions in pixels.
    subsampling : int
        Subsampling factor per pixel.

    Returns
    -------
    indices : ndarray
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

    For example usage, see :ref:`astropy:cutout_images`.

    .. warning::

        The cutout WCS object does not currently handle cases where the
        input WCS object contains distortion lookup tables described in
        the `FITS WCS distortion paper
        <https://www.atnf.csiro.au/people/mcalabre/WCS/dcs_20040422.pdf>`__.

    Parameters
    ----------
    data : ndarray
        The 2D data array from which to extract the cutout array.

    position : tuple or `~astropy.coordinates.SkyCoord`
        The position of the cutout array's center with respect to
        the ``data`` array.  The position can be specified either as
        a ``(x, y)`` tuple of pixel coordinates or a
        `~astropy.coordinates.SkyCoord`, in which case ``wcs`` is a
        required input.

    size : int, array-like, or `~astropy.units.Quantity`
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

    fill_value : float or int, optional
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

    shape : (2,) tuple
        The ``(ny, nx)`` shape of the cutout array.

    shape_input : (2,) tuple
        The ``(ny, nx)`` shape of the input (original) array.

    input_position_cutout : (2,) tuple
        The (unrounded) ``(x, y)`` position with respect to the cutout
        array.

    input_position_original : (2,) tuple
        The original (unrounded) ``(x, y)`` input position (with respect
        to the original array).

    slices_original : (2,) tuple of slice object
        A tuple of slice objects for the minimal bounding box of the
        cutout with respect to the original array.  For
        ``mode='partial'``, the slices are for the valid (non-filled)
        cutout values.

    slices_cutout : (2,) tuple of slice object
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

    wcs : `~astropy.wcs.WCS` or None
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
        if wcs is None:
            wcs = getattr(data, 'wcs', None)

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
            self.wcs.array_shape = self.data.shape
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

    @staticmethod
    def _round(a):
        """
        Round the input to the nearest integer.

        If two integers are equally close, the value is rounded up.
        Note that this is different from `np.round`, which rounds to the
        nearest even number.
        """
        return int(np.floor(a + 0.5))

    @lazyproperty
    def position_original(self):
        """
        The ``(x, y)`` position index (rounded to the nearest pixel) in
        the original array.
        """
        return (self._round(self.input_position_original[0]),
                self._round(self.input_position_original[1]))

    @lazyproperty
    def position_cutout(self):
        """
        The ``(x, y)`` position index (rounded to the nearest pixel) in
        the cutout array.
        """
        return (self._round(self.input_position_cutout[0]),
                self._round(self.input_position_cutout[1]))

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
