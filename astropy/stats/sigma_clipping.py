# Licensed under a 3-clause BSD style license - see LICENSE.rst

import warnings

import numpy as np
from numpy.core.multiarray import normalize_axis_index

from astropy.units import Quantity
from astropy.utils import isiterable
from astropy.utils.exceptions import AstropyUserWarning
from astropy.stats._fast_sigma_clip import _sigma_clip_fast
from astropy.stats.funcs import mad_std
from astropy.utils.compat.optional_deps import HAS_BOTTLENECK

if HAS_BOTTLENECK:
    import bottleneck


__all__ = ['SigmaClip', 'sigma_clip', 'sigma_clipped_stats']


def _move_tuple_axes_first(array, axis):
    """
    Bottleneck can only take integer axis, not tuple, so this function
    takes all the axes to be operated on and combines them into the
    first dimension of the array so that we can then use axis=0.
    """
    # Figure out how many axes we are operating over
    naxis = len(axis)

    # Add remaining axes to the axis tuple
    axis += tuple(i for i in range(array.ndim) if i not in axis)

    # The new position of each axis is just in order
    destination = tuple(range(array.ndim))

    # Reorder the array so that the axes being operated on are at the
    # beginning
    array_new = np.moveaxis(array, axis, destination)

    # Collapse the dimensions being operated on into a single dimension
    # so that we can then use axis=0 with the bottleneck functions
    array_new = array_new.reshape((-1,) + array_new.shape[naxis:])

    return array_new


def _nanmean(array, axis=None):
    """Bottleneck nanmean function that handle tuple axis."""

    if isinstance(axis, tuple):
        array = _move_tuple_axes_first(array, axis=axis)
        axis = 0

    if isinstance(array, Quantity):
        return array.__array_wrap__(bottleneck.nanmean(array, axis=axis))
    else:
        return bottleneck.nanmean(array, axis=axis)


def _nanmedian(array, axis=None):
    """Bottleneck nanmedian function that handle tuple axis."""

    if isinstance(axis, tuple):
        array = _move_tuple_axes_first(array, axis=axis)
        axis = 0

    if isinstance(array, Quantity):
        return array.__array_wrap__(bottleneck.nanmedian(array, axis=axis))
    else:
        return bottleneck.nanmedian(array, axis=axis)


def _nanstd(array, axis=None, ddof=0):
    """Bottleneck nanstd function that handle tuple axis."""

    if isinstance(axis, tuple):
        array = _move_tuple_axes_first(array, axis=axis)
        axis = 0

    if isinstance(array, Quantity):
        return array.__array_wrap__(bottleneck.nanstd(array, axis=axis,
                                                      ddof=ddof))
    else:
        return bottleneck.nanstd(array, axis=axis, ddof=ddof)


def _nanmadstd(array, axis=None):
    """mad_std function that ignores NaNs by default."""
    return mad_std(array, axis=axis, ignore_nan=True)


class SigmaClip:
    """
    Class to perform sigma clipping.

    The data will be iterated over, each time rejecting values that are
    less or more than a specified number of standard deviations from a
    center value.

    Clipped (rejected) pixels are those where::

        data < center - (sigma_lower * std)
        data > center + (sigma_upper * std)

    where::

        center = cenfunc(data [, axis=])
        std = stdfunc(data [, axis=])

    Invalid data values (i.e., NaN or inf) are automatically clipped.

    For a functional interface to sigma clipping, see
    :func:`sigma_clip`.

    .. note::
        `scipy.stats.sigmaclip` provides a subset of the functionality
        in this class. Also, its input data cannot be a masked array
        and it does not handle data that contains invalid values (i.e.,
        NaN or inf). Also note that it uses the mean as the centering
        function. The equivalent settings to `scipy.stats.sigmaclip`
        are::

            sigclip = SigmaClip(sigma=4., cenfunc='mean', maxiters=None)
            sigclip(data, axis=None, masked=False, return_bounds=True)

    Parameters
    ----------
    sigma : float, optional
        The number of standard deviations to use for both the lower
        and upper clipping limit. These limits are overridden by
        ``sigma_lower`` and ``sigma_upper``, if input. The default is 3.

    sigma_lower : float or None, optional
        The number of standard deviations to use as the lower bound for
        the clipping limit. If `None` then the value of ``sigma`` is
        used. The default is `None`.

    sigma_upper : float or None, optional
        The number of standard deviations to use as the upper bound for
        the clipping limit. If `None` then the value of ``sigma`` is
        used. The default is `None`.

    maxiters : int or None, optional
        The maximum number of sigma-clipping iterations to perform or
        `None` to clip until convergence is achieved (i.e., iterate
        until the last iteration clips nothing). If convergence is
        achieved prior to ``maxiters`` iterations, the clipping
        iterations will stop. The default is 5.

    cenfunc : {'median', 'mean'} or callable, optional
        The statistic or callable function/object used to compute
        the center value for the clipping. If using a callable
        function/object and the ``axis`` keyword is used, then it must
        be able to ignore NaNs (e.g., `numpy.nanmean`) and it must have
        an ``axis`` keyword to return an array with axis dimension(s)
        removed. The default is ``'median'``.

    stdfunc : {'std', 'mad_std'} or callable, optional
        The statistic or callable function/object used to compute the
        standard deviation about the center value. If using a callable
        function/object and the ``axis`` keyword is used, then it must
        be able to ignore NaNs (e.g., `numpy.nanstd`) and it must have
        an ``axis`` keyword to return an array with axis dimension(s)
        removed. The default is ``'std'``.

    grow : float or `False`, optional
        Radius within which to mask the neighbouring pixels of those
        that fall outwith the clipping limits (only applied along
        ``axis``, if specified). As an example, for a 2D image a value
        of 1 will mask the nearest pixels in a cross pattern around each
        deviant pixel, while 1.5 will also reject the nearest diagonal
        neighbours and so on.

    See Also
    --------
    sigma_clip, sigma_clipped_stats

    Notes
    -----
    The best performance will typically be obtained by setting
    ``cenfunc`` and ``stdfunc`` to one of the built-in functions
    specified as as string. If one of the options is set to a string
    while the other has a custom callable, you may in some cases see
    better performance if you have the `bottleneck`_ package installed.

    .. _bottleneck:  https://github.com/pydata/bottleneck

    Examples
    --------
    This example uses a data array of random variates from a Gaussian
    distribution. We clip all points that are more than 2 sample
    standard deviations from the median. The result is a masked array,
    where the mask is `True` for clipped data::

        >>> from astropy.stats import SigmaClip
        >>> from numpy.random import randn
        >>> randvar = randn(10000)
        >>> sigclip = SigmaClip(sigma=2, maxiters=5)
        >>> filtered_data = sigclip(randvar)

    This example clips all points that are more than 3 sigma relative
    to the sample *mean*, clips until convergence, returns an unmasked
    `~numpy.ndarray`, and modifies the data in-place::

        >>> from astropy.stats import SigmaClip
        >>> from numpy.random import randn
        >>> from numpy import mean
        >>> randvar = randn(10000)
        >>> sigclip = SigmaClip(sigma=3, maxiters=None, cenfunc='mean')
        >>> filtered_data = sigclip(randvar, masked=False, copy=False)

    This example sigma clips along one axis::

        >>> from astropy.stats import SigmaClip
        >>> from numpy.random import normal
        >>> from numpy import arange, diag, ones
        >>> data = arange(5) + normal(0., 0.05, (5, 5)) + diag(ones(5))
        >>> sigclip = SigmaClip(sigma=2.3)
        >>> filtered_data = sigclip(data, axis=0)

    Note that along the other axis, no points would be clipped, as the
    standard deviation is higher.
    """

    def __init__(self, sigma=3., sigma_lower=None, sigma_upper=None,
                 maxiters=5, cenfunc='median', stdfunc='std', grow=False):
        self.sigma = sigma
        self.sigma_lower = sigma_lower or sigma
        self.sigma_upper = sigma_upper or sigma
        self.maxiters = maxiters or np.inf
        self.cenfunc = cenfunc
        self.stdfunc = stdfunc
        self._cenfunc_parsed = self._parse_cenfunc(cenfunc)
        self._stdfunc_parsed = self._parse_stdfunc(stdfunc)
        self._min_value = np.nan
        self._max_value = np.nan
        self._niterations = 0
        self.grow = grow

        # This just checks that SciPy is available, to avoid failing
        # later than necessary if __call__ needs it:
        if self.grow:
            from scipy.ndimage import binary_dilation
            self._binary_dilation = binary_dilation

    def __repr__(self):
        return ('SigmaClip(sigma={}, sigma_lower={}, sigma_upper={}, '
                'maxiters={}, cenfunc={}, stdfunc={}, grow={})'
                .format(self.sigma, self.sigma_lower, self.sigma_upper,
                        self.maxiters, repr(self.cenfunc), repr(self.stdfunc),
                        self.grow))

    def __str__(self):
        lines = ['<' + self.__class__.__name__ + '>']
        attrs = ['sigma', 'sigma_lower', 'sigma_upper', 'maxiters', 'cenfunc',
                 'stdfunc', 'grow']
        for attr in attrs:
            lines.append(f'    {attr}: {repr(getattr(self, attr))}')
        return '\n'.join(lines)

    @staticmethod
    def _parse_cenfunc(cenfunc):
        if isinstance(cenfunc, str):
            if cenfunc == 'median':
                if HAS_BOTTLENECK:
                    cenfunc = _nanmedian
                else:
                    cenfunc = np.nanmedian  # pragma: no cover

            elif cenfunc == 'mean':
                if HAS_BOTTLENECK:
                    cenfunc = _nanmean
                else:
                    cenfunc = np.nanmean  # pragma: no cover

            else:
                raise ValueError(f'{cenfunc} is an invalid cenfunc.')

        return cenfunc

    @staticmethod
    def _parse_stdfunc(stdfunc):
        if isinstance(stdfunc, str):
            if stdfunc == 'std':
                if HAS_BOTTLENECK:
                    stdfunc = _nanstd
                else:
                    stdfunc = np.nanstd  # pragma: no cover
            elif stdfunc == 'mad_std':
                stdfunc = _nanmadstd
            else:
                raise ValueError(f'{stdfunc} is an invalid stdfunc.')

        return stdfunc

    def _compute_bounds(self, data, axis=None):
        # ignore RuntimeWarning if the array (or along an axis) has only
        # NaNs
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            self._max_value = self._cenfunc_parsed(data, axis=axis)
            std = self._stdfunc_parsed(data, axis=axis)
            self._min_value = self._max_value - (std * self.sigma_lower)
            self._max_value += std * self.sigma_upper

    def _sigmaclip_fast(self, data, axis=None,
                        masked=True, return_bounds=False,
                        copy=True):
        """
        Fast C implementation for simple use cases.
        """
        if isinstance(data, Quantity):
            data, unit = data.value, data.unit
        else:
            unit = None

        if copy is False and masked is False and data.dtype.kind != 'f':
            raise Exception("cannot mask non-floating-point array with NaN "
                            "values, set copy=True or masked=True to avoid "
                            "this.")

        if axis is None:
            axis = -1 if data.ndim == 1 else tuple(range(data.ndim))

        if not isiterable(axis):
            axis = normalize_axis_index(axis, data.ndim)
            data_reshaped = data
            transposed_shape = None
        else:
            # The gufunc implementation does not handle non-scalar axis
            # so we combine the dimensions together as the last
            # dimension and set axis=-1
            axis = tuple(normalize_axis_index(ax, data.ndim) for ax in axis)
            transposed_axes = tuple(ax for ax in range(data.ndim)
                                    if ax not in axis) + axis
            data_transposed = data.transpose(transposed_axes)
            transposed_shape = data_transposed.shape
            data_reshaped = data_transposed.reshape(
                transposed_shape[:data.ndim - len(axis)] + (-1,))
            axis = -1

        if data_reshaped.dtype.kind != 'f' or data_reshaped.dtype.itemsize > 8:
            data_reshaped = data_reshaped.astype(float)

        mask = ~np.isfinite(data_reshaped)
        if np.any(mask):
            warnings.warn('Input data contains invalid values (NaNs or '
                          'infs), which were automatically clipped.',
                          AstropyUserWarning)

        if isinstance(data_reshaped, np.ma.MaskedArray):
            mask |= data_reshaped.mask
            data = data.view(np.ndarray)
            data_reshaped = data_reshaped.view(np.ndarray)
            mask = np.broadcast_to(mask, data_reshaped.shape).copy()

        bound_lo, bound_hi = _sigma_clip_fast(
            data_reshaped, mask, self.cenfunc == 'median',
            self.stdfunc == 'mad_std',
            -1 if np.isinf(self.maxiters) else self.maxiters,
            self.sigma_lower, self.sigma_upper, axis=axis)

        with np.errstate(invalid='ignore'):
            mask |= data_reshaped < np.expand_dims(bound_lo, axis)
            mask |= data_reshaped > np.expand_dims(bound_hi, axis)

        if transposed_shape is not None:
            # Get mask in shape of data.
            mask = mask.reshape(transposed_shape)
            mask = mask.transpose(tuple(transposed_axes.index(ax)
                                        for ax in range(data.ndim)))

        if masked:
            result = np.ma.array(data, mask=mask, copy=copy)
        else:
            if copy:
                result = data.astype(float, copy=True)
            else:
                result = data
            result[mask] = np.nan

        if unit is not None:
            result = result << unit
            bound_lo = bound_lo << unit
            bound_hi = bound_hi << unit

        if return_bounds:
            return result, bound_lo, bound_hi
        else:
            return result

    def _sigmaclip_noaxis(self, data, masked=True, return_bounds=False,
                          copy=True):
        """
        Sigma clip when ``axis`` is None and ``grow`` is not >0.

        In this simple case, we remove clipped elements from the
        flattened array during each iteration.
        """
        filtered_data = data.ravel()

        # remove masked values and convert to ndarray
        if isinstance(filtered_data, np.ma.MaskedArray):
            filtered_data = filtered_data.data[~filtered_data.mask]

        # remove invalid values
        good_mask = np.isfinite(filtered_data)
        if np.any(~good_mask):
            filtered_data = filtered_data[good_mask]
            warnings.warn('Input data contains invalid values (NaNs or '
                          'infs), which were automatically clipped.',
                          AstropyUserWarning)

        nchanged = 1
        iteration = 0
        while nchanged != 0 and (iteration < self.maxiters):
            iteration += 1
            size = filtered_data.size
            self._compute_bounds(filtered_data, axis=None)
            filtered_data = filtered_data[
                (filtered_data >= self._min_value)
                & (filtered_data <= self._max_value)]
            nchanged = size - filtered_data.size

        self._niterations = iteration

        if masked:
            # return a masked array and optional bounds
            filtered_data = np.ma.masked_invalid(data, copy=copy)

            # update the mask in place, ignoring RuntimeWarnings for
            # comparisons with NaN data values
            with np.errstate(invalid='ignore'):
                filtered_data.mask |= np.logical_or(data < self._min_value,
                                                    data > self._max_value)

        if return_bounds:
            return filtered_data, self._min_value, self._max_value
        else:
            return filtered_data

    def _sigmaclip_withaxis(self, data, axis=None, masked=True,
                            return_bounds=False, copy=True):
        """
        Sigma clip the data when ``axis`` or ``grow`` is specified.

        In this case, we replace clipped values with NaNs as placeholder
        values.
        """
        # float array type is needed to insert nans into the array
        filtered_data = data.astype(float)    # also makes a copy

        # remove invalid values
        bad_mask = ~np.isfinite(filtered_data)
        if np.any(bad_mask):
            filtered_data[bad_mask] = np.nan
            warnings.warn('Input data contains invalid values (NaNs or '
                          'infs), which were automatically clipped.',
                          AstropyUserWarning)

        # remove masked values and convert to plain ndarray
        if isinstance(filtered_data, np.ma.MaskedArray):
            filtered_data = np.ma.masked_invalid(filtered_data).astype(float)
            filtered_data = filtered_data.filled(np.nan)

        if axis is not None:
            # convert negative axis/axes
            if not isiterable(axis):
                axis = (axis,)
            axis = tuple(filtered_data.ndim + n if n < 0 else n for n in axis)

            # define the shape of min/max arrays so that they can be broadcast
            # with the data
            mshape = tuple(1 if dim in axis else size
                           for dim, size in enumerate(filtered_data.shape))

        if self.grow:
            # Construct a growth kernel from the specified radius in
            # pixels (consider caching this for re-use by subsequent
            # calls?):
            cenidx = int(self.grow)
            size = 2 * cenidx + 1
            indices = np.mgrid[(slice(0, size),) * data.ndim]
            if axis is not None:
                for n, dim in enumerate(indices):
                    # For any axes that we're not clipping over, set
                    # their indices outside the growth radius, so masked
                    # points won't "grow" in that dimension:
                    if n not in axis:
                        dim[dim != cenidx] = size
            kernel = (sum(((idx - cenidx)**2 for idx in indices))
                      <= self.grow**2)
            del indices

        nchanged = 1
        iteration = 0
        while nchanged != 0 and (iteration < self.maxiters):
            iteration += 1
            self._compute_bounds(filtered_data, axis=axis)
            if not np.isscalar(self._min_value):
                self._min_value = self._min_value.reshape(mshape)
                self._max_value = self._max_value.reshape(mshape)

            with np.errstate(invalid='ignore'):
                # Since these comparisons are always False for NaNs, the
                # resulting mask contains only newly-rejected pixels and
                # we can dilate it without growing masked pixels more
                # than once.
                new_mask = ((filtered_data < self._min_value)
                            | (filtered_data > self._max_value))
            if self.grow:
                new_mask = self._binary_dilation(new_mask, kernel)
            filtered_data[new_mask] = np.nan
            nchanged = np.count_nonzero(new_mask)
            del new_mask

        self._niterations = iteration

        if masked:
            # create an output masked array
            if copy:
                filtered_data = np.ma.MaskedArray(data,
                                                  ~np.isfinite(filtered_data),
                                                  copy=True)
            else:
                # ignore RuntimeWarnings for comparisons with NaN data values
                with np.errstate(invalid='ignore'):
                    out = np.ma.masked_invalid(data, copy=False)

                    filtered_data = np.ma.masked_where(np.logical_or(
                        out < self._min_value, out > self._max_value),
                        out, copy=False)

        if return_bounds:
            return filtered_data, self._min_value, self._max_value
        else:
            return filtered_data

    def __call__(self, data, axis=None, masked=True, return_bounds=False,
                 copy=True):
        """
        Perform sigma clipping on the provided data.

        Parameters
        ----------
        data : array-like or `~numpy.ma.MaskedArray`
            The data to be sigma clipped.

        axis : None or int or tuple of int, optional
            The axis or axes along which to sigma clip the data. If
            `None`, then the flattened data will be used. ``axis`` is
            passed to the ``cenfunc`` and ``stdfunc``. The default is
            `None`.

        masked : bool, optional
            If `True`, then a `~numpy.ma.MaskedArray` is returned, where
            the mask is `True` for clipped values. If `False`, then a
            `~numpy.ndarray` is returned. The default is `True`.

        return_bounds : bool, optional
            If `True`, then the minimum and maximum clipping bounds are
            also returned.

        copy : bool, optional
            If `True`, then the ``data`` array will be copied. If
            `False` and ``masked=True``, then the returned masked array
            data will contain the same array as the input ``data`` (if
            ``data`` is a `~numpy.ndarray` or `~numpy.ma.MaskedArray`).
            If `False` and ``masked=False``, the input data is modified
            in-place. The default is `True`.

        Returns
        -------
        result : array-like
            If ``masked=True``, then a `~numpy.ma.MaskedArray` is
            returned, where the mask is `True` for clipped values and
            where the input mask was `True`.

            If ``masked=False``, then a `~numpy.ndarray` is returned.

            If ``return_bounds=True``, then in addition to the masked
            array or array above, the minimum and maximum clipping
            bounds are returned.

            If ``masked=False`` and ``axis=None``, then the output
            array is a flattened 1D `~numpy.ndarray` where the clipped
            values have been removed. If ``return_bounds=True`` then the
            returned minimum and maximum thresholds are scalars.

            If ``masked=False`` and ``axis`` is specified, then the
            output `~numpy.ndarray` will have the same shape as the
            input ``data`` and contain ``np.nan`` where values were
            clipped. If the input ``data`` was a masked array, then the
            output `~numpy.ndarray` will also contain ``np.nan`` where
            the input mask was `True`. If ``return_bounds=True`` then
            the returned minimum and maximum clipping thresholds will be
            be `~numpy.ndarray`\\s.
        """
        data = np.asanyarray(data)

        if data.size == 0:
            if masked:
                result = np.ma.MaskedArray(data)
            else:
                result = data

            if return_bounds:
                return result, self._min_value, self._max_value
            else:
                return result

        if isinstance(data, np.ma.MaskedArray) and data.mask.all():
            if masked:
                result = data
            else:
                result = np.full(data.shape, np.nan)

            if return_bounds:
                return result, self._min_value, self._max_value
            else:
                return result

        # Shortcut for common cases where a fast C implementation can be
        # used.
        if (self.cenfunc in ('mean', 'median')
                and self.stdfunc in ('std', 'mad_std')
                and axis is not None and not self.grow):
            return self._sigmaclip_fast(data, axis=axis, masked=masked,
                                        return_bounds=return_bounds,
                                        copy=copy)

        # These two cases are treated separately because when
        # ``axis=None`` we can simply remove clipped values from the
        # array. This is not possible when ``axis`` or ``grow`` is
        # specified.
        if axis is None and not self.grow:
            return self._sigmaclip_noaxis(data, masked=masked,
                                          return_bounds=return_bounds,
                                          copy=copy)
        else:
            return self._sigmaclip_withaxis(data, axis=axis, masked=masked,
                                            return_bounds=return_bounds,
                                            copy=copy)


def sigma_clip(data, sigma=3, sigma_lower=None, sigma_upper=None, maxiters=5,
               cenfunc='median', stdfunc='std', axis=None, masked=True,
               return_bounds=False, copy=True, grow=False):
    """
    Perform sigma-clipping on the provided data.

    The data will be iterated over, each time rejecting values that are
    less or more than a specified number of standard deviations from a
    center value.

    Clipped (rejected) pixels are those where::

        data < center - (sigma_lower * std)
        data > center + (sigma_upper * std)

    where::

        center = cenfunc(data [, axis=])
        std = stdfunc(data [, axis=])

    Invalid data values (i.e., NaN or inf) are automatically clipped.

    For an object-oriented interface to sigma clipping, see
    :class:`SigmaClip`.

    .. note::
        `scipy.stats.sigmaclip` provides a subset of the functionality
        in this class. Also, its input data cannot be a masked array
        and it does not handle data that contains invalid values (i.e.,
        NaN or inf). Also note that it uses the mean as the centering
        function. The equivalent settings to `scipy.stats.sigmaclip`
        are::

            sigma_clip(sigma=4., cenfunc='mean', maxiters=None, axis=None,
            ...        masked=False, return_bounds=True)

    Parameters
    ----------
    data : array-like or `~numpy.ma.MaskedArray`
        The data to be sigma clipped.

    sigma : float, optional
        The number of standard deviations to use for both the lower
        and upper clipping limit. These limits are overridden by
        ``sigma_lower`` and ``sigma_upper``, if input. The default is 3.

    sigma_lower : float or None, optional
        The number of standard deviations to use as the lower bound for
        the clipping limit. If `None` then the value of ``sigma`` is
        used. The default is `None`.

    sigma_upper : float or None, optional
        The number of standard deviations to use as the upper bound for
        the clipping limit. If `None` then the value of ``sigma`` is
        used. The default is `None`.

    maxiters : int or None, optional
        The maximum number of sigma-clipping iterations to perform or
        `None` to clip until convergence is achieved (i.e., iterate
        until the last iteration clips nothing). If convergence is
        achieved prior to ``maxiters`` iterations, the clipping
        iterations will stop. The default is 5.

    cenfunc : {'median', 'mean'} or callable, optional
        The statistic or callable function/object used to compute
        the center value for the clipping. If using a callable
        function/object and the ``axis`` keyword is used, then it must
        be able to ignore NaNs (e.g., `numpy.nanmean`) and it must have
        an ``axis`` keyword to return an array with axis dimension(s)
        removed. The default is ``'median'``.

    stdfunc : {'std', 'mad_std'} or callable, optional
        The statistic or callable function/object used to compute the
        standard deviation about the center value. If using a callable
        function/object and the ``axis`` keyword is used, then it must
        be able to ignore NaNs (e.g., `numpy.nanstd`) and it must have
        an ``axis`` keyword to return an array with axis dimension(s)
        removed. The default is ``'std'``.

    axis : None or int or tuple of int, optional
        The axis or axes along which to sigma clip the data. If `None`,
        then the flattened data will be used. ``axis`` is passed to the
        ``cenfunc`` and ``stdfunc``. The default is `None`.

    masked : bool, optional
        If `True`, then a `~numpy.ma.MaskedArray` is returned, where
        the mask is `True` for clipped values. If `False`, then a
        `~numpy.ndarray` and the minimum and maximum clipping thresholds
        are returned. The default is `True`.

    return_bounds : bool, optional
        If `True`, then the minimum and maximum clipping bounds are also
        returned.

    copy : bool, optional
        If `True`, then the ``data`` array will be copied. If `False`
        and ``masked=True``, then the returned masked array data will
        contain the same array as the input ``data`` (if ``data`` is a
        `~numpy.ndarray` or `~numpy.ma.MaskedArray`). If `False` and
        ``masked=False``, the input data is modified in-place. The
        default is `True`.

    grow : float or `False`, optional
        Radius within which to mask the neighbouring pixels of those
        that fall outwith the clipping limits (only applied along
        ``axis``, if specified). As an example, for a 2D image a value
        of 1 will mask the nearest pixels in a cross pattern around each
        deviant pixel, while 1.5 will also reject the nearest diagonal
        neighbours and so on.

    Returns
    -------
    result : array-like
        If ``masked=True``, then a `~numpy.ma.MaskedArray` is returned,
        where the mask is `True` for clipped values and where the input
        mask was `True`.

        If ``masked=False``, then a `~numpy.ndarray` is returned.

        If ``return_bounds=True``, then in addition to the masked array
        or array above, the minimum and maximum clipping bounds are
        returned.

        If ``masked=False`` and ``axis=None``, then the output array
        is a flattened 1D `~numpy.ndarray` where the clipped values
        have been removed. If ``return_bounds=True`` then the returned
        minimum and maximum thresholds are scalars.

        If ``masked=False`` and ``axis`` is specified, then the output
        `~numpy.ndarray` will have the same shape as the input ``data``
        and contain ``np.nan`` where values were clipped. If the input
        ``data`` was a masked array, then the output `~numpy.ndarray`
        will also contain ``np.nan`` where the input mask was `True`.
        If ``return_bounds=True`` then the returned minimum and maximum
        clipping thresholds will be be `~numpy.ndarray`\\s.

    See Also
    --------
    SigmaClip, sigma_clipped_stats

    Notes
    -----
    The best performance will typically be obtained by setting
    ``cenfunc`` and ``stdfunc`` to one of the built-in functions
    specified as as string. If one of the options is set to a string
    while the other has a custom callable, you may in some cases see
    better performance if you have the `bottleneck`_ package installed.

    .. _bottleneck:  https://github.com/pydata/bottleneck

    Examples
    --------
    This example uses a data array of random variates from a Gaussian
    distribution. We clip all points that are more than 2 sample
    standard deviations from the median. The result is a masked array,
    where the mask is `True` for clipped data::

        >>> from astropy.stats import sigma_clip
        >>> from numpy.random import randn
        >>> randvar = randn(10000)
        >>> filtered_data = sigma_clip(randvar, sigma=2, maxiters=5)

    This example clips all points that are more than 3 sigma relative
    to the sample *mean*, clips until convergence, returns an unmasked
    `~numpy.ndarray`, and does not copy the data::

        >>> from astropy.stats import sigma_clip
        >>> from numpy.random import randn
        >>> from numpy import mean
        >>> randvar = randn(10000)
        >>> filtered_data = sigma_clip(randvar, sigma=3, maxiters=None,
        ...                            cenfunc=mean, masked=False, copy=False)

    This example sigma clips along one axis::

        >>> from astropy.stats import sigma_clip
        >>> from numpy.random import normal
        >>> from numpy import arange, diag, ones
        >>> data = arange(5) + normal(0., 0.05, (5, 5)) + diag(ones(5))
        >>> filtered_data = sigma_clip(data, sigma=2.3, axis=0)

    Note that along the other axis, no points would be clipped, as the
    standard deviation is higher.
    """
    sigclip = SigmaClip(sigma=sigma, sigma_lower=sigma_lower,
                        sigma_upper=sigma_upper, maxiters=maxiters,
                        cenfunc=cenfunc, stdfunc=stdfunc, grow=grow)

    return sigclip(data, axis=axis, masked=masked,
                   return_bounds=return_bounds, copy=copy)


def sigma_clipped_stats(data, mask=None, mask_value=None, sigma=3.0,
                        sigma_lower=None, sigma_upper=None, maxiters=5,
                        cenfunc='median', stdfunc='std', std_ddof=0,
                        axis=None, grow=False):
    """
    Calculate sigma-clipped statistics on the provided data.

    Parameters
    ----------
    data : array-like or `~numpy.ma.MaskedArray`
        Data array or object that can be converted to an array.

    mask : `numpy.ndarray` (bool), optional
        A boolean mask with the same shape as ``data``, where a `True`
        value indicates the corresponding element of ``data`` is masked.
        Masked pixels are excluded when computing the statistics.

    mask_value : float, optional
        A data value (e.g., ``0.0``) that is ignored when computing the
        statistics. ``mask_value`` will be masked in addition to any
        input ``mask``.

    sigma : float, optional
        The number of standard deviations to use for both the lower
        and upper clipping limit. These limits are overridden by
        ``sigma_lower`` and ``sigma_upper``, if input. The default is 3.

    sigma_lower : float or None, optional
        The number of standard deviations to use as the lower bound for
        the clipping limit. If `None` then the value of ``sigma`` is
        used. The default is `None`.

    sigma_upper : float or None, optional
        The number of standard deviations to use as the upper bound for
        the clipping limit. If `None` then the value of ``sigma`` is
        used. The default is `None`.

    maxiters : int or None, optional
        The maximum number of sigma-clipping iterations to perform or
        `None` to clip until convergence is achieved (i.e., iterate
        until the last iteration clips nothing). If convergence is
        achieved prior to ``maxiters`` iterations, the clipping
        iterations will stop. The default is 5.

    cenfunc : {'median', 'mean'} or callable, optional
        The statistic or callable function/object used to compute
        the center value for the clipping. If using a callable
        function/object and the ``axis`` keyword is used, then it must
        be able to ignore NaNs (e.g., `numpy.nanmean`) and it must have
        an ``axis`` keyword to return an array with axis dimension(s)
        removed. The default is ``'median'``.

    stdfunc : {'std', 'mad_std'} or callable, optional
        The statistic or callable function/object used to compute the
        standard deviation about the center value. If using a callable
        function/object and the ``axis`` keyword is used, then it must
        be able to ignore NaNs (e.g., `numpy.nanstd`) and it must have
        an ``axis`` keyword to return an array with axis dimension(s)
        removed. The default is ``'std'``.

    std_ddof : int, optional
        The delta degrees of freedom for the standard deviation
        calculation. The divisor used in the calculation is ``N -
        std_ddof``, where ``N`` represents the number of elements. The
        default is 0.

    axis : None or int or tuple of int, optional
        The axis or axes along which to sigma clip the data. If `None`,
        then the flattened data will be used. ``axis`` is passed to the
        ``cenfunc`` and ``stdfunc``. The default is `None`.

    grow : float or `False`, optional
        Radius within which to mask the neighbouring pixels of those
        that fall outwith the clipping limits (only applied along
        ``axis``, if specified). As an example, for a 2D image a value
        of 1 will mask the nearest pixels in a cross pattern around each
        deviant pixel, while 1.5 will also reject the nearest diagonal
        neighbours and so on.

    Notes
    -----
    The best performance will typically be obtained by setting
    ``cenfunc`` and ``stdfunc`` to one of the built-in functions
    specified as as string. If one of the options is set to a string
    while the other has a custom callable, you may in some cases see
    better performance if you have the `bottleneck`_ package installed.

    .. _bottleneck:  https://github.com/pydata/bottleneck

    Returns
    -------
    mean, median, stddev : float
        The mean, median, and standard deviation of the sigma-clipped
        data.

    See Also
    --------
    SigmaClip, sigma_clip
    """
    if mask is not None:
        data = np.ma.MaskedArray(data, mask)
    if mask_value is not None:
        data = np.ma.masked_values(data, mask_value)

    if isinstance(data, np.ma.MaskedArray) and data.mask.all():
        return np.ma.masked, np.ma.masked, np.ma.masked

    sigclip = SigmaClip(sigma=sigma, sigma_lower=sigma_lower,
                        sigma_upper=sigma_upper, maxiters=maxiters,
                        cenfunc=cenfunc, stdfunc=stdfunc, grow=grow)
    data_clipped = sigclip(data, axis=axis, masked=False, return_bounds=False,
                           copy=True)

    if HAS_BOTTLENECK:
        mean = _nanmean(data_clipped, axis=axis)
        median = _nanmedian(data_clipped, axis=axis)
        std = _nanstd(data_clipped, ddof=std_ddof, axis=axis)
    else:  # pragma: no cover
        mean = np.nanmean(data_clipped, axis=axis)
        median = np.nanmedian(data_clipped, axis=axis)
        std = np.nanstd(data_clipped, ddof=std_ddof, axis=axis)

    return mean, median, std
