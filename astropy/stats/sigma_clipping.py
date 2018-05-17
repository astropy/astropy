# Licensed under a 3-clause BSD style license - see LICENSE.rst

import warnings

import numpy as np

from ..utils import isiterable
from ..utils.decorators import deprecated_renamed_argument
from ..utils.exceptions import AstropyUserWarning


__all__ = ['SigmaClip', 'sigma_clip', 'sigma_clipped_stats']


class SigmaClip:
    """
    Class to perform sigma clipping.

    The data will be iterated over, each time rejecting values that are
    less or more than a specified number of standard deviations from a
    center value.

    Clipped (rejected) pixels are those where::

        data < cenfunc(data [,axis=int]) - (sigma_lower * stdfunc(data [,axis=int]))
        data > cenfunc(data [,axis=int]) + (sigma_upper * stdfunc(data [,axis=int]))

    Invalid data values (i.e. NaN or inf) are automatically clipped.

    For a functional interface to sigma clipping, see
    :func:`sigma_clip`.

    .. note::
        `scipy.stats.sigmaclip
        <https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.sigmaclip.html>`_
        provides a subset of the functionality in this class.  Also, its
        input data cannot be a masked array and it does not handle data
        that contains invalid values (i.e.  NaN or inf).  Also note that
        it uses the mean as the centering function.

        If your data is a `~numpy.ndarray` with no invalid values and
        you want to use the mean as the centering function with
        ``axis=None``, then `scipy.stats.sigmaclip` is ~25-30% faster
        than the equivalent settings here (``SigmaClip(cenfunc=np.mean,
        maxiters=None)``).

    Parameters
    ----------
    sigma : float, optional
        The number of standard deviations to use for both the lower and
        upper clipping limit.  These limits are overridden by
        ``sigma_lower`` and ``sigma_upper``, if input.  The default is
        3.

    sigma_lower : float or `None`, optional
        The number of standard deviations to use as the lower bound for
        the clipping limit.  If `None` then the value of ``sigma`` is
        used.  The default is `None`.

    sigma_upper : float or `None`, optional
        The number of standard deviations to use as the upper bound for
        the clipping limit.  If `None` then the value of ``sigma`` is
        used.  The default is `None`.

    maxiters : int or `None`, optional
        The maximum number of sigma-clipping iterations to perform or
        `None` to clip until convergence is achieved (i.e., iterate
        until the last iteration clips nothing).  If convergence is
        achieved prior to ``maxiters`` iterations, the clipping
        iterations will stop.  The default is 5.

    cenfunc : callable, optional
        The function used to compute the center value for the clipping.
        If the ``axis`` keyword is used, then it must be callable that
        can ignore NaNs (e.g. `numpy.nanmean`) and has an ``axis``
        keyword to return an array with axis dimension(s) removed.  The
        default is `numpy.nanmedian`.

    stdfunc : callable, optional
        The function used to compute the standard deviation about the
        center value.  If the ``axis`` keyword is used, then it must be
        callable that can ignore NaNs (e.g. `numpy.nanstd`) and has an
        ``axis`` keyword to return an array with axis dimension(s)
        removed.  The default is `numpy.nanstd`.

    See Also
    --------
    sigma_clip, sigma_clipped_stats

    Examples
    --------
    This example uses a data array of random variates from a Gaussian
    distribution.  We clip all points that are more than 2 sample
    standard deviations from the median.  The result is a masked array,
    where the mask is `True` for clipped data::

        >>> from astropy.stats import SigmaClip
        >>> from numpy.random import randn
        >>> randvar = randn(10000)
        >>> sigclip = SigmaClip(sigma=2, maxiters=5)
        >>> filtered_data = sigclip(randvar)

    This example clips all points that are more than 3 sigma relative to
    the sample *mean*, clips until convergence, returns an unmasked
    `~numpy.ndarray`, and does not copy the data::

        >>> from astropy.stats import SigmaClip
        >>> from numpy.random import randn
        >>> from numpy import mean
        >>> randvar = randn(10000)
        >>> sigclip = SigmaClip(sigma=3, maxiters=None, cenfunc=mean)
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

    @deprecated_renamed_argument('iters', 'maxiters', '3.1')
    def __init__(self, sigma=3., sigma_lower=None, sigma_upper=None,
                 maxiters=5, cenfunc=np.nanmedian, stdfunc=np.nanstd):

        self.sigma = sigma
        if sigma_lower is None:
            sigma_lower = sigma
        if sigma_upper is None:
            sigma_upper = sigma
        self.sigma_lower = sigma_lower
        self.sigma_upper = sigma_upper

        self.maxiters = maxiters or np.inf
        self.cenfunc = cenfunc
        self.stdfunc = stdfunc

    def __repr__(self):
        return ('SigmaClip(sigma={0}, sigma_lower={1}, sigma_upper={2}, '
                'maxiters={3}, cenfunc={4}, stdfunc={5})'
                .format(self.sigma, self.sigma_lower, self.sigma_upper,
                        self.maxiters, self.cenfunc, self.stdfunc))

    def __str__(self):
        lines = ['<' + self.__class__.__name__ + '>']
        attrs = ['sigma', 'sigma_lower', 'sigma_upper', 'maxiters', 'cenfunc',
                 'stdfunc']
        for attr in attrs:
            lines.append('    {0}: {1}'.format(attr, getattr(self, attr)))
        return '\n'.join(lines)

    def _compute_bounds(self, data, axis=None):
        # ignore RuntimeWarning if the array has only NaNs (or along an
        # axis)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            self._max_value = self.cenfunc(data, axis=axis)
            std = self.stdfunc(data, axis=axis)
            self._min_value = self._max_value - (std * self.sigma_lower)
            self._max_value += std * self.sigma_upper

    def _sigmaclip_noaxis(self, data, masked=True, copy=True):
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
            filtered_data = filtered_data[(filtered_data >= self._min_value) &
                                          (filtered_data <= self._max_value)]
            nchanged = size - filtered_data.size

        self._niterations = iteration

        if masked:
            # return a masked array
            data = np.ma.masked_invalid(data, copy=copy)

            # update the mask in place, ignoring RuntimeWarnings for
            # comparisons with NaN data values
            with np.errstate(invalid='ignore'):
                data.mask |= np.logical_or(data < self._min_value,
                                           data > self._max_value)

            return data
        else:
            # return the truncated array with the mix/max clipping thresholds
            return filtered_data, self._min_value, self._max_value

    def _sigmaclip_withaxis(self, data, axis=None, masked=True, copy=True):
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

        # convert negative axis/axes
        if not isiterable(axis):
            axis = (axis,)
        axis = tuple(filtered_data.ndim + n if n < 0 else n for n in axis)

        # define the shape of min/max arrays so that they can be broadcast
        # with the data
        mshape = tuple(1 if dim in axis else size
                       for dim, size in enumerate(filtered_data.shape))

        nchanged = 1
        iteration = 0
        while nchanged != 0 and (iteration < self.maxiters):
            iteration += 1
            size = filtered_data.size

            self._compute_bounds(filtered_data, axis=axis)
            self._min_value = self._min_value.reshape(mshape)
            self._max_value = self._max_value.reshape(mshape)

            with np.errstate(invalid='ignore'):
                filtered_data[(filtered_data < self._min_value) |
                              (filtered_data > self._max_value)] = np.nan

            nchanged = size - filtered_data.size

        self._niterations = iteration

        if masked:
            if copy:
                return np.ma.masked_invalid(filtered_data)
            else:
                # ignore RuntimeWarnings for comparisons with NaN data values
                with np.errstate(invalid='ignore'):
                    out = np.ma.masked_invalid(data, copy=False)

                    return np.ma.masked_where(np.logical_or(
                        out < self._min_value, out > self._max_value),
                        out, copy=False)
        else:
            # return the truncated array with the mix/max clipping thresholds
            return filtered_data, self._min_value, self._max_value

    def __call__(self, data, axis=None, masked=True, copy=True):
        """
        Perform sigma clipping on the provided data.

        Parameters
        ----------
        data : array-like or `~np.ma.MaskedArray`
            The data to be sigma clipped.

        axis : `None` or int or tuple of int, optional
            The axis or axes along which to sigma clip the data.  If `None`,
            then the flattened data will be used.  ``axis`` is passed
            to the ``cenfunc`` and ``stdfunc``.  The default is `None`.

        masked : bool, optional
            If `True`, then a `~numpy.ma.MaskedArray` is returned, where
            the mask is `True` for clipped values.  If `False`, then a
            `~numpy.ndarray` and the minimum and maximum clipping
            thresholds are returned.  The default is `True`.

        copy : bool, optional
            If `True`, then the ``data`` array will be copied.  If
            `False` and ``masked=True``, then the returned masked array
            data will contain the same array as the input ``data`` (if
            ``data`` is a `~numpy.ndarray` or `~numpy.ma.MaskedArray`).
            The default is `True`.

        Returns
        -------
        result : flexible
            If ``masked=True``, then a `~numpy.ma.MaskedArray` is
            returned, where the mask is `True` for clipped values.

            If ``masked=False``, then a `~numpy.ndarray` and the minimum
            and maximum clipping thresholds are returned.

            If ``masked=False`` and ``axis=None``, then the output array
            is a flattened 1D `~numpy.ndarray` where the clipped values
            have been removed and the output minimum and maximum
            thresholds are scalars.

            If ``masked=False`` and ``axis`` is specified, then the
            output `~numpy.ndarray` will have the same shape as the
            input ``data`` and contain ``np.nan`` where values were
            clipped.  In this case the returned minimum and maximum
            clipping thresholds will be also be `~numpy.ndarray`\\s.
        """

        data = np.asanyarray(data)

        if data.size == 0:
            return data

        if isinstance(data, np.ma.MaskedArray) and data.mask.all():
            return data

        if axis is None:
            return self._sigmaclip_noaxis(data, masked=masked, copy=copy)
        else:
            return self._sigmaclip_withaxis(data, axis=axis, masked=masked,
                                            copy=copy)


@deprecated_renamed_argument('iters', 'maxiters', '3.1')
def sigma_clip(data, sigma=3, sigma_lower=None, sigma_upper=None, maxiters=5,
               cenfunc=np.ma.median, stdfunc=np.std, axis=None, masked=True,
               copy=True):
    """
    Perform sigma-clipping on the provided data.

    The data will be iterated over, each time rejecting values that are
    less or more than a specified number of standard deviations from a
    center value.

    Clipped (rejected) pixels are those where::

        data < cenfunc(data [,axis=int]) - (sigma_lower * stdfunc(data [,axis=int]))
        data > cenfunc(data [,axis=int]) + (sigma_upper * stdfunc(data [,axis=int]))

    Invalid data values (i.e. NaN or inf) are automatically clipped.

    For an object-oriented interface to sigma clipping, see
    :class:`SigmaClip`.

    .. note::
        `scipy.stats.sigmaclip
        <https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.sigmaclip.html>`_
        provides a subset of the functionality in this class.  Also, its
        input data cannot be a masked array and it does not handle data
        that contains invalid values (i.e.  NaN or inf).  Also note that
        it uses the mean as the centering function.

        If your data is a `~numpy.ndarray` with no invalid values and
        you want to use the mean as the centering function with
        ``axis=None``, then `scipy.stats.sigmaclip` is ~25-30% faster
        than the equivalent settings here (``SigmaClip(cenfunc=np.mean,
        maxiters=None)``).

    Parameters
    ----------
    data : array-like or `~np.ma.MaskedArray`
        The data to be sigma clipped.

    sigma : float, optional
        The number of standard deviations to use for both the lower and
        upper clipping limit.  These limits are overridden by
        ``sigma_lower`` and ``sigma_upper``, if input.  The default is
        3.

    sigma_lower : float or `None`, optional
        The number of standard deviations to use as the lower bound for
        the clipping limit.  If `None` then the value of ``sigma`` is
        used.  The default is `None`.

    sigma_upper : float or `None`, optional
        The number of standard deviations to use as the upper bound for
        the clipping limit.  If `None` then the value of ``sigma`` is
        used.  The default is `None`.

    maxiters : int or `None`, optional
        The maximum number of sigma-clipping iterations to perform or
        `None` to clip until convergence is achieved (i.e., iterate
        until the last iteration clips nothing).  If convergence is
        achieved prior to ``maxiters`` iterations, the clipping
        iterations will stop.  The default is 5.

    cenfunc : callable, optional
        The function used to compute the center value for the clipping.
        If the ``axis`` keyword is used, then it must be callable that
        can ignore NaNs (e.g. `numpy.nanmean`) and has an ``axis``
        keyword to return an array with axis dimension(s) removed.  The
        default is `numpy.nanmedian`.

    stdfunc : callable, optional
        The function used to compute the standard deviation about the
        center value.  If the ``axis`` keyword is used, then it must be
        callable that can ignore NaNs (e.g. `numpy.nanstd`) and has an
        ``axis`` keyword to return an array with axis dimension(s)
        removed.  The default is `numpy.nanstd`.

    axis : `None` or int or tuple of int, optional
        The axis or axes along which to sigma clip the data.  If `None`,
        then the flattened data will be used.  ``axis`` is passed to the
        ``cenfunc`` and ``stdfunc``.  The default is `None`.

    masked : bool, optional
        If `True`, then a `~numpy.ma.MaskedArray` is returned, where the
        mask is `True` for clipped values.  If `False`, then a
        `~numpy.ndarray` and the minimum and maximum clipping thresholds
        are returned.  The default is `True`.

    copy : bool, optional
        If `True`, then the ``data`` array will be copied.  If `False`
        and ``masked=True``, then the returned masked array data will
        contain the same array as the input ``data`` (if ``data`` is a
        `~numpy.ndarray` or `~numpy.ma.MaskedArray`).  The default is
        `True`.

    Returns
    -------
    result : flexible
        If ``masked=True``, then a `~numpy.ma.MaskedArray` is returned,
        where the mask is `True` for clipped values.

        If ``masked=False``, then a `~numpy.ndarray` and the minimum and
        maximum clipping thresholds are returned.

        If ``masked=False`` and ``axis=None``, then the output array is
        a flattened 1D `~numpy.ndarray` where the clipped values have
        been removed and the output minimum and maximum thresholds are
        scalars.

        If ``masked=False`` and ``axis`` is specified, then the output
        `~numpy.ndarray` will have the same shape as the input ``data``
        and contain ``np.nan`` where values were clipped.  In this case
        the returned minimum and maximum clipping thresholds will be
        also be `~numpy.ndarray`\\s.

    See Also
    --------
    SigmaClip, sigma_clipped_stats

    Examples
    --------
    This example uses a data array of random variates from a Gaussian
    distribution.  We clip all points that are more than 2 sample
    standard deviations from the median.  The result is a masked array,
    where the mask is `True` for clipped data::

        >>> from astropy.stats import sigma_clip
        >>> from numpy.random import randn
        >>> randvar = randn(10000)
        >>> filtered_data = sigma_clip(randvar, sigma=2, maxiters=5)

    This example clips all points that are more than 3 sigma relative to
    the sample *mean*, clips until convergence, returns an unmasked
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
                        cenfunc=cenfunc, stdfunc=stdfunc)

    return sigclip(data, axis=axis, masked=masked, copy=copy)


@deprecated_renamed_argument('iters', 'maxiters', '3.1')
def sigma_clipped_stats(data, mask=None, mask_value=None, sigma=3.0,
                        sigma_lower=None, sigma_upper=None, maxiters=5,
                        cenfunc=np.ma.median, stdfunc=np.std, std_ddof=0,
                        axis=None):
    """
    Calculate sigma-clipped statistics on the provided data.

    Parameters
    ----------
    data : array-like or `~np.ma.MaskedArray`
        Data array or object that can be converted to an array.

    mask : `numpy.ndarray` (bool), optional
        A boolean mask with the same shape as ``data``, where a `True`
        value indicates the corresponding element of ``data`` is masked.
        Masked pixels are excluded when computing the statistics.

    mask_value : float, optional
        A data value (e.g., ``0.0``) that is ignored when computing the
        statistics.  ``mask_value`` will be masked in addition to any
        input ``mask``.

    sigma : float, optional
        The number of standard deviations to use for both the lower and
        upper clipping limit.  These limits are overridden by
        ``sigma_lower`` and ``sigma_upper``, if input.  The default is
        3.

    sigma_lower : float or `None`, optional
        The number of standard deviations to use as the lower bound for
        the clipping limit.  If `None` then the value of ``sigma`` is
        used.  The default is `None`.

    sigma_upper : float or `None`, optional
        The number of standard deviations to use as the upper bound for
        the clipping limit.  If `None` then the value of ``sigma`` is
        used.  The default is `None`.

    maxiters : int or `None`, optional
        The maximum number of sigma-clipping iterations to perform or
        `None` to clip until convergence is achieved (i.e., iterate
        until the last iteration clips nothing).  If convergence is
        achieved prior to ``maxiters`` iterations, the clipping
        iterations will stop.  The default is 5.

    cenfunc : callable, optional
        The function used to compute the center value for the clipping.
        If the ``axis`` keyword is used, then it must be callable that
        can ignore NaNs (e.g. `numpy.nanmean`) and has an ``axis``
        keyword to return an array with axis dimension(s) removed.  The
        default is `numpy.nanmedian`.

    stdfunc : callable, optional
        The function used to compute the standard deviation about the
        center value.  If the ``axis`` keyword is used, then it must be
        callable that can ignore NaNs (e.g. `numpy.nanstd`) and has an
        ``axis`` keyword to return an array with axis dimension(s)
        removed.  The default is `numpy.nanstd`.

    std_ddof : int, optional
        The delta degrees of freedom for the standard deviation
        calculation.  The divisor used in the calculation is ``N -
        std_ddof``, where ``N`` represents the number of elements.  The
        default is 0.

    axis : `None` or int or tuple of int, optional
        The axis or axes along which to sigma clip the data.  If `None`,
        then the flattened data will be used.  ``axis`` is passed
        to the ``cenfunc`` and ``stdfunc``.  The default is `None`.

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

    sigclip = SigmaClip(sigma=sigma, sigma_lower=sigma_lower,
                        sigma_upper=sigma_upper, maxiters=maxiters,
                        cenfunc=cenfunc, stdfunc=stdfunc)
    data_clipped = sigclip(data, axis=axis, masked=False, copy=False)[0]

    mean = np.nanmean(data_clipped, axis=axis)
    median = np.nanmedian(data_clipped, axis=axis)
    std = np.nanstd(data_clipped, ddof=std_ddof, axis=axis)

    return mean, median, std
