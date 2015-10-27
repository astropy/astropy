# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import numpy as np
import warnings
from ..utils.exceptions import AstropyDeprecationWarning, AstropyUserWarning

try:
    from scipy.stats import sigmaclip
    HAS_SCIPY=True
except ImportError:
    warnings.warn("Scipy is required for some functionality of"
            "sigma_clip. The `masked` keyword will be forced to `True`.",
            AstropyUserWarning)
    HAS_SCIPY=False


__all__ = ['sigma_clip', 'sigma_clipped_stats']


def sigma_clip(data, **kwargs):
    # temporary function to handle deprecated and removed keywords
    if 'sig' in kwargs:
        warnings.warn('The "sig" keyword is now deprecated, use the '
                      '"sigma" keyword instead.', AstropyDeprecationWarning)

        if 'sigma' not in kwargs:
            kwargs['sigma'] = kwargs['sig']
        else:
            warnings.warn('Both the "sig" and "sigma" keywords were set. '
                          'Using the value of "sigma".', AstropyUserWarning)
        del kwargs['sig']

    if 'varfunc' in kwargs:
        raise SyntaxError('The "varfunc" keyword is no longer supported. '
                          'Please use the "stdfunc" keyword instead.')

    return _sigma_clip(data, **kwargs)


def _sigma_clip(data, sigma=3, sigma_lower=None, sigma_upper=None, iters=5,
                cenfunc=np.ma.median, stdfunc=np.std, axis=None, copy=True,
                masked=True):
    """
    sigma_clip(data, sigma=3, sigma_lower=None, sigma_upper=None, iters=5, cenfunc=np.ma.median, stdfunc=np.std, axis=None, copy=True, masked=True)

    Perform sigma-clipping on the provided data.

    The data will be iterated over, each time rejecting points that are
    discrepant by more than a specified number of standard deviations from a
    center value. If the data contains invalid values (NaNs or infs),
    they are automatically masked before performing the sigma clipping.

    .. note::
        `scipy.stats.sigmaclip
        <http://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.sigmaclip.html>`_
        provides a subset of the functionality in this function.

    Parameters
    ----------
    data : array-like
        The data to be sigma clipped.
    sigma : float, optional
        The number of standard deviations to use for both the lower and
        upper clipping limit. These limits are overridden by
        ``sigma_lower`` and ``sigma_upper``, if input. Defaults to 3.
    sigma_lower : float or `None`, optional
        The number of standard deviations to use as the lower bound for
        the clipping limit. If `None` then the value of ``sigma`` is
        used. Defaults to `None`.
    sigma_upper : float or `None`, optional
        The number of standard deviations to use as the upper bound for
        the clipping limit. If `None` then the value of ``sigma`` is
        used. Defaults to `None`.
    iters : int or `None`, optional
        The number of iterations to perform sigma clipping, or `None` to
        clip until convergence is achieved (i.e., continue until the
        last iteration clips nothing). Defaults to 5.
    cenfunc : callable, optional
        The function used to compute the center for the clipping. Must
        be a callable that takes in a masked array and outputs the
        central value. Defaults to the median (`numpy.ma.median`).
    stdfunc : callable, optional
        The function used to compute the standard deviation about the
        center. Must be a callable that takes in a masked array and
        outputs a width estimator. Masked (rejected) pixels are those
        where::

             deviation < (-sigma_lower * stdfunc(deviation))
             deviation > (sigma_upper * stdfunc(deviation))

        where::

            deviation = data - cenfunc(data [,axis=int])

        Defaults to the standard deviation (`numpy.std`).
    axis : int or `None`, optional
        If not `None`, clip along the given axis.  For this case,
        ``axis`` will be passed on to ``cenfunc`` and ``stdfunc``, which
        are expected to return an array with the axis dimension removed
        (like the numpy functions).  If `None`, clip over all axes.
        Defaults to `None`.
    copy : bool, optional
        If `True`, the ``data`` array will be copied.  If `False`, the
        returned masked array data will contain the same array as
        ``data``.  Defaults to `True`.
    masked: bool, optional
        if `True` will use Numpy masked arrays to mask points rejected by
        sigma clipping. If `False` will use the `scipy.stats.sigmaclip` 
        function which throws away points rejected by sigma clipping. Setting
        it to `False` can result in significant speedups. Defaults to `True`.

    Returns
    -------
    filtered_data : `numpy.ma.MaskedArray`
        A masked array with the same shape as ``data`` input, where the
        points rejected by the algorithm have been masked.

    Notes
    -----
     1. The routine works by calculating::

            deviation = data - cenfunc(data [,axis=int])

        and then setting a mask for points outside the range::

           deviation < (-sigma_lower * stdfunc(deviation))
           deviation > (sigma_upper * stdfunc(deviation))

        It will iterate a given number of times, or until no further
        data are rejected.

     2. Most numpy functions deal well with masked arrays, but if one
        would like to have an array with just the good (or bad) values, one
        can use::

            good_only = filtered_data.data[~filtered_data.mask]
            bad_only = filtered_data.data[filtered_data.mask]

        However, for multidimensional data, this flattens the array,
        which may not be what one wants (especially if filtering was
        done along an axis).

    Examples
    --------
    This example generates random variates from a Gaussian distribution
    and returns a masked array in which all points that are more than 2
    sample standard deviations from the median are masked::

        >>> from astropy.stats import sigma_clip
        >>> from numpy.random import randn
        >>> randvar = randn(10000)
        >>> filtered_data = sigma_clip(randvar, sigma=2, iters=5)

    This example sigma clips on a similar distribution, but uses 3 sigma
    relative to the sample *mean*, clips until convergence, and does not
    copy the data::

        >>> from astropy.stats import sigma_clip
        >>> from numpy.random import randn
        >>> from numpy import mean
        >>> randvar = randn(10000)
        >>> filtered_data = sigma_clip(randvar, sigma=3, iters=None,
        ...                            cenfunc=mean, copy=False)

    This example sigma clips along one axis on a similar distribution
    (with bad points inserted)::

        >>> from astropy.stats import sigma_clip
        >>> from numpy.random import normal
        >>> from numpy import arange, diag, ones
        >>> data = arange(5) + normal(0., 0.05, (5, 5)) + diag(ones(5))
        >>> filtered_data = sigma_clip(data, sigma=2.3, axis=0)

    Note that along the other axis, no points would be masked, as the
    variance is higher.
    """
    def perform_clip(_filtered_data, _kwargs):
        """
        Perform sigma clip by comparing the data to the minimum and maximum
        values (median + sig * standard deviation). Use sigma_lower and
        sigma_upper to get the correct limits. Data values less or greater
        than the minimum / maximum values will have True set in the mask array.
        """
        max_value = cenfunc(_filtered_data, **_kwargs)
        std = stdfunc(_filtered_data, **_kwargs)
        min_value = max_value - std * sigma_lower
        max_value += std * sigma_upper
        if axis is not None:
            if axis > 0:
                min_value = np.expand_dims(min_value, axis=axis)
                max_value = np.expand_dims(max_value, axis=axis)
        _filtered_data.mask |= _filtered_data > max_value
        _filtered_data.mask |= _filtered_data < min_value

    if sigma_lower is None:
        sigma_lower = sigma
    if sigma_upper is None:
        sigma_upper = sigma

    kwargs = dict()

    if axis is not None:
        kwargs['axis'] = axis

    if np.any(~np.isfinite(data)):
        data = np.ma.masked_invalid(data)
        warnings.warn("Input data contains invalid values (NaNs or infs), "
                      "which were automatically masked.", AstropyUserWarning)

    if HAS_SCIPY is True:
        if masked is False:
            if cenfunc is not np.ma.mean or stdfunc is not np.std:
                warnings.warn("Using np.mean and np.std to calculate the"
                "sigma limits, ignoring user defined functions..", AstropyUserWarning)

            filtered_data = sigmaclip(data, sigma_lower, sigma_upper)[0]
            return filtered_data
    elif HAS_SCIPY is False:
        if masked is False:
            warnings.warn("Scipy is not present - Cannot use the",
                    "scipy.stats.sigmaclip function. Falling back on ",
                    "astropy.stats.sigma_clip function.", AstropyUserWarning)
            # Will continue to use masked arrays for the rest of the code, 
            # no need to explicitly set masked to True since the only place it
            # is checked is here

    filtered_data = np.ma.array(data, copy=copy)

    if iters is None:
        i = -1
        lastrej = filtered_data.count() + 1
        while filtered_data.count() != lastrej:
            i += 1
            lastrej = filtered_data.count()
            perform_clip(filtered_data, kwargs)
    else:
        for i in range(iters):
            perform_clip(filtered_data, kwargs)

    # prevent filtered_data.mask = False (scalar) if no values are clipped
    if filtered_data.mask.shape == ():
        filtered_data.mask = False   # .mask shape will now match .data shape

    return filtered_data


sigma_clip.__doc__ = _sigma_clip.__doc__


def sigma_clipped_stats(data, mask=None, mask_value=None, sigma=3.0,
                        sigma_lower=None, sigma_upper=None, iters=5,
                        cenfunc=np.ma.median, stdfunc=np.std, axis=None,
                        masked=True):
    """
    Calculate sigma-clipped statistics from data.

    For example, sigma-clipped statistics can be used to estimate the
    background and background noise in an image.

    Parameters
    ----------
    data : array-like
        Data array or object that can be converted to an array.

    mask : `numpy.ndarray` (bool), optional
        A boolean mask with the same shape as ``data``, where a `True`
        value indicates the corresponding element of ``data`` is masked.
        Masked pixels are excluded when computing the image statistics.

    mask_value : float, optional
        An image data value (e.g., ``0.0``) that is ignored when
        computing the image statistics.  ``mask_value`` will be masked
        in addition to any input ``mask``.

    sigma : float, optional
        The number of standard deviations to use as the lower and upper
        clipping limit.  These limits are overridden by ``sigma_lower``
        and ``sigma_upper``, if input. Defaults to 3.

    sigma_lower : float, optional
        The number of standard deviations to use as the lower bound for
        the clipping limit.  If `None` then the value of ``sigma`` is used.
        Defaults to `None`.

    sigma_upper : float, optional
        The number of standard deviations to use as the upper bound for
        the clipping limit.  If `None` then the value of ``sigma`` is used.
        Defaults to `None`.

    iters : int, optional
        The number of iterations to perform sigma clipping, or `None` to
        clip until convergence is achieved (i.e., continue until the
        last iteration clips nothing) when calculating the image
        statistics. Defaults to 5.

    cenfunc : callable, optional
        The function used to compute the center for the clipping. Must
        be a callable that takes in a masked array and outputs the
        central value. Defaults to the median (`numpy.ma.median`).

    stdfunc : callable, optional
        The function used to compute the standard deviation about the
        center. Must be a callable that takes in a masked array and
        outputs a width estimator. Masked (rejected) pixels are those
        where::

             deviation < (-sigma_lower * stdfunc(deviation))
             deviation > (sigma_upper * stdfunc(deviation))

        where::

            deviation = data - cenfunc(data [,axis=int])

        Defaults to the standard deviation (`numpy.std`).

    axis : int or `None`, optional
        If not `None`, clip along the given axis.  For this case,
        ``axis`` will be passed on to ``cenfunc`` and ``stdfunc``, which
        are expected to return an array with the axis dimension removed
        (like the numpy functions).  If `None`, clip over all axes.
        Defaults to `None`.

    masked: bool, optional
        if `True` will use Numpy masked arrays to mask points rejected by
        sigma clipping. If `False` will use the `scipy.stats.sigmaclip` 
        function which throws away points rejected by sigma clipping. Setting

    Returns
    -------
    mean, median, stddev : float
        The mean, median, and standard deviation of the sigma-clipped
        image.
    """

    if mask is not None:
        data = np.ma.MaskedArray(data, mask)
    if mask_value is not None:
        data = np.ma.masked_values(data, mask_value)
    data_clip = sigma_clip(data, sigma=sigma, sigma_lower=sigma_lower,
                           sigma_upper=sigma_upper, iters=iters,
                           cenfunc=cenfunc, stdfunc=stdfunc, axis=axis,
                           masked=masked)
    goodvals = np.ma.compressed(data_clip)
    return np.mean(goodvals), np.median(goodvals), np.std(goodvals)
