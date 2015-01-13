# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np


__all__ = ['sigma_clip', 'sigma_clipped_stats']


def sigma_clip(data, sig=3, iters=1, cenfunc=np.ma.median, varfunc=np.var,
               axis=None, copy=True):
    """Perform sigma-clipping on the provided data.

    This performs the sigma clipping algorithm - i.e. the data will be iterated
    over, each time rejecting points that are more than a specified number of
    standard deviations discrepant.

    .. note::
        `scipy.stats.sigmaclip
        <http://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.sigmaclip.html>`_
        provides a subset of the functionality in this function.

    Parameters
    ----------
    data : array-like
        The data to be sigma-clipped (any shape).
    sig : float
        The number of standard deviations (*not* variances) to use as the
        clipping limit.
    iters : int or `None`
        The number of iterations to perform clipping for, or `None` to clip
        until convergence is achieved (i.e. continue until the last
        iteration clips nothing).
    cenfunc : callable
        The technique to compute the center for the clipping. Must be a
        callable that takes in a masked array and outputs the central value.
        Defaults to the median (numpy.median).
    varfunc : callable
        The technique to compute the standard deviation about the center. Must
        be a callable that takes in a masked array and outputs a width
        estimator::

             deviation**2 > sig**2 * varfunc(deviation)

        Defaults to the variance (numpy.var).

    axis : int or `None`
        If not `None`, clip along the given axis.  For this case, axis=int will
        be passed on to cenfunc and varfunc, which are expected to return an
        array with the axis dimension removed (like the numpy functions).
        If `None`, clip over all values.  Defaults to `None`.
    copy : bool
        If `True`, the data array will be copied.  If `False`, the masked array
        data will contain the same array as ``data``.  Defaults to `True`.

    Returns
    -------
    filtered_data : `numpy.ma.MaskedArray`
        A masked array with the same shape as ``data`` input, where the points
        rejected by the algorithm have been masked.

    Notes
    -----
     1. The routine works by calculating::

            deviation = data - cenfunc(data [,axis=int])

        and then setting a mask for points outside the range::

            data.mask = deviation**2 > sig**2 * varfunc(deviation)

        It will iterate a given number of times, or until no further points are
        rejected.

     2. Most numpy functions deal well with masked arrays, but if one would
        like to have an array with just the good (or bad) values, one can use::

            good_only = filtered_data.data[~filtered_data.mask]
            bad_only = filtered_data.data[filtered_data.mask]

        However, for multidimensional data, this flattens the array, which may
        not be what one wants (especially is filtering was done along an axis).

    Examples
    --------

    This will generate random variates from a Gaussian distribution and return
    a masked array in which all points that are more than 2 *sample* standard
    deviation from the median are masked::

        >>> from astropy.stats import sigma_clip
        >>> from numpy.random import randn
        >>> randvar = randn(10000)
        >>> filtered_data = sigma_clip(randvar, 2, 1)

    This will clipping on a similar distribution, but for 3 sigma relative to
    the sample *mean*, will clip until converged, and does not copy the data::

        >>> from astropy.stats import sigma_clip
        >>> from numpy.random import randn
        >>> from numpy import mean
        >>> randvar = randn(10000)
        >>> filtered_data = sigma_clip(randvar, 3, None, mean, copy=False)

    This will clip along one axis on a similar distribution with bad points
    inserted::

        >>> from astropy.stats import sigma_clip
        >>> from numpy.random import normal
        >>> from numpy import arange, diag, ones
        >>> data = arange(5)+normal(0.,0.05,(5,5))+diag(ones(5))
        >>> filtered_data = sigma_clip(data, axis=0, sig=2.3)

    Note that along the other axis, no points would be masked, as the variance
    is higher.

    """

    if axis is not None:
        cenfunc_in = cenfunc
        varfunc_in = varfunc
        cenfunc = lambda d: np.expand_dims(cenfunc_in(d, axis=axis), axis=axis)
        varfunc = lambda d: np.expand_dims(varfunc_in(d, axis=axis), axis=axis)

    filtered_data = np.ma.array(data, copy=copy)

    if iters is None:
        i = -1
        lastrej = filtered_data.count() + 1
        while filtered_data.count() != lastrej:
            i += 1
            lastrej = filtered_data.count()
            do = filtered_data - cenfunc(filtered_data)
            filtered_data.mask |= do * do > varfunc(filtered_data) * sig ** 2
    else:
        for i in range(iters):
            do = filtered_data - cenfunc(filtered_data)
            filtered_data.mask |= do * do > varfunc(filtered_data) * sig ** 2

    return filtered_data


def sigma_clipped_stats(data, mask=None, mask_val=None, sigma=3.0, iters=None):
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

    mask_val : float, optional
        An image data value (e.g., ``0.0``) that is ignored when
        computing the image statistics.  ``mask_val`` will be masked in
        addition to any input ``mask``.

    sigma : float, optional
        The number of standard deviations to use as the clipping limit.

    iters : int, optional
        The number of iterations to perform sigma clipping, or `None` to
        clip until convergence is achieved (i.e., continue until the
        last iteration clips nothing) when calculating the image
        statistics.

    Returns
    -------
    mean, median, stddev : float
        The mean, median, and standard deviation of the sigma-clipped
        image.
    """

    if mask is not None:
        data = np.ma.MaskedArray(data, mask)
    if mask_val is not None:
        data = np.ma.masked_values(data, mask_val)
    data_clip = sigma_clip(data, sig=sigma, iters=iters)
    goodvals = data_clip.data[~data_clip.mask]
    return np.mean(goodvals), np.median(goodvals), np.std(goodvals)
