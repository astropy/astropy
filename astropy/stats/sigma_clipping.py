# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import numpy as np
from copy import deepcopy
from collections import OrderedDict
from astropy.utils.compat.funcsigs import Parameter, Signature


__all__ = ['sigma_clip', 'sigma_clipped_stats']


_sigma_clip_keywords = OrderedDict([('sigma', 3.), ('sigma_lower', None),
                                    ('sigma_upper', None), ('iters', 1),
                                    ('cenfunc', np.ma.median),
                                    ('stdfunc', np.std), ('axis', None),
                                    ('copy', True)])
_sigma_clip_keywords_old = deepcopy(_sigma_clip_keywords)
_sigma_clip_keywords_old['sig'] = 3.
_sigma_clip_keywords_old['varfunc'] = np.var


def _make_sigma_clip_signature(_sigma_clip_keywords):
    params = ([Parameter('data', Parameter.POSITIONAL_ONLY)] +
              [Parameter(name, Parameter.KEYWORD_ONLY, default=default)
               for name, default in _sigma_clip_keywords.items()])
    return Signature(params)


def sigma_clip(*args, **kwargs):
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
    sigma : float, optional
        The number of standard deviations to use as the lower and upper
        clipping limit.  These limits are overridden by ``sigma_lower``
        and ``sigma_upper``, if input.
    sigma_lower : float, optional
        The number of standard deviations to use as the lower bound for
        the clipping limit.  If `None` then the value of ``sigma`` is used.
    sigma_upper : float, optional
        The number of standard deviations to use as the upper bound for
        the clipping limit.  If `None` then the value of ``sigma`` is used.
    iters : int or `None`, optional
        The number of iterations to perform clipping for, or `None` to clip
        until convergence is achieved (i.e. continue until the last
        iteration clips nothing).
    cenfunc : callable, optional
        The technique to compute the center for the clipping. Must be a
        callable that takes in a masked array and outputs the central value.
        Defaults to the median (`numpy.ma.median`).
    stdfunc : callable, optional
        The technique to compute the standard deviation about the center. Must
        be a callable that takes in a masked array and outputs a width
        estimator.  Masked (rejected) pixels are those where::

             deviation < (-sigma_lower * stdfunc(deviation))
             deviation > (sigma_upper * stdfunc(deviation))

        Defaults to the standard deviation (`numpy.std`).

    axis : int or `None`, optional
        If not `None`, clip along the given axis.  For this case, axis=int will
        be passed on to cenfunc and stdfunc, which are expected to return an
        array with the axis dimension removed (like the numpy functions).
        If `None`, clip over all values.  Defaults to `None`.
    copy : bool, optional
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

           deviation < (-sigma_lower * stdfunc(deviation))
           deviation > (sigma_upper * stdfunc(deviation))

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
        >>> filtered_data = sigma_clip(randvar, 2, iters=1)

    This will clipping on a similar distribution, but for 3 sigma relative to
    the sample *mean*, will clip until converged, and does not copy the data::

        >>> from astropy.stats import sigma_clip
        >>> from numpy.random import randn
        >>> from numpy import mean
        >>> randvar = randn(10000)
        >>> filtered_data = sigma_clip(randvar, 3, iters=None, cenfunc=mean,
        ...                            copy=False)

    This will clip along one axis on a similar distribution with bad points
    inserted::

        >>> from astropy.stats import sigma_clip
        >>> from numpy.random import normal
        >>> from numpy import arange, diag, ones
        >>> data = arange(5)+normal(0.,0.05,(5,5))+diag(ones(5))
        >>> filtered_data = sigma_clip(data, axis=0, sigma=2.3)

    Note that along the other axis, no points would be masked, as the variance
    is higher.
    """

    signature = _make_sigma_clip_signature(_sigma_clip_keywords_old)
    bound_args = signature.bind(*args, **kwargs)
    for name, value in bound_args.arguments.items():
        print(name, value)

    for param in signature.parameters.values():
        if (param.name not in bound_args.arguments and
            param.default is not param.empty):
                bound_args.arguments[param.name] = param.default

    for name, value in bound_args.arguments.items():
        print(name, value)

    if sigma_lower is None:
        sigma_lower = sigma
    if sigma_upper is None:
        sigma_upper = sigma

    if axis is not None:
        cenfunc_in = cenfunc
        stdfunc_in = stdfunc
        cenfunc = lambda d: np.expand_dims(cenfunc_in(d, axis=axis), axis=axis)
        stdfunc = lambda d: np.expand_dims(stdfunc_in(d, axis=axis), axis=axis)

    filtered_data = np.ma.array(data, copy=copy)

    if iters is None:
        i = -1
        lastrej = filtered_data.count() + 1
        while filtered_data.count() != lastrej:
            i += 1
            lastrej = filtered_data.count()
            deviation = filtered_data - cenfunc(filtered_data)
            std = stdfunc(filtered_data)
            filtered_data.mask |= np.ma.masked_less(deviation,
                                                    -std * sigma_lower).mask
            filtered_data.mask |= np.ma.masked_greater(deviation,
                                                       std * sigma_upper).mask
    else:
        for i in range(iters):
            deviation = filtered_data - cenfunc(filtered_data)
            std = stdfunc(filtered_data)
            filtered_data.mask |= np.ma.masked_less(deviation,
                                                    -std * sigma_lower).mask
            filtered_data.mask |= np.ma.masked_greater(deviation,
                                                       std * sigma_upper).mask

    return filtered_data


sigma_clip.__signature__ = _make_sigma_clip_signature(_sigma_clip_keywords)


def sigma_clipped_stats(data, mask=None, mask_value=None, sigma=3.0,
                        sigma_lower=3.0, sigma_upper=3.0, iters=None):
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
        and ``sigma_upper``, if input.

    sigma_lower : float, optional
        The number of standard deviations to use as the lower bound for
        the clipping limit.  If `None` then the value of ``sigma`` is used.

    sigma_upper : float, optional
        The number of standard deviations to use as the upper bound for
        the clipping limit.  If `None` then the value of ``sigma`` is used.

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
    if mask_value is not None:
        data = np.ma.masked_values(data, mask_value)
    data_clip = sigma_clip(data, sigma=sigma, sigma_lower=sigma_lower,
                           sigma_upper=sigma_upper, iters=iters)
    goodvals = data_clip.data[~data_clip.mask]
    return np.mean(goodvals), np.median(goodvals), np.std(goodvals)
