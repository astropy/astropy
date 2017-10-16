# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import numpy as np
import warnings
from ..utils.exceptions import AstropyUserWarning
from ..extern.six.moves import range


__all__ = ['SigmaClip', 'sigma_clip', 'sigma_clipped_stats']


class SigmaClip(object):
    """
    Class to perform sigma clipping.

    The data will be iterated over, each time rejecting points that are
    discrepant by more than a specified number of standard deviations
    from a center value. If the data contains invalid values (NaNs or
    infs), they are automatically masked before performing the sigma
    clipping.

    For a functional interface to sigma clipping, see
    :func:`sigma_clip`.

    .. note::
        `scipy.stats.sigmaclip
        <https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.sigmaclip.html>`_
        provides a subset of the functionality in this class.

    Parameters
    ----------
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

    See Also
    --------
    sigma_clip

    Examples
    --------
    This example generates random variates from a Gaussian distribution
    and returns a masked array in which all points that are more than 2
    sample standard deviations from the median are masked::

        >>> from astropy.stats import SigmaClip
        >>> from numpy.random import randn
        >>> randvar = randn(10000)
        >>> sigclip = SigmaClip(sigma=2, iters=5)
        >>> filtered_data = sigclip(randvar)

    This example sigma clips on a similar distribution, but uses 3 sigma
    relative to the sample *mean*, clips until convergence, and does not
    copy the data::

        >>> from astropy.stats import SigmaClip
        >>> from numpy.random import randn
        >>> from numpy import mean
        >>> randvar = randn(10000)
        >>> sigclip = SigmaClip(sigma=3, iters=None, cenfunc=mean)
        >>> filtered_data = sigclip(randvar, copy=False)

    This example sigma clips along one axis on a similar distribution
    (with bad points inserted)::

        >>> from astropy.stats import SigmaClip
        >>> from numpy.random import normal
        >>> from numpy import arange, diag, ones
        >>> data = arange(5) + normal(0., 0.05, (5, 5)) + diag(ones(5))
        >>> sigclip = SigmaClip(sigma=2.3)
        >>> filtered_data = sigclip(data, axis=0)

    Note that along the other axis, no points would be masked, as the
    variance is higher.
    """

    def __init__(self, sigma=3., sigma_lower=None, sigma_upper=None, iters=5,
                 cenfunc=np.ma.median, stdfunc=np.std):
        self.sigma = sigma
        self.sigma_lower = sigma_lower
        self.sigma_upper = sigma_upper
        self.iters = iters
        self.cenfunc = cenfunc
        self.stdfunc = stdfunc

    def __repr__(self):
        return ('SigmaClip(sigma={0}, sigma_lower={1}, sigma_upper={2}, '
                'iters={3}, cenfunc={4}, stdfunc={5})'
                .format(self.sigma, self.sigma_lower, self.sigma_upper,
                        self.iters, self.cenfunc, self.stdfunc))

    def __str__(self):
        lines = ['<' + self.__class__.__name__ + '>']
        attrs = ['sigma', 'sigma_lower', 'sigma_upper', 'iters', 'cenfunc',
                 'stdfunc']
        for attr in attrs:
            lines.append('    {0}: {1}'.format(attr, getattr(self, attr)))
        return '\n'.join(lines)

    def _perform_clip(self, _filtered_data, axis=None):
        """
        Perform sigma clip by comparing the data to the minimum and
        maximum values (median + sig * standard deviation). Use
        sigma_lower and sigma_upper to get the correct limits. Data
        values less or greater than the minimum / maximum values
        will have True set in the mask array.
        """

        if _filtered_data.size == 0:
            return _filtered_data

        max_value = self.cenfunc(_filtered_data, axis=axis)
        std = self.stdfunc(_filtered_data, axis=axis)
        min_value = max_value - std * self.sigma_lower
        max_value += std * self.sigma_upper

        if axis is not None:
            if axis != 0:
                min_value = np.expand_dims(min_value, axis=axis)
                max_value = np.expand_dims(max_value, axis=axis)
        if max_value is np.ma.masked:
            max_value = np.ma.MaskedArray(np.nan, mask=True)
            min_value = np.ma.MaskedArray(np.nan, mask=True)

        _filtered_data.mask |= _filtered_data > max_value
        _filtered_data.mask |= _filtered_data < min_value

        return _filtered_data

    def __call__(self, data, axis=None, copy=True):
        """
        Perform sigma clipping on the provided data.

        Parameters
        ----------
        data : array-like
            The data to be sigma clipped.
        axis : int or `None`, optional
            If not `None`, clip along the given axis.  For this case,
            ``axis`` will be passed on to ``cenfunc`` and ``stdfunc``,
            which are expected to return an array with the axis
            dimension removed (like the numpy functions).  If `None`,
            clip over all axes.  Defaults to `None`.
        copy : bool, optional
            If `True`, the ``data`` array will be copied.  If `False`,
            the returned masked array data will contain the same array
            as ``data``.  Defaults to `True`.

        Returns
        -------
        filtered_data : `numpy.ma.MaskedArray`
            A masked array with the same shape as ``data`` input, where
            the points rejected by the algorithm have been masked.
        """

        if self.sigma_lower is None:
            self.sigma_lower = self.sigma
        if self.sigma_upper is None:
            self.sigma_upper = self.sigma

        if np.any(~np.isfinite(data)):
            data = np.ma.masked_invalid(data)
            warnings.warn('Input data contains invalid values (NaNs or '
                          'infs), which were automatically masked.',
                          AstropyUserWarning)

        filtered_data = np.ma.array(data, copy=copy)

        if self.iters is None:
            lastrej = filtered_data.count() + 1
            while filtered_data.count() != lastrej:
                lastrej = filtered_data.count()
                self._perform_clip(filtered_data, axis=axis)
        else:
            for i in range(self.iters):
                self._perform_clip(filtered_data, axis=axis)

        # prevent filtered_data.mask = False (scalar) if no values are clipped
        if filtered_data.mask.shape == ():
            # make .mask shape match .data shape
            filtered_data.mask = False

        return filtered_data


def sigma_clip(data, sigma=3, sigma_lower=None, sigma_upper=None, iters=5,
               cenfunc=np.ma.median, stdfunc=np.std, axis=None, copy=True):
    """
    Perform sigma-clipping on the provided data.

    The data will be iterated over, each time rejecting points that are
    discrepant by more than a specified number of standard deviations from a
    center value. If the data contains invalid values (NaNs or infs),
    they are automatically masked before performing the sigma clipping.

    For an object-oriented interface to sigma clipping, see
    :func:`SigmaClip`.

    .. note::
        `scipy.stats.sigmaclip
        <https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.sigmaclip.html>`_
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

    See Also
    --------
    SigmaClip

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

    sigclip = SigmaClip(sigma=sigma, sigma_lower=sigma_lower,
                        sigma_upper=sigma_upper, iters=iters,
                        cenfunc=cenfunc, stdfunc=stdfunc)
    return sigclip(data, axis=axis, copy=copy)


def sigma_clipped_stats(data, mask=None, mask_value=None, sigma=3.0,
                        sigma_lower=None, sigma_upper=None, iters=5,
                        cenfunc=np.ma.median, stdfunc=np.std, std_ddof=0,
                        axis=None):
    """
    Calculate sigma-clipped statistics on the provided data.

    Parameters
    ----------
    data : array-like
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
        last iteration clips nothing) when calculating the statistics.
        Defaults to 5.

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

    std_ddof : int, optional
        The delta degrees of freedom for the standard deviation
        calculation.  The divisor used in the calculation is ``N -
        std_ddof``, where ``N`` represents the number of elements.  The
        default is zero.

    axis : int or `None`, optional
        If not `None`, clip along the given axis.  For this case,
        ``axis`` will be passed on to ``cenfunc`` and ``stdfunc``, which
        are expected to return an array with the axis dimension removed
        (like the numpy functions).  If `None`, clip over all axes.
        Defaults to `None`.

    Returns
    -------
    mean, median, stddev : float
        The mean, median, and standard deviation of the sigma-clipped
        data.
    """

    if mask is not None:
        data = np.ma.MaskedArray(data, mask)
    if mask_value is not None:
        data = np.ma.masked_values(data, mask_value)

    data_clip = sigma_clip(data, sigma=sigma, sigma_lower=sigma_lower,
                           sigma_upper=sigma_upper, iters=iters,
                           cenfunc=cenfunc, stdfunc=stdfunc, axis=axis)

    mean = np.ma.mean(data_clip, axis=axis)
    median = np.ma.median(data_clip, axis=axis)
    std = np.ma.std(data_clip, ddof=std_ddof, axis=axis)

    if axis is None and np.ma.isMaskedArray(median):
        # With Numpy 1.10 np.ma.median always return a MaskedArray, even with
        # one element. So for compatibility with previous versions, we take the
        # scalar value
        median = median.item()

    return mean, median, std
