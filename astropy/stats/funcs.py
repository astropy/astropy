# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module contains simple statistical algorithms that are straightforwardly
implemented as a single python function (or family of functions).

This package should generally not be used directly.  Everything in `__all__` is
imported into `astropy.stats`, and hence that package should be used for
access.
"""

import numpy as np  # needed for some defaults

__all__ = ['sigma_clip']


def sigma_clip(data, sig=3, iters=1, cenfunc=np.median, varfunc=np.var,
               maout=False):
    """ Perform sigma-clipping on the provided data.

    This performs the sigma clipping algorithm - i.e. the data will be iterated
    over, each time rejecting points that are more than a specified number of
    standard deviations discrepant.

    .. note::
        `scipy.stats.sigmaclip` provides a subset of the functionality in this
        function.

    Parameters
    ----------
    data : array-like
        The data to be sigma-clipped (any shape).
    sig : float
        The number of standard deviations (*not* variances) to use as the
        clipping limit.
    iters : int or None
        The number of iterations to perform clipping for, or None to clip until
        convergence is achieved (i.e. continue until the last iteration clips
        nothing).
    cenfunc : callable
        The technique to compute the center for the clipping. Must be a
        callable that takes in a 1D data array and outputs the central value.
        Defaults to the median.
    varfunc : callable
        The technique to compute the variance about the center. Must be a
        callable that takes in a 1D data array and outputs the width estimator
        that will be interpreted as a variance. Defaults to the variance.
    maout : bool or 'copy'
        If True, a masked array will be returned. If the special string
        'inplace', the masked array will contain the same array as `data`,
        otherwise the array data will be copied.

    Returns
    -------
    filtereddata : `numpy.ndarray` or `numpy.masked.MaskedArray`
        If `maout` is True, this is a masked array with a shape matching the
        input that is masked where the algorithm has rejected those values.
        Otherwise, a 1D array of values including only those that are not
        clipped.
    mask : boolean array
        Only present if `maout` is False. A boolean array with a shape matching
        the input `data` that is False for rejected values and True for all
        others.

    Examples
    --------

    This will generate random variates from a Gaussian distribution and return
    only the points that are within 2 *sample* standard deviation from the
    median::

        >>> from astropy.stats import sigma_clip
        >>> from numpy.random import randn
        >>> randvar = randn(10000)
        >>> data,mask = sigma_clip(randvar, 2, 1)

    This will clipping on a similar distribution, but for 3 sigma relative to
    the sample *mean*, will clip until converged, and produces a
    `numpy.masked.MaskedArray`::

        >>> from astropy.stats import sigma_clip
        >>> from numpy.random import randn
        >>> from numpy import mean
        >>> randvar = randn(10000)
        >>> maskedarr = sigma_clip(randvar, 3, None, mean, maout=True)

    """

    data = np.array(data, copy=False)
    oldshape = data.shape
    data = data.ravel()

    mask = np.ones(data.size, bool)
    if iters is None:
        i = -1
        lastrej = sum(mask) + 1
        while(sum(mask) != lastrej):
            i += 1
            lastrej = sum(mask)
            do = data - cenfunc(data[mask])
            mask = do * do <= varfunc(data[mask]) * sig ** 2
        iters = i + 1
        #TODO: ?print iters to the log if iters was None?
    else:
        for i in range(iters):
            do = data - cenfunc(data[mask])
            mask = do * do <= varfunc(data[mask]) * sig ** 2

    if maout:
        return np.ma.MaskedArray(data, ~mask, copy=maout != 'inplace')
    else:
        return data[mask], mask.reshape(oldshape)
