# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import numpy as np
from astropy.stats import sigma_clip

__all__ = ['sigmaclip_stats']


def sigmaclip_stats(data, mask=None, mask_val=None, sigma=3.0, iters=None):
    """
    Calculate sigma-clipped statistics of an image.

    For example, sigma-clipped statistics can be used to estimate the
    background and background noise in an image.

    Parameters
    ----------
    data : array_like
        The 2D array of the image.

    mask : `numpy.ndarray` (bool), optional
        A boolean mask with the same shape as ``data``, where a `True`
        value indicates the corresponding element of ``data`` is masked.
        Masked pixels are excluded when computing the image statistics.

    mask_val : float, optional
        An image data value (e.g., ``0.0``) that is ignored when
        computing the image statistics.  ``mask_val`` will be ignored if
        ``mask`` is input.

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
        if mask.dtype != np.bool:
            raise TypeError('mask must be a boolean ndarray')
        data = data[~mask]
    if mask_val is not None and mask is None:
        idx = (data != mask_val).nonzero()
        data = data[idx]
    data_clip = sigma_clip(data, sig=sigma, iters=iters)
    goodvals = data_clip.data[~data_clip.mask]
    return np.mean(goodvals), np.median(goodvals), np.std(goodvals)
