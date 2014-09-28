# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import numpy as np
from astropy.stats import sigma_clip

__all__ = ['sigmaclip_stats']


def sigmaclip_stats(image, image_mask=None, mask_val=None, sigma=3.0,
                    iters=None):
    """
    Calculate sigma-clipped statistics of an image.  For example,
    sigma-clipped statistics can be used to estimate the background and
    background noise in an image.

    Parameters
    ----------
    image : array_like
        The 2D array of the image.

    image_mask : array_like, bool, optional
        A boolean mask with the same shape as ``image``, where a `True`
        value indicates the corresponding element of ``image`` is
        invalid.  Masked pixels are ignored when computing the image
        statistics.

    mask_val : float, optional
        An image data value (e.g., ``0.0``) that is ignored when
        computing the image statistics.  ``mask_val`` will be ignored if
        ``image_mask`` is input.

    sigma : float, optional
        The number of standard deviations to use as the clipping limit.

    iters : float, optional
       The number of clipping iterations to perform, or `None` to clip
       until convergence is achieved (i.e. continue until the last
       iteration clips nothing).

    Returns
    -------
    mean, median, stddev : floats
        The mean, median, and standard deviation of the sigma-clipped
        image.
    """

    if image_mask is not None:
        assert image_mask.dtype == np.bool, \
            'image_mask must be a boolean ndarray'
        image = image[~image_mask]
    if mask_val is not None and image_mask is None:
        idx = (image != mask_val).nonzero()
        image = image[idx]
    image_clip = sigma_clip(image, sig=sigma, iters=iters)
    goodvals = image_clip.data[~image_clip.mask]
    return np.mean(goodvals), np.median(goodvals), np.std(goodvals)
