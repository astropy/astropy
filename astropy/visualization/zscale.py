# Licensed under a 3-clause BSD style license - see LICENSE.rst

# Adapted from code from the Numdisplay package (see NUMDISPLAY_LICENSE.rst)

"""
Zscale implementation based on the one from the STScI Numdisplay package.

Original implementation from Numdisplay is BSD licensed:
https://trac.stsci.edu/ssb/stsci_python/browser/stsci_python/trunk/numdisplay/lib/stsci/numdisplay/zscale.py?rev=19347
https://trac.stsci.edu/ssb/stsci_python/browser/stsci_python/trunk/numdisplay/LICENSE.txt?rev=19347

"""

from __future__ import absolute_import, division

import numpy as np

__all__ = ['zscale']


def zscale(image, nsamples=1000, contrast=0.25, max_reject=0.5, min_npixels=5,
           krej=2.5, max_iterations=5):
    """Implement IRAF zscale algorithm.

    Parameters
    ----------
    image : array_like
        Input array.
    nsamples : int, optional
        Number of points in array to sample for determining scaling factors.
        Default to 1000.
    contrast : float, optional
        Scaling factor (between 0 and 1) for determining min and max. Larger
        values increase the difference between min and max values used for
        display. Default to 0.25.
    max_reject : float, optional
        If more than ``max_reject * npixels`` pixels are rejected, then the
        returned values are the min and max of the data. Default to 0.5.
    min_npixels : int, optional
        If less than ``min_npixels`` pixels are rejected, then the
        returned values are the min and max of the data. Default to 5.
    krej : float, optional
        Number of sigma used for the rejection. Default to 2.5.
    max_iterations : int, optional
        Maximum number of iterations for the rejection. Default to 5.

    Returns
    -------
    zmin, zmax: float
        Computed min and max values.

    """

    # Sample the image
    image = np.asarray(image)
    image = image[np.isfinite(image)]
    stride = int(max(1.0, image.size / nsamples))
    samples = image[::stride][:nsamples]
    samples.sort()

    npix = len(samples)
    zmin = samples[0]
    zmax = samples[-1]

    # Fit a line to the sorted array of samples
    minpix = max(min_npixels, int(npix * max_reject))
    x = np.arange(npix)
    ngoodpix = npix
    last_ngoodpix = npix + 1

    # Bad pixels mask used in k-sigma clipping
    badpix = np.zeros(npix, dtype=bool)

    # Kernel used to dilate the bad pixels mask
    ngrow = max(1, int(npix * 0.01))
    kernel = np.ones(ngrow, dtype=bool)

    for niter in range(max_iterations):
        if ngoodpix >= last_ngoodpix or ngoodpix < minpix:
            break

        fit = np.polyfit(x, samples, deg=1, w=(~badpix).astype(int))
        fitted = np.poly1d(fit)(x)

        # Subtract fitted line from the data array
        flat = samples - fitted

        # Compute the k-sigma rejection threshold
        threshold = krej * flat[~badpix].std()

        # Detect and reject pixels further than k*sigma from the fitted line
        badpix[(flat < - threshold) | (flat > threshold)] = True

        # Convolve with a kernel of length ngrow
        badpix = np.convolve(badpix, kernel, mode='same')

        last_ngoodpix = ngoodpix
        ngoodpix = np.sum(~badpix)

    slope, intercept = fit

    if ngoodpix >= minpix:
        if contrast > 0:
            slope = slope / contrast
        center_pixel = (npix - 1) // 2
        median = np.median(samples)
        zmin = max(zmin, median - (center_pixel - 1) * slope)
        zmax = min(zmax, median + (npix - center_pixel) * slope)
    return zmin, zmax
