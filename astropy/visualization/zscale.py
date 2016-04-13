# Licensed under a 3-clause BSD style license - see NUMDISPLAY_LICENSE.rst

"""
Zscale implementation from the STScI numdisplay package, BSD licensed.

https://trac.stsci.edu/ssb/stsci_python/browser/stsci_python/trunk/numdisplay/lib/stsci/numdisplay/zscale.py?rev=19347
https://trac.stsci.edu/ssb/stsci_python/browser/stsci_python/trunk/numdisplay/LICENSE.txt?rev=19347

"""

from __future__ import absolute_import, division

import numpy as np

__all__ = ['zscale']

MAX_REJECT = 0.5
MIN_NPIXELS = 5
KREJ = 2.5
MAX_ITERATIONS = 5


def zscale(image, nsamples=1000, contrast=0.25):
    """Implement IRAF zscale algorithm

    Parameters
    ----------
    image : arr
        2-d numpy array

    nsamples : int (Default: 1000)
        Number of points in array to sample for determining scaling factors

    contrast : float (Default: 0.25)
        Scaling factor for determining min and max. Larger values increase the
        difference between min and max values used for display.

    Returns
    -------
    (z1, z2)
    """

    # Sample the image
    stride = int(image.size / nsamples)
    samples = image.flatten()[::stride][:nsamples]
    samples.sort()

    npix = len(samples)
    zmin = samples[0]
    zmax = samples[-1]

    # Fit a line to the sorted array of samples
    minpix = max(MIN_NPIXELS, int(npix * MAX_REJECT))
    ngoodpix, zslope = zsc_fit_line(samples, KREJ, MAX_ITERATIONS)

    if ngoodpix >= minpix:
        if contrast > 0:
            zslope = zslope / contrast
        center_pixel = (npix - 1) // 2
        median = np.median(samples)
        zmin = max(zmin, median - (center_pixel - 1) * zslope)
        zmax = min(zmax, median + (npix - center_pixel) * zslope)
    return zmin, zmax


def zsc_fit_line(samples, krej, maxiter):
    npix = len(samples)
    x = np.arange(npix)

    ngoodpix = npix
    minpix = max(MIN_NPIXELS, int(npix * MAX_REJECT))
    last_ngoodpix = npix + 1

    ngrow = max(1, int(npix * 0.01))
    kernel = np.ones(ngrow, dtype=bool)

    # This is the mask used in k-sigma clipping
    badpix = np.zeros(npix, dtype=bool)

    for niter in range(maxiter):
        if (ngoodpix >= last_ngoodpix) or (ngoodpix < minpix):
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
    return ngoodpix, slope
