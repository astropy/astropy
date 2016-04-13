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
    ngrow = max(1, int(npix * 0.01))
    ngoodpix, zstart, zslope = zsc_fit_line(samples, npix, KREJ, ngrow,
                                            MAX_ITERATIONS)

    if ngoodpix >= minpix:
        if contrast > 0:
            zslope = zslope / contrast
        center_pixel = (npix - 1) // 2
        median = np.median(samples)
        zmin = max(zmin, median - (center_pixel - 1) * zslope)
        zmax = min(zmax, median + (npix - center_pixel) * zslope)
    return zmin, zmax


def zsc_fit_line(samples, npix, krej, ngrow, maxiter):
    # First re-map indices from -1.0 to 1.0
    xscale = 2.0 / (npix - 1)
    xnorm = np.linspace(-1, 1, npix)

    ngoodpix = npix
    minpix = max(MIN_NPIXELS, int(npix * MAX_REJECT))
    last_ngoodpix = npix + 1

    kernel = np.ones(ngrow, dtype=bool)

    # This is the mask used in k-sigma clipping
    badpix = np.zeros(npix, dtype=bool)

    for niter in range(maxiter):

        if (ngoodpix >= last_ngoodpix) or (ngoodpix < minpix):
            break

        # Accumulate sums to calculate straight line fit
        goodpixels = ~badpix
        sumx = xnorm[goodpixels].sum()
        sumxx = (xnorm[goodpixels] * xnorm[goodpixels]).sum()
        sumxy = (xnorm[goodpixels] * samples[goodpixels]).sum()
        sumy = samples[goodpixels].sum()
        sum = goodpixels.sum()

        delta = sum * sumxx - sumx * sumx
        # Slope and intercept
        intercept = (sumxx * sumy - sumx * sumxy) / delta
        slope = (sum * sumxy - sumx * sumy) / delta

        # Subtract fitted line from the data array
        fitted = xnorm * slope + intercept
        flat = samples - fitted

        # Compute the k-sigma rejection threshold
        threshold = krej * flat[goodpixels].std()

        # Detect and reject pixels further than k*sigma from the fitted line
        badpix[flat < - threshold] = True
        badpix[flat > threshold] = True

        # Convolve with a kernel of length ngrow
        badpix = np.convolve(badpix, kernel, mode='same')

        last_ngoodpix = ngoodpix
        ngoodpix = np.sum(~badpix)

    # Transform the line coefficients back to the X range [0:npix-1]
    zstart = intercept - slope
    zslope = slope * xscale

    return ngoodpix, zstart, zslope
