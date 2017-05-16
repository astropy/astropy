# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np

def azimuthalAverage(image, center=None, stddev=False, returnradii=False,
                     return_nr=False, binsize=0.5, weights=None, steps=False,
                     interpnan=False, left=None, right=None, mask=None):
    """
    Calculate the azimuthally averaged radial profile.

    Parameters
    ----------
    image : `~numpy.ndarray`
        The 2D image
    center : (int, int)
        The [x,y] pixel coordinates used as the center. The default is ``None``,
        which then uses the center of the image (including fractional pixels).
    stddev : bool
        if specified, return the azimuthal standard deviation instead of the
        average
    returnradii : bool
        if specified, return (radii_array,radial_profile)
    return_nr : bool
        if specified, return number of pixels per radius *and* radius
    binsize : bool
        size of the averaging bin.  Can lead to strange results if non-binsize
        factors are used to specify the center and the binsize is too large
    weights : `~numpy.ndarray` or ``None``
        can do a weighted average instead of a simple average if this keyword
        parameter is set.  ``weights.shape`` must = ``image.shape``.  weighted
        stddev is undefined, so don't set weights and stddev.
    steps : bool
        if specified, will return a double-length bin array and radial profile
        so you can plot a step-form radial profile (which more accurately
        represents what's going on)
    interpnan : bool
        Interpolate over NAN values, i.e. bins where there is no data?
    left,right : float,float
        passed to interpnan; they set the extrapolated values
    mask : `~numpy.ndarray`
        can supply a mask (boolean array same size as image with True for OK
        and False for not) to average over only select data.

    Returns
    -------
    radial_profile : `~numpy.ndarray`
        The radial profile.
        If a bin contains NO DATA, it will have a NAN value because of the
        divide-by-sum-of-weights component.

    """
    # Calculate the indices from the image
    yy, xx = np.indices(image.shape)

    if center is None:
        center = np.array([(xx.max()-xx.min())/2.0, (yy.max()-yy.min())/2.0])

    rr = np.hypot(xx - center[0], yy - center[1])

    if weights is None:
        weights = np.ones(image.shape)
    elif stddev:
        raise ValueError("Weighted standard deviation is not defined.")

    if mask is None:
        mask = np.ones(image.shape, dtype='bool')
    # obsolete elif len(mask.shape) > 1:
    # obsolete     mask = mask.ravel()

    # the 'bins' as initially defined are lower/upper bounds for each bin
    # so that values will be in [lower,upper)
    nbins = int(np.round(rr.max() / binsize)+1)
    maxbin = nbins * binsize
    bins = np.linspace(0, maxbin, nbins+1)
    # but we're probably more interested in the bin centers than their left or right sides...
    bin_centers = (bins[1:]+bins[:-1])/2.0

    # how many per bin (i.e., histogram)?
    # there are never any in bin 0, because the lowest index returned by digitize is 1
    #nr = np.bincount(whichbin)[1:]
    nr = np.histogram(rr, bins, weights=mask.astype('int'))[0]

    # recall that bins are from 1 to nbins (which is expressed in array terms by arange(nbins)+1 or xrange(1,nbins+1) )
    # radial_prof.shape = bin_centers.shape
    if stddev:
        # Find out which radial bin each point in the map belongs to
        whichbin = np.digitize(rr.flat, bins)
        # This method is still very slow; is there a trick to do this with histograms?
        radial_prof = np.array([image.flat[mask.flat*(whichbin==b)].std()
                                for b in xrange(1, nbins+1)])
    else:
        radial_prof = np.histogram(rr, bins, weights=(image*weights*mask))[0] / np.histogram(rr, bins, weights=(mask*weights))[0]

    if interpnan:
        radial_prof = np.interp(bin_centers,
                                bin_centers[radial_prof==radial_prof],
                                radial_prof[radial_prof==radial_prof],
                                left=left, right=right)

    if steps:
        xarr = np.array(zip(bins[:-1], bins[1:])).ravel()
        yarr = np.array(zip(radial_prof, radial_prof)).ravel()
        return xarr, yarr
    elif returnradii:
        return bin_centers, radial_prof
    elif return_nr:
        return nr, bin_centers, radial_prof
    else:
        return radial_prof

def azimuthalAverageBins(image, azbins, symmetric=None,  center=None,
                         **kwargs):
    """ Compute the azimuthal average over a limited range of angles
    kwargs are passed to azimuthalAverage """
    y, x = np.indices(image.shape)
    if center is None:
        center = np.array([(x.max()-x.min())/2.0, (y.max()-y.min())/2.0])
    r = np.hypot(x - center[0], y - center[1])
    theta = np.arctan2(x - center[0], y - center[1])
    theta[theta < 0] += 2*np.pi
    theta_deg = theta*180.0/np.pi

    if isinstance(azbins, np.ndarray):
        pass
    elif isinstance(azbins, int):
        if symmetric == 2:
            azbins = np.linspace(0, 90, azbins)
            theta_deg = theta_deg % 90
        elif symmetric == 1:
            azbins = np.linspace(0, 180, azbins)
            theta_deg = theta_deg % 180
        elif azbins == 1:
            return azbins, azimuthalAverage(image, center=center,
                                            returnradii=True, **kwargs)
        else:
            azbins = np.linspace(0, 359.9999999999999, azbins)
    else:
        raise ValueError("azbins must be an ndarray or an integer")

    azavlist = []
    for blow, bhigh in zip(azbins[:-1], azbins[1:]):
        mask = (theta_deg > (blow % 360)) * (theta_deg < (bhigh % 360))
        rr, zz = azimuthalAverage(image, center=center, mask=mask,
                                  returnradii=True, **kwargs)
        azavlist.append(zz)

    return azbins, rr, azavlist
