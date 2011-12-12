# Licensed under a 3-clause BSD style license - see LICENSE.rst
#This module implements the base NDData class.

__all__ = ['NDData']


class NDData(object):
    """Class to store array type data.
    `NDData` provides a superclass for all array based data.

    Initialize an NDData object with data
    Parameters
    ----------
    data : `~numpy.ndarray`
        n-dimensional array containing the n-dimensional data
    error: `~numpy.ndarray`, optional
        n-dimensional array containing the error of the data
    mask: `~numpy.ndarray`, optional
        n-dimensional array masking the data (suggested dtype=`bool`)
    wcs: `astropy.wcs`, optional
        WCS-object containing the world coordinate algorithms for the adta
    meta: `dict`-like object, optional
        contains meta data for the n-dimensional data
    units: `astropy.units`, optional
        describing the units of the data

    """
    def __init__(self, data, error=None, mask=None, wcs=None, meta=None,
                 units=None):
        self.data = data
        self.error = error
        self.mask = mask
        self.wcs = wcs
        self.meta = meta
        self.units = units
