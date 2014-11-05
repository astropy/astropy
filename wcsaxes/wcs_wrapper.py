# Licensed under a 3-clause BSD style license - see LICENSE.rst

from astropy.wcs import WCS as AstropyWCS

from .core import WCSAxes


# The purpose of the following class is to implement the ``_as_mpl_axes`` method
# which allows us to use the ``projection=`` API for Matplotlib. Once WCSAxes is
# merged into Astropy, we can just add this method directly to the 'real' WCS
# class.

class WCS(AstropyWCS):

    def _as_mpl_axes(self):
        return WCSAxes, {'wcs': self}

    def __iter__(self):
        # Backport of astropy/astropy#3066
        raise TypeError("'{0}' object is not iterable".format(self.__class__.__name__))
