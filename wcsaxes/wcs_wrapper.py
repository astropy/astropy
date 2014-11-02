from astropy.wcs import WCS as AstropyWCS

from .core import WCSAxes


class WCS(AstropyWCS):

    def _as_mpl_axes(self):
        return WCSAxes, {'wcs': self}

    def __iter__(self):
        # Backport of astropy/astropy#3066
        raise TypeError("'{0}' object is not iterable".format(self.__class__.__name__))
