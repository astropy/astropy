from astropy.wcs import WCS as AstropyWCS

from .core import WCSAxes


class WCS(AstropyWCS):
    
    def _as_mpl_axes(self):
        return WCSAxes, {'wcs': self}