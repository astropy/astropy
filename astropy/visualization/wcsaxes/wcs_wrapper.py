# Licensed under a 3-clause BSD style license - see LICENSE.rst

import warnings

from astropy.wcs import WCS as AstropyWCS
from astropy.utils.exceptions import AstropyDeprecationWarning

from .core import WCSAxes


class WCS(AstropyWCS):

    def __init__(self, *args, **kwargs):
        warnings.warn("The wcsaxes.WCS class has been deprecated - use "
                      "astropy.wcs.WCS instead", AstropyDeprecationWarning)
