import warnings
warnings.warn(
    "The pywcs namespace is deprecated.  Use astropy.wcs instead.",
    DeprecationWarning)

from astropy.wcs import *
