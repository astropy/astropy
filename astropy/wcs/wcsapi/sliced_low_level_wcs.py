import warnings

from astropy.utils.exceptions import AstropyDeprecationWarning

from .wrappers.sliced_wcs import SlicedLowLevelWCS, sanitize_slices

warnings.warn(
    "SlicedLowLevelWCS has been moved to"
    " astropy.wcs.wcsapi.wrappers.sliced_wcs.SlicedLowLevelWCS, or can be"
    " imported from astropy.wcs.wcsapi.",
    AstropyDeprecationWarning)
