import warnings

from astropy.utils.exceptions import AstropyDeprecationWarning

from .wrappers.sliced_wcs import SlicedLowLevelWCS, sanitize_slices  # noqa: F401

warnings.warn(
    "SlicedLowLevelWCS has been moved to"
    " astropy.wcs.wcsapi.wrappers.sliced_wcs.SlicedLowLevelWCS, or can be"
    " imported from astropy.wcs.wcsapi.",
    AstropyDeprecationWarning,
)
