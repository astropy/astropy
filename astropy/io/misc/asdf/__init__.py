import warnings

from astropy.utils.exceptions import AstropyDeprecationWarning

warnings.warn(
    "ASDF functionality for astropy is being moved out of the astropy"
    " package to the new asdf-astropy package. Please use this package"
    " instead. astropy.io.misc.asdf is deprecated since astropy 5.1 and will be removed in a future release.",
    AstropyDeprecationWarning
)
