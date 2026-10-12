import warnings
from astropy.utils.exceptions import AstropyDeprecationWarning

warnings.warn(
    "astropy.extern.configobj is deprecated, import from the "
    "configobj or validate modules directly instead",
    AstropyDeprecationWarning
)
