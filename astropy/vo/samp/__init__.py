import warnings

from ...utils.exceptions import AstropyDeprecationWarning
from ...samp import *  # noqa

warnings.warn('The astropy.samp module has now been moved to astropy.samp',
              AstropyDeprecationWarning)
