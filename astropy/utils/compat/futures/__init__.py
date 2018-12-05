
import warnings
from concurrent.futures import *

from astropy.utils.exceptions import AstropyDeprecationWarning

warnings.warn("astropy.utils.compat.futures is now deprecated - "
              "use concurrent.futures instead", AstropyDeprecationWarning)
