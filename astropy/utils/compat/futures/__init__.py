from concurrent.futures import *

import warnings
from ...exceptions import AstropyDeprecationWarning

warnings.warn("astropy.utils.compat.futures is now deprecated - "
              "use concurrent.futures instead", AstropyDeprecationWarning)
