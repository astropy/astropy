from __future__ import absolute_import, print_function

import warnings
from subprocess import *

from ..exceptions import AstropyDeprecationWarning

warnings.warn("astropy.utils.compat.subprocess is now deprecated - use the Python subprocess module directly instead", AstropyDeprecationWarning)
