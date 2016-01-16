# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import absolute_import

import warnings
from gzip import *

from ..exceptions import AstropyDeprecationWarning

warnings.warn("astropy.utils.compat.gzip is now deprecated - use the gzip module directly instead", AstropyDeprecationWarning)
