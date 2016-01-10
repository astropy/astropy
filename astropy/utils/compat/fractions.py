# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import absolute_import

import warnings
from fractions import *

from ..exceptions import AstropyDeprecationWarning

warnings.warn("astropy.utils.compat.fractions is now deprecated - use the fractions module directly instead", AstropyDeprecationWarning)
