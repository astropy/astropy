# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import absolute_import

import warnings
from argparse import *

from ..exceptions import AstropyDeprecationWarning

warnings.warn("astropy.utils.compat.argparse is now deprecated - use the argparse module directly instead", AstropyDeprecationWarning)
