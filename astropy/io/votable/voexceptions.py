# Licensed under a 3-clause BSD style license - see LICENSE.rst
import warnings
warnings.warn(
    "astropy.io.votable.voexceptions is deprecated. "
    "Use astropy.io.votable.exceptions",
    DeprecationWarning)

from .exceptions import *
