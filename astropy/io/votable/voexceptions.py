import warnings
warnings.warn(
    "astropy.io.votable.voexceptions is deprecated. "
    "Use astropy.io.votable.exceptions",
    DeprecationWarning)

from .exceptions import *
