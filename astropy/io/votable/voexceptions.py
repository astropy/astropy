from ...utils.exceptions import AstropyDeprecationWarning

# Licensed under a 3-clause BSD style license - see LICENSE.rst
if not _ASTROPY_SETUP_:
    import warnings
    warnings.warn(
        "astropy.io.votable.voexceptions is deprecated. "
        "Use astropy.io.votable.exceptions",
        AstropyDeprecationWarning)

from .exceptions import *
