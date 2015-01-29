# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This subpackage contains utility modules for compatibility with older/newer
versions of python, as well as including some bugfixes for the stdlib that are
important for Astropy.

Note that all public functions in the `astropy.utils.compat.misc` module are
imported here for easier access.
"""

from .misc import *

if not _ASTROPY_SETUP_:
    # Importing this module will also install monkey-patches defined in it
    from .numpycompat import *
