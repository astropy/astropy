# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This subpackage contains utility modules for compatibility with older/newer
versions of python, as well as including some bugfixes for the stdlib that are
important for Astropy.

Note that all public functions in the `astropy.utils.compat.misc` module are
imported here for easier access.

The content of this module is solely for internal use of ``astropy``
and subject to changes without deprecations. Do not use it in external
packages or code.

"""

from lazy_loader import attach_stub

__getattr__, __dir__, __all__ = attach_stub(__name__, __file__)
