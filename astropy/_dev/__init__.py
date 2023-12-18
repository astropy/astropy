# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This package contains utilities that are only used when developing astropy
in a copy of the source repository.

These files are not installed, and should not be assumed to exist at runtime.
"""

from lazy_loader import attach_stub

__getattr__, __dir__, __all__ = attach_stub(__name__, __file__)
