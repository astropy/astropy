# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This sub-package contains classes and functions for creating distributions that
work similar to `~astropy.units.Quantity` or array objects, but can propagate
uncertainties.
"""

from lazy_loader import attach_stub

__getattr__, __dir__, __all__ = attach_stub(__name__, __file__)
