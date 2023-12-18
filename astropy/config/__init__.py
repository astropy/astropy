# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module contains configuration and setup utilities for the
Astropy project.
"""

from lazy_loader import attach_stub

__getattr__, __dir__, __all__ = attach_stub(__name__, __file__)
