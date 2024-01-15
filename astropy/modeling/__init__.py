# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This subpackage provides a framework for representing models and
performing model evaluation and fitting. It supports 1D and 2D models
and fitting with parameter constraints. It has some predefined models
and fitting routines.
"""

from lazy_loader import attach_stub

__getattr__, __dir__, __all__ = attach_stub(__name__, __file__)
