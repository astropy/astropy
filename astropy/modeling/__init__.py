# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This subpackage provides a framework for representing models and
performing model evaluation and fitting. It supports 1D and 2D models
and fitting with parameter constraints. It has some predefined models
and fitting routines.
"""

from . import fitting
from . import models
from .core import *
from .parameters import *
from .separable import *
