# Licensed under a 3-clause BSD style license - see LICENSE.rst
""" astropy.cosmology contains classes and functions for cosmological
distance measures and other cosmology-related calculations.

See the `Astropy documentation
<https://docs.astropy.org/en/latest/cosmology/index.html>`_ for more
detailed usage examples and references.

"""

from . import core, funcs, realizations
from .core import *
from .funcs import *
from .realizations import *

__all__ = core.__all__ + realizations.__all__ + funcs.__all__
