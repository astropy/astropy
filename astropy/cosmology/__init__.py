# Licensed under a 3-clause BSD style license - see LICENSE.rst
""" astropy.cosmology contains classes and functions for cosmological
distance measures and other cosmology-related calculations.

See the `Astropy documentation
<https://docs.astropy.org/en/latest/cosmology/index.html>`_ for more
detailed usage examples and references.

"""

from . import core, flrw, funcs, units, utils

from . import io  # needed before 'realizations'  # isort: split
from . import realizations
from .core import *
from .flrw import *
from .funcs import *
from .realizations import *
from .utils import *

__all__ = (core.__all__ + flrw.__all__       # cosmology classes
           + realizations.__all__            # instances thereof
           + funcs.__all__ + utils.__all__)  # utils
