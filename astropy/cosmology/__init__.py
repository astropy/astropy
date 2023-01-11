# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""``astropy.cosmology`` contains classes and functions for cosmological
distance measures and other cosmology-related calculations.

See the `Astropy documentation
<https://docs.astropy.org/en/latest/cosmology/index.html>`_ for more
detailed usage examples and references.

"""

from . import core, flrw, funcs, parameter, units, utils

from . import io  # needed before 'realizations'  # isort: split
from . import realizations
from .core import *
from .flrw import *
from .funcs import *
from .parameter import *
from .realizations import available, default_cosmology
from .utils import *

__all__ = (
    core.__all__
    + flrw.__all__  # cosmology classes
    + realizations.__all__  # instances thereof
    + ["units"]
    # utils
    + funcs.__all__
    + parameter.__all__
    + utils.__all__
)


def __getattr__(name):
    """Get realizations using lazy import from
    `PEP 562 <https://www.python.org/dev/peps/pep-0562/>`_.

    Raises
    ------
    AttributeError
        If "name" is not in :mod:`astropy.cosmology.realizations`
    """
    if name not in available:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}.")

    return getattr(realizations, name)


def __dir__():
    """Directory, including lazily-imported objects."""
    return __all__
