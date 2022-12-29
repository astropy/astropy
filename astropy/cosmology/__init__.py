# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""``astropy.cosmology`` contains classes and functions for cosmological
distance measures and other cosmology-related calculations.

See the `Astropy documentation
<https://docs.astropy.org/en/latest/cosmology/index.html>`_ for more
detailed usage examples and references.

"""

from . import units

from . import io  # needed before 'realizations'  # isort: split
from . import realizations
from .core import Cosmology, CosmologyError, FlatCosmologyMixin
from .flrw.base import FLRW, FlatFLRWMixin
from .flrw.lambdacdm import FlatLambdaCDM, LambdaCDM
from .flrw.w0cdm import FlatwCDM, wCDM
from .flrw.w0wacdm import Flatw0waCDM, w0waCDM
from .flrw.w0wzcdm import Flatw0wzCDM, w0wzCDM
from .flrw.wpwazpcdm import FlatwpwaCDM, wpwaCDM
from .funcs import cosmology_equal, z_at_value
from .parameter import Parameter
from .realizations import available, default_cosmology

__all__ = [
    "Cosmology",
    "FlatCosmologyMixin",
    "Parameter",
    "CosmologyError",
    # FLRW
    "FLRW",
    "FlatFLRWMixin",
    "LambdaCDM",
    "FlatLambdaCDM",
    "wCDM",
    "FlatwCDM",
    "w0waCDM",
    "Flatw0waCDM",
    "w0wzCDM",
    "Flatw0wzCDM",
    "wpwaCDM",
    "FlatwpwaCDM",
    # realizations
    "available",
    "default_cosmology",
    "WMAP1",  # In `__getattr__`
    "WMAP3",  # In `__getattr__`
    "WMAP5",  # In `__getattr__`
    "WMAP7",  # In `__getattr__`
    "WMAP9",  # In `__getattr__`
    "Planck13",  # In `__getattr__`
    "Planck15",  # In `__getattr__`
    "Planck18",  # In `__getattr__`
    # funcs
    "z_at_value",
    "cosmology_equal",
    # public modules
    "realizations",
    "units",
]


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
    return sorted(__all__ + list(available))
