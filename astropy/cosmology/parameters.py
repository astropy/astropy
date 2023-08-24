# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""This module contains dictionaries with sets of parameters for a given cosmology.

The list of cosmologies available are given by the tuple `available`.
"""

# STDLIB
import sys
from types import MappingProxyType

# LOCAL
from .realizations import available

__all__ = ["available"]
__all__ += [  # noqa: F822
    "WMAP1",
    "WMAP3",
    "WMAP5",
    "WMAP7",
    "WMAP9",
    "Planck13",
    "Planck15",
    "Planck18",
]


def __getattr__(name):
    """Get parameters of cosmology representations with lazy import from ``PEP 562``."""
    from astropy.cosmology import realizations

    cosmo = getattr(realizations, name)
    m = cosmo.to_format("mapping", cosmology_as_str=True, move_from_meta=True)
    proxy = MappingProxyType(m)

    # Cache in this module so `__getattr__` is only called once per `name`.
    setattr(sys.modules[__name__], name, proxy)

    return proxy


def __dir__():
    """Directory, including lazily-imported objects."""
    return __all__
