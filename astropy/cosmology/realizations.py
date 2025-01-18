# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Built-in cosmologies.

See :attr:`~astropy.cosmology.realizations.available` for a full list.
"""

__all__ = [  # noqa: F822, RUF022
    "available",
    "default_cosmology",
    # ----
    "WMAP1",
    "WMAP3",
    "WMAP5",
    "WMAP7",
    "WMAP9",
    "Planck13",
    "Planck15",
    "Planck18",
]

import pathlib
import sys

from astropy.utils.data import get_pkg_data_path

from ._src.core import Cosmology
from ._src.default import default_cosmology

__doctest_requires__ = {"*": ["scipy"]}

_COSMOLOGY_DATA_DIR = pathlib.Path(
    get_pkg_data_path("cosmology", "data", package="astropy")
)
available = (
    "WMAP1",
    "WMAP3",
    "WMAP5",
    "WMAP7",
    "WMAP9",
    "Planck13",
    "Planck15",
    "Planck18",
)


def __getattr__(name: str) -> Cosmology:
    """Make specific realizations from data files with lazy import from ``PEP 562``.

    Raises
    ------
    AttributeError
        If "name" is not in :mod:`astropy.cosmology.realizations`
    """
    if name not in available:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}.")

    cosmo = Cosmology.read(
        str(_COSMOLOGY_DATA_DIR / name) + ".ecsv", format="ascii.ecsv"
    )
    object.__setattr__(
        cosmo,
        "__doc__",
        f"{name} instance of {cosmo.__class__.__qualname__} "
        f"cosmology\n(from {cosmo.meta['reference']})",
    )

    # Cache in this module so `__getattr__` is only called once per `name`.
    setattr(sys.modules[__name__], name, cosmo)

    return cosmo


def __dir__() -> list[str]:
    """Directory, including lazily-imported objects."""
    return __all__
