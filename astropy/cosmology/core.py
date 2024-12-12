# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Astropy cosmology core module."""

__all__ = ["Cosmology", "CosmologyError", "FlatCosmologyMixin"]

from ._src.core import (
    _COSMOLOGY_CLASSES,  # noqa: F401
    Cosmology,
    CosmologyError,
    FlatCosmologyMixin,
    dataclass_decorator,  # noqa: F401
)
