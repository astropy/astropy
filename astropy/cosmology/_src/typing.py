"""Static typing for :mod:`astropy.cosmology`. PRIVATE API."""
# Licensed under a 3-clause BSD style license - see LICENSE.rst

__all__ = ["_CosmoT"]

from typing import TYPE_CHECKING, TypeVar

if TYPE_CHECKING:
    import astropy.cosmology

_CosmoT = TypeVar("_CosmoT", bound="astropy.cosmology.Cosmology")
"""Type variable for :class:`~astropy.cosmology.Cosmology` and subclasses."""
