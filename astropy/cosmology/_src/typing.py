"""Static typing for :mod:`astropy.cosmology`. PRIVATE API."""
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import annotations

__all__ = ["_CosmoT"]

from typing import TYPE_CHECKING, TypeVar

if TYPE_CHECKING:
    from astropy.cosmology import Cosmology

_CosmoT = TypeVar("_CosmoT", bound="Cosmology")
"""Type variable for :class:`~astropy.cosmology.Cosmology` and subclasses."""
