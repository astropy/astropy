# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Static typing for :mod:`astropy.cosmology`."""

from __future__ import annotations

from typing import TYPE_CHECKING, TypeVar

if TYPE_CHECKING:
    from astropy.cosmology import Cosmology

_CosmoT = TypeVar("_CosmoT", bound="Cosmology")  # noqa: PYI018
"""Type variable for :class:`~astropy.cosmology.Cosmology` and subclasses."""
