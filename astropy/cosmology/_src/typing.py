"""Static typing for :mod:`astropy.cosmology`. PRIVATE API."""
# Licensed under a 3-clause BSD style license - see LICENSE.rst

__all__ = ["FArray", "_CosmoT"]

from typing import TYPE_CHECKING, TypeAlias, TypeVar

import numpy as np
from numpy.typing import NDArray

if TYPE_CHECKING:
    import astropy.cosmology

_CosmoT = TypeVar("_CosmoT", bound="astropy.cosmology.Cosmology")
"""Type variable for :class:`~astropy.cosmology.Cosmology` and subclasses."""

FArray: TypeAlias = NDArray[np.floating]
