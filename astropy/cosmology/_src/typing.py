"""Static typing for :mod:`astropy.cosmology`. PRIVATE API."""
# Licensed under a 3-clause BSD style license - see LICENSE.rst

__all__ = ["CosmoMeta", "FArray", "_CosmoT"]

from collections.abc import Mapping
from typing import TYPE_CHECKING, Any, TypeAlias, TypeVar

import numpy as np
from numpy.typing import NDArray

if TYPE_CHECKING:
    import astropy.cosmology

_CosmoT = TypeVar("_CosmoT", bound="astropy.cosmology.Cosmology")
"""Type variable for :class:`~astropy.cosmology.Cosmology` and subclasses."""

CosmoMeta: TypeAlias = Mapping[Any, Any]
"""Type alias for cosmology metadata."""

FArray: TypeAlias = NDArray[np.floating]
"""Type alias for numpy array of floating dtype."""
