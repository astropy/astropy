# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Hubble parameter trait.

This is private API. See `~astropy.cosmology.traits` for public API.
"""

__all__ = ["_HubbleParameter"]

from collections.abc import Callable
from typing import Any

from numpy.typing import ArrayLike, NDArray

from astropy.cosmology._src.utils import deprecated_keywords
from astropy.units import Quantity


class _HubbleParameter:
    """The object has attributes and methods for the Hubble parameter."""

    H0: Quantity
    """Hubble constant at redshift 0."""

    efunc: Callable[[Any], NDArray[Any]]

    @deprecated_keywords("z", since="7.0")
    def hubble_parameter(self, z: Quantity | ArrayLike) -> Quantity:
        """Hubble parameter at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

        Returns
        -------
        H : Quantity ['frequency']
            Hubble parameter at each input redshift.
        """
        return self.H0 * self.efunc(z)
