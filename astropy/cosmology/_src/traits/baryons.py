# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Baryon component.

This is private API. See `~astropy.cosmology.traits` for public API.

"""

__all__ = ["_BaryonComponent"]

from collections.abc import Callable
from typing import Any

from numpy.typing import ArrayLike, NDArray

from astropy.cosmology._src.utils import aszarr, deprecated_keywords
from astropy.units import Quantity


class _BaryonComponent:
    """The cosmology has attributes and methods for the baryon density."""

    Ob0: float
    """Omega baryons: density of baryonic matter in units of the critical density at z=0."""

    inv_efunc: Callable[[NDArray[Any]], NDArray[Any]]

    @deprecated_keywords("z", since="7.0")
    def Ob(self, z: Quantity | ArrayLike) -> NDArray[Any] | float:
        """Return the density parameter for baryonic matter at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

            .. versionchanged:: 7.0
                Passing z as a keyword argument is deprecated.

        Returns
        -------
        Ob : ndarray or float
            The density of baryonic matter relative to the critical density at
            each redshift.
            Returns `float` if the input is scalar.

        """
        z = aszarr(z)
        return self.Ob0 * (z + 1.0) ** 3 * self.inv_efunc(z) ** 2
