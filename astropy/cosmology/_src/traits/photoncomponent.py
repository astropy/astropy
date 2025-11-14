# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Photon component."""

__all__ = ["PhotonComponent"]

from collections.abc import Callable
from typing import Any

from numpy.typing import ArrayLike, NDArray

from astropy.cosmology._src.typing import FArray
from astropy.cosmology._src.utils import aszarr, deprecated_keywords
from astropy.units import Quantity


class PhotonComponent:
    """The cosmology has attributes and methods for the photon density."""

    Ogamma0: float
    """Omega gamma; the density/critical density of photons at z=0."""

    inv_efunc: Callable[[NDArray[Any]], NDArray[Any]]

    @deprecated_keywords("z", since="7.0")
    def Ogamma(self, z: Quantity | ArrayLike) -> FArray:
        """Return the density parameter for photons at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

            .. versionchanged:: 7.0
                Passing z as a keyword argument is deprecated.

        Returns
        -------
        Ogamma : array
            The energy density of photons relative to the critical density at
            each redshift.
        """
        z = aszarr(z)
        return self.Ogamma0 * (z + 1.0) ** 4 * self.inv_efunc(z) ** 2
