"""Trait for dark matter component of cosmology."""

__all__ = ("DarkMatterComponent",)

from collections.abc import Callable
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray

from astropy.cosmology._src.typing import FArray
from astropy.cosmology._src.utils import aszarr
from astropy.units import Quantity


class DarkMatterComponent:
    """The cosmology has attributes and methods for the dark matter density.

    This trait provides an ``Odm`` method that returns the dark matter
    density parameter (i.e., total matter minus baryons) at redshift ``z``.
    """

    Odm0: float | np.floating
    """Omega dark matter: dark matter density/critical density at z=0."""

    inv_efunc: Callable[[NDArray[Any]], NDArray[Any]]

    def Odm(self, z: Quantity | ArrayLike, /) -> FArray:
        """Return the density parameter for dark matter at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

            .. versionchanged:: 7.0
                Passing z as a keyword argument is deprecated.

            .. versionchanged:: 8.0
               z must be a positional argument.

        Returns
        -------
        Odm : ndarray
            The density of dark matter relative to the critical density at
            each redshift.

        """
        z = aszarr(z)
        return self.Odm0 * (z + 1.0) ** 3 * self.inv_efunc(z) ** 2
