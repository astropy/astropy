# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Baryon component."""

__all__ = ("BaryonComponent",)

from collections.abc import Callable
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray

from astropy.cosmology._src.typing import FArray
from astropy.cosmology._src.utils import aszarr
from astropy.units import Quantity


class BaryonComponent:
    """The cosmology has attributes and methods for the baryon density.

    Examples
    --------
    For an example of a real cosmology that implements this trait, see
    :class:`~astropy.cosmology.LambdaCDM`. Here we will define an illustrative example
    class that meets the minimum API requirements, but is not cosmologically meaningful:

    >>> from astropy.cosmology.traits import BaryonComponent
    >>> import dataclasses

    >>> @dataclasses.dataclass(frozen=True)
    ... class ExampleHasBaryons(BaryonComponent):
    ...     Ob0: float
    ...     def inv_efunc(self, z): return 1.0  # necessary for Ob(z)

    >>> cosmo = ExampleHasBaryons(Ob0=0.05)
    >>> cosmo.Ob0
    0.05

    >>> cosmo.Ob([0.0, 1.0, 2.0])
    array([0.05, 0.4 , 1.35])

    """

    Ob0: float | np.floating
    """Omega baryons: density of baryonic matter in units of the critical density at z=0."""

    inv_efunc: Callable[[NDArray[Any]], NDArray[Any]]

    def Ob(self, z: Quantity | ArrayLike, /) -> FArray:
        """Return the density parameter for baryonic matter at redshift ``z``.

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
        Ob : ndarray
            The density of baryonic matter relative to the critical density at
            each redshift.

        """
        z = aszarr(z)
        return self.Ob0 * (z + 1.0) ** 3 * self.inv_efunc(z) ** 2
