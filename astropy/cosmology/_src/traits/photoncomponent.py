# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Photon component."""

__all__ = ("PhotonComponent",)

from collections.abc import Callable
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray

from astropy.cosmology._src.typing import FArray
from astropy.cosmology._src.utils import aszarr
from astropy.units import Quantity


class PhotonComponent:
    """The cosmology has attributes and methods for the photon density.

    This is a trait class; it is not meant to be instantiated directly, but instead to
    be used as a mixin to other classes.

    Examples
    --------
    For an example of a real cosmology that implements this trait, see
    :class:`~astropy.cosmology.LambdaCDM`. Here we will define an illustrative example
    class that meets the minimum API requirements, but is not cosmologically meaningful:

    >>> from astropy.cosmology.traits import MatterComponent
    >>> import dataclasses

    >>> @dataclasses.dataclass(frozen=True)
    ... class ExampleHasPhotons(PhotonComponent):
    ...     Ogamma0: float
    ...     def inv_efunc(self, z): return 1.0  # necessary for Ogamma(z)

    >>> cosmo = ExampleHasPhotons(Ogamma0=5e-5)
    >>> cosmo.Ogamma0
    5e-05

    >>> cosmo.Ogamma([0.0, 1.0, 2.0])
    array([5.00e-05, 8.00e-04, 4.05e-03])

    """

    Ogamma0: float | np.floating
    """Omega gamma; the density/critical density of photons at z=0."""

    inv_efunc: Callable[[NDArray[Any]], NDArray[Any]]

    def Ogamma(self, z: Quantity | ArrayLike, /) -> FArray:
        """Return the density parameter for photons at redshift ``z``.

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
        Ogamma : array
            The energy density of photons relative to the critical density at
            each redshift.
        """
        z = aszarr(z)
        return self.Ogamma0 * (z + 1.0) ** 4 * self.inv_efunc(z) ** 2
