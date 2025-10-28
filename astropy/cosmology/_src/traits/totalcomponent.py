# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Total density component.

This is private API. See `~astropy.cosmology.traits` for public API.

"""

__all__ = ["_TotalComponent"]

from collections.abc import Callable

from numpy.typing import ArrayLike

from astropy.cosmology._src.typing import FArray
from astropy.cosmology._src.utils import aszarr, deprecated_keywords
from astropy.units import Quantity


class _TotalComponent:
    """The cosmology has attributes and methods for the total density.

    This trait provides an ``Otot`` method that returns the total density
    parameter (i.e., sum of all components: matter, photons, neutrinos, dark
    energy, and curvature) at redshift ``z``.
    """

    Om0: float
    """Omega matter; matter density/critical density at z=0."""

    Ogamma0: float
    """Omega gamma; photon density/critical density at z=0."""

    Onu0: float
    """Omega nu; neutrino density/critical density at z=0."""

    Ode0: float
    """Omega dark energy; dark energy density/critical density at z=0."""

    Ok0: float
    """Omega curvature; curvature density/critical density at z=0."""

    Om: Callable[[Quantity | ArrayLike], FArray]
    Ogamma: Callable[[Quantity | ArrayLike], FArray]
    Onu: Callable[[Quantity | ArrayLike], FArray]
    Ode: Callable[[Quantity | ArrayLike], FArray]
    Ok: Callable[[Quantity | ArrayLike], FArray]

    @property
    def Otot0(self) -> float:
        """Omega total; the total density/critical density at z=0."""
        return self.Om0 + self.Ogamma0 + self.Onu0 + self.Ode0 + self.Ok0

    @deprecated_keywords("z", since="7.0")
    def Otot(self, z: Quantity | ArrayLike) -> FArray:
        """The total density parameter at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshifts.

            .. versionchanged:: 7.0
                Passing z as a keyword argument is deprecated.

        Returns
        -------
        Otot : array
            The total density relative to the critical density at each
            redshift.
        """
        z = aszarr(z)
        return self.Om(z) + self.Ogamma(z) + self.Onu(z) + self.Ode(z) + self.Ok(z)
