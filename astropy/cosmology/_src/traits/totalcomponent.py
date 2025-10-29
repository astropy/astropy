# Licensed under a 3-clause BSD style license - see LICENSE.rst

__all__ = ["TotalComponent"]

from abc import abstractmethod

from numpy.typing import ArrayLike

from astropy.cosmology._src.typing import FArray
from astropy.cosmology._src.utils import aszarr, deprecated_keywords
from astropy.units import Quantity


class TotalComponent:
    """The cosmology has attributes and methods for the total density.

    This trait has the abstract ``Otot`` method that returns the total density
    parameter at redshift ``z``. It should be the sum of all other components.
    """

    @property
    @abstractmethod
    def Otot0(self) -> float:
        """Omega total; the total density/critical density at z=0."""
        raise NotImplementedError  # pragma: no cover

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
        # Sum all density components at redshift z
        # This relies on the cosmology having the appropriate component traits
        return self.Om(z) + self.Ogamma(z) + self.Onu(z) + self.Ode(z) + self.Ok(z)
