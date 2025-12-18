"""Matter component."""

import numpy as np
from numpy.typing import ArrayLike

from astropy.cosmology._src.typing import FArray
from astropy.cosmology._src.utils import aszarr
from astropy.units import Quantity

__all__ = ("MatterComponent",)


class MatterComponent:
    """The cosmology has attributes and methods for the matter density."""

    Om0: float | np.floating
    """Omega matter; matter density/critical density at z=0."""

    def Om(self, z: Quantity | ArrayLike, /) -> FArray:
        """Return the density parameter for non-relativistic matter at redshift ``z``.

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
        Om : ndarray
            The density of non-relativistic matter relative to the critical
            density at each redshift.

        Notes
        -----
        This does not include neutrinos, even if non-relativistic at the
        redshift of interest.
        """
        z = aszarr(z)
        return self.Om0 * (z + 1.0) ** 3 * self.inv_efunc(z) ** 2
