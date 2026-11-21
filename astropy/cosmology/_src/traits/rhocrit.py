# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Critical density component."""

__all__ = ("CriticalDensity",)


from numpy.typing import ArrayLike

from astropy.units import Quantity

from .hubble import _HasHoverH0


class CriticalDensity(_HasHoverH0):
    """The object has attributes and methods for the critical density."""

    critical_density0: Quantity
    """Critical density at redshift 0."""

    def critical_density(self, z: Quantity | ArrayLike, /) -> Quantity:
        """Critical density in grams per cubic cm at redshift ``z``.

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
        rho : Quantity ['mass density']
            Critical density at each input redshift.
        """
        return self.critical_density0 * self.efunc(z) ** 2
