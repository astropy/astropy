# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Critical density component."""

__all__ = ("CriticalDensity",)

from collections.abc import Callable
from typing import Any

from numpy.typing import ArrayLike, NDArray

from astropy.units import Quantity


class CriticalDensity:
    """The object has attributes and methods for the critical density.

    This is a trait class; it is not meant to be instantiated directly, but instead to
    be used as a mixin to other classes.

    Examples
    --------
    For an example of a real cosmology that implements this trait, see
    :class:`~astropy.cosmology.LambdaCDM`. Here we will define an illustrative example
    class that meets the minimum API requirements, but is not cosmologically meaningful:

    >>> import numpy as np
    >>> import astropy.units as u
    >>> from astropy.cosmology.traits import MatterComponent
    >>> import dataclasses

    >>> @dataclasses.dataclass(frozen=True)
    ... class ExampleHasCriticalDensity(CriticalDensity):
    ...     critical_density0: u.Quantity
    ...     def efunc(self, z): return np.ones_like(z)

    >>> cosmo = ExampleHasCriticalDensity(critical_density0=1e-29 * u.g / u.cm**3)
    >>> cosmo.critical_density0
    <Quantity 1.e-29 g / cm3>

    >>> cosmo.critical_density([0.0, 1.0, 2.0])
    <Quantity [1.e-29, 1.e-29, 1.e-29] g / cm3>

    """

    critical_density0: Quantity
    """Critical density at redshift 0."""

    efunc: Callable[[Any], NDArray[Any]]

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
