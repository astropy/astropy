# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Hubble parameter trait.

This is private API. See `~astropy.cosmology.traits` for public API.
"""

__all__ = ("HubbleParameter",)

from collections.abc import Callable
from functools import cached_property
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray

import astropy.units as u
from astropy import constants as const
from astropy.cosmology._src.typing import FArray
from astropy.units import Quantity


class HubbleParameter:
    """The object has attributes and methods for the Hubble parameter.

    This is a trait class; it is not meant to be instantiated directly, but instead to
    be used as a mixin to other classes.

    Examples
    --------
    For an example of a real cosmology that implements this trait, see
    :class:`~astropy.cosmology.LambdaCDM`. Here we will define an illustrative example
    class that meets the minimum API requirements, but is not cosmologically meaningful:

    >>> import astropy.units as u
    >>> from astropy.cosmology.traits import HubbleParameter
    >>> import dataclasses

    >>> @dataclasses.dataclass(frozen=True)
    ... class ExampleHasHubble(HubbleParameter):
    ...     H0: u.Quantity
    ...     def efunc(self, z): return np.ones_like(z)
    ...     def inv_efunc(self, z): return np.ones_like(z)

    >>> cosmo = ExampleHasHubble(H0=70 * u.km / u.s / u.Mpc)
    >>> cosmo.H0
    <Quantity 70. km / (Mpc s)>

    >>> cosmo.H([0.0, 1.0, 2.0])
    <Quantity [70., 70., 70.] km / (Mpc s)>

    >>> float(cosmo.h)
    0.7

    >>> cosmo.hubble_time.round(2)
    <Quantity 13.97 Gyr>

    >>> cosmo.hubble_distance.round(2)
    <Quantity 4282.75 Mpc>

    """

    H0: Quantity
    """Hubble Parameter at redshift 0."""

    efunc: Callable[[Any], NDArray[Any]]

    inv_efunc: Callable[[Any], FArray | float]

    def H(self, z: Quantity | ArrayLike, /) -> Quantity:
        """Hubble parameter at redshift ``z``.

        Parameters
        ----------
        z : Quantity-like ['redshift'], array-like
            Input redshift.

        Returns
        -------
        H : Quantity ['frequency']
            Hubble parameter at each input redshift.
        """
        return self.H0 * self.efunc(z)

    @cached_property
    def h(self) -> np.floating:
        """Dimensionless Hubble constant: h = H_0 / 100 [km/sec/Mpc]."""
        return self.H0.to_value("km/(s Mpc)") / 100.0

    @cached_property
    def hubble_time(self) -> u.Quantity:
        """Hubble time."""
        return (1 / self.H0).to(u.Gyr)

    @cached_property
    def hubble_distance(self) -> u.Quantity:
        """Hubble distance."""
        return (const.c / self.H0).to(u.Mpc)
