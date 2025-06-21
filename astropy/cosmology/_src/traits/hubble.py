# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Hubble parameter trait.

This is private API. See `~astropy.cosmology.traits` for public API.
"""

__all__ = ["HubbleParameter"]

from collections.abc import Callable
from functools import cached_property
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray

import astropy.units as u
from astropy import constants as const
from astropy.cosmology._src.utils import deprecated_keywords
from astropy.units import Quantity


class HubbleParameter:
    """The object has attributes and methods for the Hubble parameter."""

    H0: Quantity
    """Hubble Parameter at redshift 0."""

    efunc: Callable[[Any], NDArray[Any]]

    @deprecated_keywords("z", since="7.0")
    def H(self, z: Quantity | ArrayLike) -> Quantity:
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
    def h(self) -> np.floating[Any]:
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
