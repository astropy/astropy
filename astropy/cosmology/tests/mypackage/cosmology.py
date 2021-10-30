# -*- coding: utf-8 -*-

"""Cosmology classes and instances for package ``mypackage``."""

# STDLIB
from collections import UserDict

# THIRD PARTY
import numpy as np


class MyCosmology(UserDict):
    """Cosmology, from ``mypackage``.

    Parameters
    ----------
    **kw
        values for `MyCosmology`.
    """

    def __init__(self, **kw):
        super().__init__()
        self.update(kw)

    def __getattr__(self, key):
        """Get attributes like items."""
        if key in self:
            return self[key]
        super().__getattribute__(key)

    def age_in_Gyr(z):
        """Cosmology age in Gyr. Returns `NotImplemented`."""
        return NotImplemented

    def __eq__(self, other):
        """Equality check."""
        if not isinstance(other, MyCosmology):
            return NotImplemented
        return all(np.all(v == other[k]) for k, v in self.items())


myplanck = MyCosmology(
    name="planck",
    hubble_parameter=67.66,  # km/s/Mpc
    initial_dm_density=0.2607,  # X/rho_critical
    initial_baryon_density=0.04897,  # X/rho_critical
    initial_matter_density=0.30966,  # X/rho_critical
    n=0.9665,
    sigma8=0.8102,
    z_reionization=7.82,  # redshift
    current_age=13.787,  # Gyr
    initial_temperature=2.7255,  # Kelvin
    Neff=3.046,
    neutrino_masses=[0., 0., 0.06],  # eV
)
