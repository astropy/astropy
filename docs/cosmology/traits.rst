.. _astropy-cosmology-traits:

****************
Cosmology Traits
****************

.. currentmodule:: astropy.cosmology.traits

The :mod:`~astropy.cosmology.traits` module provides reusable components, called
:term:`traits <trait type>`, that encapsulate specific cosmological properties or
behaviors. For example, the :class:`~astropy.cosmology.traits.HubbleParameter` trait
provides the Hubble constant (``H0``) and related methods, while
:class:`~astropy.cosmology.traits.ScaleFactor`,
:class:`~astropy.cosmology.traits.TemperatureCMB`,
:class:`~astropy.cosmology.traits.DarkEnergyComponent`
:class:`~astropy.cosmology.traits.DarkMatterComponent`
:class:`~astropy.cosmology.traits.BaryonComponent`,
:class:`~astropy.cosmology.traits.PhotonComponent`,
:class:`~astropy.cosmology.traits.TotalComponent`,
:class:`~astropy.cosmology.traits.MatterComponent`, and
:class:`~astropy.cosmology.traits.CriticalDensity` traits provide the scale factor, the temperature or the CMB, the Dark Energy component, and the Dark Matter component,
respectively.

By combining these traits, you can easily construct custom cosmology classes with
precisely the features you need, without having to reimplement common functionality.


Here is an example of how to use the
:class:`~astropy.cosmology.traits.HubbleParameter`,
:class:`~astropy.cosmology.traits.ScaleFactor`,
:class:`~astropy.cosmology.traits.TemperatureCMB`,
:class:`~astropy.cosmology.traits.DarkEnergyComponent`
:class:`~astropy.cosmology.traits.DarkMatterComponent`
:class:`~astropy.cosmology.traits.BaryonComponent`,
:class:`~astropy.cosmology.traits.PhotonComponent`,
:class:`~astropy.cosmology.traits.TotalComponent`,
:class:`~astropy.cosmology.traits.MatterComponent`, and
:class:`~astropy.cosmology.traits.CriticalDensity` traits in custom cosmology classes:

>>> import numpy as np
>>> from astropy import units as u
>>> from astropy.cosmology import Cosmology
>>> from astropy.cosmology.traits import (
...     HubbleParameter,
...     ScaleFactor,
...     TemperatureCMB,
...     DarkEnergyComponent,
...     DarkMatterComponent,
...     BaryonComponent,
...     MatterComponent,
...     CriticalDensity,
...     PhotonComponent,
...     TotalComponent,
... )
>>> class CustomStandardCosmology(Cosmology, HubbleParameter, ScaleFactor, TemperatureCMB,
...                              DarkEnergyComponent, DarkMatterComponent, BaryonComponent,
...                              MatterComponent, CriticalDensity, PhotonComponent, TotalComponent):
...     """Mimics standard ΛCDM cosmology with Planck 2018 parameters."""
...     is_flat = False
...     def __init__(self, H0=67.66, Om0=0.3111, Ode0=0.6889, Ob0=0.0490, Tcmb0=2.7255, Ogamma0=5.38e-5):
...         self.H0 = H0 << (u.km / u.s / u.Mpc)
...         self.Om0 = Om0
...         self.Ode0 = Ode0
...         self.Ob0 = Ob0
...         # Dark matter density at z=0
...         self.Odm0 = Om0 - Ob0
...         # neutrinos and curvature are not provided in this simple example
...         self.Onu0 = 0.0
...         self.Ok0 = 0.0
...         self.Tcmb0 = u.Quantity(Tcmb0, "K")
...         self.Ogamma0 = Ogamma0  # Photon density parameter
...         self.is_flat = bool(abs(self.Otot0 - 1.0) < 1e-12)
...         super().__init__()
...
...     def w(self, z):
...         """Standard cosmological constant."""
...         return -1.0
...
...     def inv_efunc(self, z):
...         """Standard cosmological constant."""
...         zp1 = np.asarray(z) + 1.0
...         return 1.0 / np.sqrt(self.Om0 * zp1**3 + self.Ogamma0 * zp1**4 + self.Ode0)
...
...     @property
...     def Otot0(self):
...         """Total density parameter at z=0 as the sum of the defined components."""
...         return np.float64(1.0 - self.Ok0)
...
...     def Otot(self, z):
...         """Total density parameter at redshift z as the sum of all components."""
...         z = np.asarray(z)
...         return (
...             self.Om(z) + self.Ogamma(z) + self.Ode(z)
...         )

>>> # Create and test standard cosmology
>>> std_cosmo = CustomStandardCosmology()
>>> std_cosmo.H0
<Quantity 67.66 km / (Mpc s)>
>>> std_cosmo.scale_factor(0)
<Quantity 1.>
>>> std_cosmo.Tcmb(1)
<Quantity 5.451 K>
>>> std_cosmo.w(0.5)
-1.0
>>> std_cosmo.Ogamma(0)  # Photon density at z=0
np.float64(5.37...e-05)
>>> std_cosmo.Otot0  # Total density at z=0
np.float64(1.0)
>>> std_cosmo.is_flat
True

.. doctest-requires:: scipy

    >>> # Calculate total density at different redshifts
    >>> std_cosmo.Otot([0, 1, 2])  # Total density at z=0, 1, and 2
    array([1., 1., 1.])
    >>> std_cosmo.Otot(1.5)  # Total density at z=1.5
    np.float64(1...)

Reference/API
*************

.. py:currentmodule:: astropy.cosmology.traits

.. automodapi:: astropy.cosmology.traits
   :inherited-members:
