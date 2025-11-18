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
:class:`~astropy.cosmology.traits.NeutrinoComponent`,
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
:class:`~astropy.cosmology.traits.NeutrinoComponent`,
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
...     CurvatureComponent,
...     PhotonComponent,
...     NeutrinoComponent,
...     TotalComponent,
... )
>>> NEUTRINO_FERMI_DIRAC_CORRECTION = 0.22710731766  # 7/8 (4/11)^4/3
>>> class CustomStandardCosmology(Cosmology, HubbleParameter, ScaleFactor, TemperatureCMB,
...                              DarkEnergyComponent, DarkMatterComponent, BaryonComponent,
...                              MatterComponent, CriticalDensity, CurvatureComponent,
...                              PhotonComponent, TotalComponent, NeutrinoComponent):
...     """Mimics standard ΛCDM cosmology with Planck 2018 parameters."""
...     def __init__(self, H0=67.66, Om0=0.3111, Ob0=0.0490, Tcmb0=2.7255, Neff=3.046):
...         import astropy.constants as const
...         from astropy.cosmology._src.flrw.base import a_B_c2
...         from math import pi
...         # Set basic parameters
...         self.H0 = H0 << (u.km / u.s / u.Mpc)
...         self.Om0 = Om0
...         self.Ob0 = Ob0
...         self.Odm0 = Om0 - Ob0
...         self.Neff = Neff
...         self.m_nu = None  # massless neutrinos
...         self.Tcmb0 = u.Quantity(Tcmb0, "K")
...         # Compute critical_density0 (needed by CriticalDensity trait)
...         self.critical_density0 = (3.0 * self.H0**2 / (8.0 * pi * const.G)).cgs
...         # Compute Ogamma0 (needed by PhotonComponent trait)
...         self.Ogamma0 = a_B_c2 * self.Tcmb0.value**4 / self.critical_density0.value
...         self.Ode0 = 0.0  # temporary
...         super().__init__()
...         # Compute Ode0 to ensure flatness (Ok0 will be ~0)
...         self.Ode0 = 1.0 - (self.Om0 + self.Ogamma0 + self.Onu0)
...
...     @property
...     def Ok0(self):
...         """Omega curvature (0 for flat cosmology)."""
...         return 0.0
...
...     def w(self, z):
...         """Standard cosmological constant."""
...         return -1.0
...
...     def inv_efunc(self, z):
...         """Inverse E-function for standard ΛCDM."""
...         # Total radiation density (photons + neutrinos)
...         Or = self.Ogamma0 + self.Onu0
...         zp1 = np.asarray(z) + 1.0
...         return 1.0 / np.sqrt(Or * zp1**4 + self.Om0 * zp1**3 + self.Ode0)
...
...     @property
...     def Otot0(self):
...         """Total density parameter at z=0."""
...         return self.Om0 + self.Ogamma0 + self.Onu0 + self.Ode0 + self.Ok0
...
...     def Otot(self, z):
...         """Total density parameter at redshift z."""
...         return np.ones_like(np.asarray(z), dtype=float)
...
...     @property
...     def has_massive_nu(self):
...         return False
...     @property
...     def Onu0(self):
...         return NEUTRINO_FERMI_DIRAC_CORRECTION * self.Neff * self.Ogamma0
...     def nu_relative_density(self, z):
...         return NEUTRINO_FERMI_DIRAC_CORRECTION * self.Neff * np.ones_like(np.asarray(z))
...     def Ogamma(self, z):
...         return self.Ogamma0 * (np.asarray(z) + 1.0) ** 4
...
...     @property
...     def is_flat(self):
...         """Return True if the cosmology is flat."""
...         return bool(abs(self.Ok0) < 1e-12)

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
np.float64(5.40...e-05)
>>> std_cosmo.Otot0  # Total density at z=0
np.float64(1.0)
>>> std_cosmo.is_flat
True
>>> std_cosmo.Neff  # Effective number of neutrino species
3.046
>>> std_cosmo.Onu(0)  # Neutrino density at z=0
np.float64(3.7...e-05)

.. doctest-requires:: scipy

    >>> # Calculate total density at different redshifts
    >>> std_cosmo.Otot([0, 1, 2])  # Total density at z=0, 1, and 2
    array([1., 1., 1.])
    >>> std_cosmo.Otot(1.5)  # Total density at z=1.5
    array(1...)

Reference/API
*************

.. py:currentmodule:: astropy.cosmology.traits

.. automodapi:: astropy.cosmology.traits
   :inherited-members:
