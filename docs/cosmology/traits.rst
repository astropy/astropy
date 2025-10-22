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
:class:`~astropy.cosmology.traits.DarkEnergyComponent` and
:class:`~astropy.cosmology.traits.DarkMatterComponent` provide the scale factor, the
temperature or the CMB, the Dark Energy component, and the Dark Matter component,
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
:class:`~astropy.cosmology.traits.MatterComponent`, and
:class:`~astropy.cosmology.traits.CriticalDensity` traits in custom cosmology classes:

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
... )
>>> class CustomStandardCosmology(Cosmology, HubbleParameter, ScaleFactor, TemperatureCMB,
...                              DarkEnergyComponent, DarkMatterComponent, BaryonComponent,
...                              MatterComponent, CriticalDensity):
...     """Mimics standard ΛCDM cosmology with Planck 2018 parameters."""
...     def __init__(self, H0=67.66, Om0=0.3111, Ode0=0.6889, Ob0=0.0490, Tcmb0=2.7255):
...         self.H0 = H0 << (u.km / u.s / u.Mpc)
...         self.Om0 = Om0
...         self.Ode0 = Ode0
...         self.Ob0 = Ob0
...         self.Odm = Om0 - Ob0  # Dark matter density
...         self.Tcmb0 = u.Quantity(Tcmb0, "K")
...         super().__init__()
...
...     def w(self, z):
...         """Standard cosmological constant."""
...         return -1.0
...
...     is_flat = True  # Standard ΛCDM is flat

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

Reference/API
*************

.. py:currentmodule:: astropy.cosmology.traits

.. automodapi:: astropy.cosmology.traits
   :inherited-members:
