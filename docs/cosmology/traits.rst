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
:class:`~astropy.cosmology.traits.TemperatureCMB` and
:class:`~astropy.cosmology.traits.DarkEnergyComponent` provide the scale factor, the
temperature or the CMB, and the Dark Energy component, respectively.

By combining these traits, you can easily construct custom cosmology classes with
precisely the features you need, without having to reimplement common functionality.


Here is an example of how to use the
:class:`~astropy.cosmology.traits.HubbleParameter`,
:class:`~astropy.cosmology.traits.ScaleFactor`,
:class:`~astropy.cosmology.traits.TemperatureCMB` and
:class:`~astropy.cosmology.traits.DarkEnergyComponent` traits in a custom cosmology class:

>>> import astropy.units as u
>>> from astropy.cosmology.traits import HubbleParameter, ScaleFactor, TemperatureCMB, DarkEnergyComponent
>>> from astropy.cosmology import Cosmology
>>> import numpy as np
>>>
>>> class CustomCosmology(Cosmology, HubbleParameter, ScaleFactor, TemperatureCMB, DarkEnergyComponent):
...     def __init__(self, Om0, Ode0, H0=70, Tcmb0=2.725):
...         self.H0 = H0 << (u.km / u.s / u.Mpc)
...         self.Om0 = Om0
...         self.Ode0 = Ode0
...         self.Tcmb0 = u.Quantity(Tcmb0, "K")
...         super().__init__()
...
...     def w(self, z):
...         # Example: equation of state varies with redshift
...         return -1 + 0.1 * np.sin(z)
...
...     is_flat = False

>>> cosmo = CustomCosmology(H0=70, Om0=0.3, Ode0=0.7)
>>> cosmo.H0
<Quantity 70. km / (Mpc s)>
>>> cosmo.scale_factor(0)
<Quantity 1.>
>>> cosmo.Tcmb(1)
<Quantity 5.45 K>
>>> cosmo.hubble_time
<Quantity 13.96846031 Gyr>
>>> cosmo.w(0.5)
np.float64(-0.9520574461395797)

Reference/API
*************

.. py:currentmodule:: astropy.cosmology.traits

.. automodapi:: astropy.cosmology.traits
   :inherited-members:
