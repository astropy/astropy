.. _astropy-cosmology-traits:

****************
Cosmology Traits
****************

.. currentmodule:: astropy.cosmology.traits

The :mod:`~astropy.cosmology.traits` module provides reusable components, called
:term:`traits <trait type>`, that encapsulate specific cosmological properties or
behaviors. For example, the :class:`~astropy.cosmology.traits.HubbleParameter` trait
provides the Hubble constant (``H0``) and related methods, while
:class:`~astropy.cosmology.traits.ScaleFactor` and
:class:`~astropy.cosmology.traits.TemperatureCMB` provide the scale factor and the
temperature or the CMB, respectively.

By combining these traits, you can easily construct custom cosmology classes with
precisely the features you need, without having to reimplement common functionality.


Here is an example of how to use the
:class:`~astropy.cosmology.traits.HubbleParameter`,
:class:`~astropy.cosmology.traits.ScaleFactor`, and
:class:`~astropy.cosmology.traits.TemperatureCMB` traits in a custom cosmology class:

>>> import astropy.units as u
>>> from astropy.cosmology.traits import HubbleParameter, ScaleFactor, TemperatureCMB
>>> from astropy.cosmology import Cosmology
>>>
>>> class CustomCosmology(Cosmology, HubbleParameter, ScaleFactor, TemperatureCMB):
...     def __init__(self, Om0, Ode0, H0=70, Tcmb0=2.725):
...         self.H0 = H0 << (u.km / u.s / u.Mpc)
...         self.Om0 = Om0
...         self.Ode0 = Ode0
...         self.Tcmb0 = u.Quantity(Tcmb0, "K")
...         super().__init__()
...
...     is_flat = False
...     # Additional custom methods and properties can be added here

>>> cosmo = CustomCosmology(H0=70, Om0=0.3, Ode0=0.7)
>>> cosmo.H0
<Quantity 70. km / (Mpc s)>
>>> cosmo.scale_factor(0)
<Quantity 1.>
>>> cosmo.Tcmb(1)
<Quantity 5.45 K>
>>> cosmo.hubble_time
<Quantity 13.96846031 Gyr>

Reference/API
*************

.. py:currentmodule:: astropy.cosmology.traits

.. automodapi:: astropy.cosmology.traits
   :inherited-members:
