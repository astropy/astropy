.. _astropy-cosmology-traits:

****************
Cosmology Traits
****************

.. currentmodule:: astropy.cosmology.traits

The :mod:`~astropy.cosmology.traits` module hosts various parts of cosmologies, such as
the :class:`~astropy.cosmology.traits.ScaleFactor` or
:class:`~astropy.cosmology.traits.TemperatureCMB`. These :term:`traits <trait type>` can
be used to more easily construct custom cosmologies by combining different components.

As a simple example, the :class:`~astropy.cosmology.traits.TemperatureCMB` trait
provides the ``Tcmb0`` property and
:meth:`~astropy.cosmology.traits.TemperatureCMB.Tcmb` method for computing the
cosmological CMB temperature at specified redshifts. By using this trait, you can add
temperature-related  functionality to your custom cosmology class without having to
implement it from scratch.

Here is an example of how to use the :class:`~astropy.cosmology.traits.ScaleFactor` and
:class:`~astropy.cosmology.traits.TemperatureCMB` traits in a custom cosmology class:

>>> import astropy.units as u
>>> from astropy.cosmology.traits import ScaleFactor, TemperatureCMB
>>> from astropy.cosmology import Cosmology
>>>
>>> class CustomCosmology(Cosmology, ScaleFactor, TemperatureCMB):
...     def __init__(self, H0, Om0, Ode0, Tcmb0=2.725):
...         self.H0 = H0
...         self.Om0 = Om0
...         self.Ode0 = Ode0
...         self.Tcmb0 = u.Quantity(Tcmb0, "K")
...         super().__init__()
...
...     is_flat = False
...     # Additional custom methods and properties can be added here

>>> cosmo = CustomCosmology(H0=70, Om0=0.3, Ode0=0.7)
>>> cosmo.scale_factor(0)
<Quantity 1.>
>>> cosmo.Tcmb(1)
<Quantity 5.45 K>

By combining different traits, you can create fully-featured cosmology classes with
minimal effort.


Reference/API
*************

.. py:currentmodule:: astropy.cosmology.traits

.. automodapi:: astropy.cosmology.traits
   :inherited-members:
