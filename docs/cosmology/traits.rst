.. _astropy-cosmology-traits:

****************
Cosmology Traits
****************

.. currentmodule:: astropy.cosmology.traits

The :mod:`~astropy.cosmology.traits` module provides reusable components, called
:term:`traits <trait type>`, that encapsulate specific cosmological properties or
behaviors. For example, the :class:`~astropy.cosmology.traits.NeutrinoComponent` trait
provides the neutrino component and related methods.
By combining these traits, you can easily construct custom cosmology classes with
precisely the features you need, without having to reimplement common functionality.
Here is an example of how to use the
:class:`~astropy.cosmology.traits.NeutrinoComponent` trait in custom cosmology classes:

>>> import numpy as np
>>> from astropy import units as u
>>> from astropy.cosmology import Cosmology
>>> from astropy.cosmology.traits import TemperatureCMB, NeutrinoComponent
>>> NEUTRINO_FERMI_DIRAC_CORRECTION = 0.22710731766  # 7/8 (4/11)^4/3
>>>
>>> class SimpleNeutrinoCosmology(Cosmology, TemperatureCMB, NeutrinoComponent):
...     '''Minimal cosmology demonstrating NeutrinoComponent trait.'''
...     def __init__(self, Tcmb0=2.7255, Neff=3.046, Ogamma0=5e-5):
...         self.Tcmb0 = u.Quantity(Tcmb0, "K")
...         self.Neff = Neff
...         self.Ogamma0 = Ogamma0
...         super().__init__()
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
...     @property
...     def is_flat(self):
...         return True
>>>
>>> cosmo = SimpleNeutrinoCosmology()
>>> cosmo.Tnu0
<Quantity 1.94... K>
>>> cosmo.has_massive_nu
False
>>> cosmo.Onu0
3.4...e-05
>>> cosmo.Onu(0)
np.float64(3.4...e-05)
>>> cosmo.Tnu(1)
<Quantity 3.89... K>

Reference/API
*************

.. py:currentmodule:: astropy.cosmology.traits

.. automodapi:: astropy.cosmology.traits
   :inherited-members:
