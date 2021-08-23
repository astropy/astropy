.. _astropy-cosmology-units-and-equivalencies:

************************************
Cosmological Units and Equivalencies
************************************

.. currentmodule:: astropy.cosmology.units

This package defines and collects cosmological units and equivalencies.
Many of the units and equivalencies are also available in the
:mod:`astropy.units` namespace.

We suggest importing this units package as

    >>> import astropy.cosmology.units as cu


To enable the main :mod:`astropy.units` to access these units when searching
for unit conversions and equivalencies, use
:func:`~astropy.units.add_enabled_units`.

    >>> u.add_enabled_units(cu)  # doctest: +SKIP


Reference/API
=============

.. automodapi:: astropy.cosmology.units
   :inherited-members:
   :include-all-objects:
   :allowed-package-names: astropy.units
