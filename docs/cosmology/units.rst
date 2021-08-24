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


About the Units
===============

.. _littleh-and-H0-equivalency:

Reduced Hubble Constant and "little-h" Equivalency
--------------------------------------------------

The dimensionless version of the Hubble constant — often known as "little h" —
is a frequently used quantity in extragalactic astrophysics. It is also widely
known as the bane of beginners' existence in such fields (See e.g., the title
of `this paper <https://doi.org/10.1017/pasa.2013.31>`__, which also provides
valuable advice on the use of little h). ``astropy`` provides the
:func:`~astropy.cosmology.units.with_H0` equivalency that helps keep this
straight in at least some of these cases, by providing a way to convert to/from
physical to "little h" units.

Examples
^^^^^^^^

.. EXAMPLE START: Using the "little h" Equivalency

To convert to or from physical to "little h" units:

.. code-block:: python

    >>> import astropy.units as u
    >>> import astropy.cosmology.units as cu
    >>> H0_70 = 70 * u.km/u.s/u.Mpc
    >>> distance = 70 * (u.Mpc/cu.littleh)
    >>> distance.to(u.Mpc, cu.with_H0(H0_70))  # doctest: +FLOAT_CMP
    <Quantity 100.0 Mpc>
    >>> luminosity = 0.49 * u.Lsun * cu.littleh**-2
    >>> luminosity.to(u.Lsun, cu.with_H0(H0_70))  # doctest: +FLOAT_CMP
    <Quantity 1.0 solLum>

Note the unit name ``littleh``: while this unit is usually expressed in the
literature as just ``h``, here it is ``littleh`` to avoid confusion with
"hours."

If no argument is given (or the argument is `None`), this equivalency assumes
the ``H0`` from the current default :ref:`cosmology <astropy-cosmology>`:

.. code-block:: python

    >>> distance = 100 * (u.Mpc/cu.littleh)
    >>> distance.to(u.Mpc, cu.with_H0())  # doctest: +FLOAT_CMP
    <Quantity 147.79781259 Mpc>

This equivalency also allows a common magnitude formulation of little h
scaling:

.. code-block:: python

    >>> mag_quantity = 12 * (u.mag - u.MagUnit(cu.littleh**2))
    >>> mag_quantity  # doctest: +FLOAT_CMP
    <Magnitude 12. mag(1 / littleh2)>
    >>> mag_quantity.to(u.mag, cu.with_H0(H0_70))  # doctest: +FLOAT_CMP
    <Quantity 11.2254902 mag>

.. EXAMPLE END


Reference/API
=============

.. automodapi:: astropy.cosmology.units
   :inherited-members:
