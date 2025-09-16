.. _astropy-cosmology-units-and-equivalencies:

************************************
Cosmological Units and Equivalencies
************************************

.. currentmodule:: astropy.cosmology.units

This package defines and collects cosmological units and equivalencies.
We suggest importing this units package as

    >>> import astropy.cosmology.units as cu


To enable the main :mod:`astropy.units` to access these units when searching
for unit conversions and equivalencies, use
:func:`~astropy.units.add_enabled_units`.

    >>> import astropy.units as u
    >>> u.add_enabled_units(cu)  # doctest: +SKIP


About the Units
===============

.. doctest::
   :hide:

   >>> import astropy.units as u


.. _cosmological-redshift:

Cosmological Redshift and Dimensionless Equivalency
---------------------------------------------------

There are numerous measures of distance in cosmology -- luminosities, CMB
temperature, the universe's age, etc. -- but redshift is the principal measure
from which others are defined. In cosmology, distance measures are commonly
exasperating to follow in a derivation, because they are used interchangeably.
``astropy`` provides the ``redshift`` unit and associated equivalencies to
assist in these derivations and unify the distance measures.

Examples
^^^^^^^^

.. EXAMPLE START: Using redshift-dimensionless equivalency

To convert to or from dimensionless to "redshift" units:

    >>> import astropy.units as u
    >>> import astropy.cosmology.units as cu
    >>> z = 1100 * cu.redshift
    >>> z.to(u.dimensionless_unscaled, equivalencies=cu.dimensionless_redshift())
    <Quantity 1100.>

The equivalency works as part of a quantity with composite units

    >>> q = (2.7 * u.K) * z
    >>> q.to(u.K, equivalencies=cu.dimensionless_redshift())
    <Quantity 2970. K>

Since the redshift is not a true unit and is used so frequently, the
redshift / dimensionless equivalency is actually enabled by default.

    >>> z == 1100 * u.dimensionless_unscaled
    np.True_

    >>> q.to(u.K)
    <Quantity 2970. K>

To temporarily remove the equivalency and enforce unit strictness, use
:func:`astropy.units.set_enabled_equivalencies` as a context.

    >>> with u.set_enabled_equivalencies([]):
    ...     try:
    ...         z.to(u.dimensionless_unscaled)
    ...     except u.UnitConversionError:
    ...         print("equivalency disabled")
    equivalency disabled

.. EXAMPLE END


.. EXAMPLE START: Using `with_redshift` equivalency

The other redshift equivalency is `~astropy.cosmology.units.with_redshift`,
enabling redshift to be converted to other units, like CMB temperature:

    >>> from astropy.cosmology import WMAP9
    >>> z = 1100 * cu.redshift
    >>> z.to(u.K, cu.with_redshift(WMAP9))
    <Quantity 3000.225 K>

or the Hubble parameter:

    >>> z.to(u.km / u.s / u.Mpc, cu.with_redshift(WMAP9))  # doctest: +FLOAT_CMP
    <Quantity 1565637.40154275 km / (Mpc s)>

    >>> z.to(cu.littleh, cu.with_redshift(WMAP9))  # doctest: +FLOAT_CMP
    <Quantity 15656.37401543 littleh>

or a physical distance (comoving, lookback, or luminosity):

.. doctest-requires:: scipy

    >>> z.to(u.Mpc, cu.with_redshift(WMAP9, distance="luminosity"))  # doctest: +FLOAT_CMP
    <Quantity 15418438.76317008 Mpc>

These conversions are cosmology dependent, so if the cosmology changes,
so too will the conversions.

    >>> excosmo = WMAP9.clone(Tcmb0=3.0)
    >>> z.to(u.K, cu.with_redshift(excosmo))
    <Quantity 3303. K>

If no argument is given (or the argument is `None`), this equivalency assumes
the current default |Cosmology|:

    >>> z.to(u.K, cu.with_redshift())
    <Quantity 3000.7755 K>

To use this equivalency in a larger block of code:

    >>> with u.add_enabled_equivalencies(cu.with_redshift()):
    ...     # long derivation here
    ...     z.to(u.K)
    <Quantity 3000.7755 K>

.. EXAMPLE END


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
the ``H0`` from the current default |Cosmology|:

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
