Using Prior Versions of Constants
*********************************

By default, `astropy.units` are initialized upon first import to use
the current versions of `astropy.constants`. For units to initialize
properly to a prior version of constants, the constants versions must
be set before the first import of `astropy.units` or `astropy.constants`.

This is accomplished using :class:`~astropy.utils.state.ScienceState` classes
in the top-level package. Setting the prior versions at the start of a Python
session will allow consistent units.

Example
=======

.. EXAMPLE START: Using Prior Versions of Constants

To initialize units to a prior version of constants:

>>> import astropy
>>> astropy.physical_constants.set('codata2010')  # doctest: +SKIP
<ScienceState physical_constants: 'codata2010'>
>>> astropy.astronomical_constants.set('iau2012')  # doctest: +SKIP
<ScienceState astronomical_constants: 'iau2012'>
>>> import astropy.units as u
>>> import astropy.constants as const
>>> (const.M_sun / u.M_sun).to(u.dimensionless_unscaled) - 1  # doctest: +SKIP
<Quantity 0.>
>>> print(const.M_sun)  # doctest: +SKIP
  Name   = Solar mass
  Value  = 1.9891e+30
  Uncertainty  = 5e+25
  Unit  = kg
  Reference = Allen's Astrophysical Quantities 4th Ed.

If :mod:`astropy.units` has already been imported, a :class:`RuntimeError` is
raised.

.. EXAMPLE END
