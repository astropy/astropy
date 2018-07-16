.. _astropy-constants:

*******************************
Constants (`astropy.constants`)
*******************************

.. currentmodule:: astropy.constants

Introduction
============

`astropy.constants` contains a number of physical constants useful in
Astronomy. Constants are `~astropy.units.Quantity` objects with
additional meta-data describing their provenance and uncertainties.

Getting Started
===============

To use the constants in S.I. units, you can import the constants directly from
the `astropy.constants` sub-package::

    >>> from astropy.constants import G

or, if you want to avoid having to explicitly import all the constants you
need, you can simply do:

    >>> from astropy import constants as const

and then subsequently use for example ``const.G``. Constants are fully-fledged
`~astropy.units.Quantity` objects, so you can easily convert them to
different units for example::

    >>> print(const.c)
      Name   = Speed of light in vacuum
      Value  = 299792458.0
      Uncertainty  = 0.0
      Unit  = m / s
      Reference = CODATA 2014

    >>> print(const.c.to('km/s'))
    299792.458 km / s

    >>> print(const.c.to('pc/yr'))  # doctest: +FLOAT_CMP
    0.306601393788 pc / yr

and you can use them in conjunction with unit and other non-constant
`~astropy.units.Quantity` objects::

    >>> from astropy import units as u
    >>> F = (const.G * 3. * const.M_sun * 100 * u.kg) / (2.2 * u.au) ** 2
    >>> print(F.to(u.N))  # doctest: +FLOAT_CMP
    0.3675671602160826 N

It is possible to convert most constants to cgs using e.g.::

    >>> const.c.cgs  # doctest: +FLOAT_CMP
    <Quantity   2.99792458e+10 cm / s>

However, some constants are defined with different physical dimensions in cgs
and cannot be directly converted. Because of this ambiguity, such constants
cannot be used in expressions without specifying a system::

    >>> 100 * const.e
    Traceback (most recent call last):
        ...
    TypeError: Constant u'e' does not have physically compatible units
    across all systems of units and cannot be combined with other
    values without specifying a system (eg. e.emu)
    >>> 100 * const.e.esu  # doctest: +FLOAT_CMP
    <Quantity 4.8032045057134676e-08 Fr>

.. _astropy-constants-prior:

Collections of constants (and prior versions)
=============================================

Constants are organized into version modules. The constants for
Astropy 1.3 can be accessed in the ``astropyconst13`` module.
For example:

    >>> from astropy.constants import astropyconst13 as const
    >>> print(const.e)
      Name   = Electron charge
      Value  = 1.602176565e-19
      Uncertainty  = 3.5e-27
      Unit  = C
      Reference = CODATA 2010

Physical CODATA constants are in modules with names like ``codata2010`` or
``codata2014``:

    >>> from astropy.constants import codata2010 as const
    >>> print(const.h)
      Name   = Planck constant
      Value  = 6.62606957e-34
      Uncertainty  = 2.9e-41
      Unit  = J s
      Reference = CODATA 2010

Astronomical constants defined (primarily) by the IAU are collected in
modules with names like ``iau2012`` or ``iau2015``:

    >>> from astropy.constants import iau2012 as const
    >>> print(const.L_sun)
      Name   = Solar luminosity
      Value  = 3.846e+26
      Uncertainty  = 5e+22
      Unit  = W
      Reference = Allen's Astrophysical Quantities 4th Ed.

    >>> from astropy.constants import iau2015 as const
    >>> print(const.L_sun)
      Name   = Nominal solar luminosity
      Value  = 3.828e+26
      Uncertainty  = 0.0
      Unit  = W
      Reference = IAU 2015 Resolution B 3

The astronomical and physical constants are combined into modules with
names like ``astropyconst13`` and ``astropyconst20``. To temporarily set
constants to an older version (e.g., for regression testing), a context
manager is available, as follows:

    >>> from astropy import constants as const
    >>> with const.set_enabled_constants('astropyconst13'):
    ...     print(const.h)
      Name   = Planck constant
      Value  = 6.62606957e-34
      Uncertainty  = 2.9e-41
      Unit  = J s
      Reference = CODATA 2010
    >>> print(const.h)
      Name   = Planck constant
      Value  = 6.62607004e-34
      Uncertainty  = 8.1e-42
      Unit  = J s
      Reference = CODATA 2014

.. warning::

    Units such as ``u.M_sun`` will use the current version of the
    corresponding constant. When using prior versions of the constants,
    quantities should be constructed with constants instead of units.

.. note that if this section gets too long, it should be moved to a separate
   doc page - see the top of performance.inc.rst for the instructions on how to do
   that
.. include:: performance.inc.rst

Reference/API
=============

.. automodapi:: astropy.constants
