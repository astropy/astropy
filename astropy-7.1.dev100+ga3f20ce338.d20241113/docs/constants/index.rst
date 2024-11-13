.. _astropy-constants:

*******************************
Constants (`astropy.constants`)
*******************************

.. currentmodule:: astropy.constants

Introduction
============

`astropy.constants` contains a number of physical constants useful in
Astronomy. A `~astropy.constants.Constant` is a |Quantity| object with
additional metadata describing its provenance and uncertainty.

Getting Started
===============

You can import a :class:`~astropy.constants.Constant` directly from the
:mod:`astropy.constants` sub-package::

    >>> from astropy.constants import G
    >>> print(G)
      Name   = Gravitational constant
      Value  = 6.6743e-11
      Uncertainty  = 1.5e-15
      Unit  = m3 / (kg s2)
      Reference = CODATA 2018

Or, if you want to avoid having to explicitly import all of the constants you
need, you can do::

    >>> from astropy import constants as const
    >>> print(const.G)
      Name   = Gravitational constant
      ...

Constants can be used in :ref:`quantity_arithmetic` operations and
:ref:`quantity_and_numpy` just like any other |Quantity|::

    >>> from astropy import units as u
    >>> F = (const.G * 3. * const.M_sun * 100 * u.kg) / (2.2 * u.au) ** 2
    >>> print(F.to(u.N))  # doctest: +FLOAT_CMP
    0.3675671602160826 N

Unit Conversion
===============

..
  EXAMPLE START
  Converting Constants to Different Units

Explicitly :ref:`quantity_unit_conversion` is often not necessary, but can be
done if needed::

    >>> print(const.c)
      Name   = Speed of light in vacuum
      Value  = 299792458.0
      Uncertainty  = 0.0
      Unit  = m / s
      Reference = CODATA 2018

    >>> print(const.c.to('km/s'))
    299792.458 km / s

    >>> print(const.c.to('pc/yr'))  # doctest: +FLOAT_CMP
    0.306601393788 pc / yr

It is possible to convert most constants to `Centimeter-Gram-Second (CGS)
<https://en.wikipedia.org/wiki/Centimetre-gram-second_system_of_units>`_ units
using, for example::

    >>> const.c.cgs  # doctest: +FLOAT_CMP
    <Quantity   2.99792458e+10 cm / s>

However, some constants are defined with `different physical dimensions in CGS
<https://en.wikipedia.org/wiki/Centimetre-gram-second_system_of_units#Alternative_derivations_of_CGS_units_in_electromagnetism>`_
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

..
  EXAMPLE END

.. _astropy-constants-prior:

Collections of Constants (and Prior Versions)
=============================================

Constants are organized into version modules. The constants for
``astropy`` 2.0 can be accessed in the ``astropyconst20`` module.
For example::

    >>> from astropy.constants import astropyconst20 as const
    >>> print(const.e)
      Name   = Electron charge
      Value  = 1.6021766208e-19
      Uncertainty  = 9.8e-28
      Unit  = C
      Reference = CODATA 2014

The version modules contain physical and astronomical constants, and both sets
can also be chosen independently from each other. Physical `CODATA constants
<https://physics.nist.gov/cuu/Constants/index.html>`_ are in modules with names
like ``codata2010``, ``codata2014``, or ``codata2018``::

    >>> from astropy.constants import codata2014 as const
    >>> print(const.h)
      Name   = Planck constant
      Value  = 6.62607004e-34
      Uncertainty  = 8.1e-42
      Unit  = J s
      Reference = CODATA 2014

Astronomical constants defined (primarily) by the International Astronomical
Union (IAU) are collected in modules with names like ``iau2012`` or ``iau2015``::

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

However, importing these prior version modules directly will lead to
inconsistencies with other subpackages that have already imported
`astropy.constants`. Notably, `astropy.units` will have already used
the default version of constants. When using prior versions of the constants
in this manner, quantities should be constructed with constants instead of units.

To ensure consistent use of a prior version of constants in other ``astropy``
packages (such as :mod:`astropy.units`) that import :mod:`astropy.constants`,
the physical and astronomical constants versions should be set via
:class:`~astropy.utils.state.ScienceState` classes. These must be set before
the first import of either :mod:`astropy.constants` or :mod:`astropy.units`.
For example, you can use the CODATA2010 physical constants together with the
IAU 2012 astronomical constants::

    >>> from astropy import physical_constants, astronomical_constants
    >>> physical_constants.set('codata2010')  # doctest: +SKIP
    <ScienceState physical_constants: 'codata2010'>
    >>> physical_constants.get()  # doctest: +SKIP
    'codata2010'
    >>> astronomical_constants.set('iau2012')  # doctest: +SKIP
    <ScienceState astronomical_constants: 'iau2012'>
    >>> astronomical_constants.get()  # doctest: +SKIP
    'iau2012'

Then all other packages that import `astropy.constants` will self-consistently
initialize with these prior versions of constants.

The versions may also be set using values referring to the version modules::

    >>> from astropy import physical_constants, astronomical_constants
    >>> physical_constants.set('astropyconst13')  # doctest: +SKIP
    <ScienceState physical_constants: 'codata2010'>
    >>> physical_constants.get()  # doctest: +SKIP
    'codata2010'
    >>> astronomical_constants.set('astropyconst13')  # doctest: +SKIP
    <ScienceState astronomical_constants: 'iau2012'>
    >>> astronomical_constants.get()  # doctest: +SKIP
    'iau2012'

.. The doctest should not be skipped, ideally. See https://github.com/astropy/astropy/issues/8781

If :mod:`astropy.constants` or :mod:`astropy.units` have already been imported,
a :class:`RuntimeError` will be raised::

    >>> import astropy.units
    >>> from astropy import physical_constants, astronomical_constants
    >>> astronomical_constants.set('astropyconst13')
    Traceback (most recent call last):
        ...
    RuntimeError: astropy.units is already imported

.. note that if this section gets too long, it should be moved to a separate
   doc page - see the top of performance.inc.rst for the instructions on how to
   do that
.. include:: performance.inc.rst

Reference/API
=============

.. automodapi:: astropy.constants
