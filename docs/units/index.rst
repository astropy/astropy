.. _astropy-units:

**************************************
Units and Quantities (`astropy.units`)
**************************************

.. currentmodule:: astropy.units

Introduction
============

`astropy.units` handles defining, converting between, and performing
arithmetic with physical quantities, such as meters, seconds, Hz,
etc. It also handles logarithmic units such as magnitude and decibel.

`astropy.units` does not know spherical geometry or sexagesimal
(hours, min, sec): if you want to deal with celestial coordinates,
see the `astropy.coordinates` package.

Getting Started
===============

Most users of the `astropy.units` package will work with :ref:`Quantity objects
<quantity>`: the combination of a value and a unit. The most convenient way to
create a |Quantity| is to multiply or divide a value by one of the built-in
units. It works with scalars, sequences, and ``numpy`` arrays.

Examples
--------

.. EXAMPLE START: Creating and Combining Quantities with Units

To create a |Quantity| object::

    >>> from astropy import units as u
    >>> 42.0 * u.meter  # doctest: +FLOAT_CMP
    <Quantity  42. m>
    >>> [1., 2., 3.] * u.m  # doctest: +FLOAT_CMP
    <Quantity [1., 2., 3.] m>
    >>> import numpy as np
    >>> np.array([1., 2., 3.]) * u.m  # doctest: +FLOAT_CMP
    <Quantity [1., 2., 3.] m>

You can get the unit and value from a |Quantity| using the unit and
value members::

    >>> q = 42.0 * u.meter
    >>> q.value
    np.float64(42.0)
    >>> q.unit
    Unit("m")

From this basic building block, it is possible to start combining
quantities with different units::

    >>> 15.1 * u.meter / (32.0 * u.second)  # doctest: +FLOAT_CMP
    <Quantity 0.471875 m / s>
    >>> 3.0 * u.kilometer / (130.51 * u.meter / u.second)  # doctest: +FLOAT_CMP
    <Quantity 0.022986744310780783 km s / m>
    >>> (3.0 * u.kilometer / (130.51 * u.meter / u.second)).decompose()  # doctest: +FLOAT_CMP
    <Quantity 22.986744310780782 s>

Unit conversion is done using the
:meth:`~astropy.units.quantity.Quantity.to` method, which returns a new
|Quantity| in the given unit::

    >>> x = 1.0 * u.parsec
    >>> x.to(u.km)  # doctest: +FLOAT_CMP
    <Quantity 30856775814671.914 km>

.. EXAMPLE END

.. EXAMPLE START: Creating Custom Units for Quantity Objects

It is also possible to work directly with units at a lower level, for
example, to create custom units::

    >>> from astropy.units import imperial

    >>> cms = u.cm / u.s
    >>> # ...and then use some imperial units
    >>> mph = imperial.mile / u.hour

    >>> # And do some conversions
    >>> q = 42.0 * cms
    >>> q.to(mph)  # doctest: +FLOAT_CMP
    <Quantity 0.939513242662849 mi / h>

Units that "cancel out" become a special unit called the
"dimensionless unit":

    >>> u.m / u.m
    Unit(dimensionless)

To create a basic :ref:`dimensionless quantity <doc_dimensionless_unit>`,
multiply a value by the unscaled dimensionless unit::

    >>> q = 1.0 * u.dimensionless_unscaled
    >>> q.unit
    Unit(dimensionless)

.. EXAMPLE END

.. EXAMPLE START: Matching and Converting Between Units

`astropy.units` is able to match compound units against the units it already
knows about::

    >>> (u.s ** -1).compose()  # doctest: +SKIP
    [Unit("Bq"), Unit("Hz"), Unit("2.7027e-11 Ci")]

And it can convert between unit systems, such as `SI
<https://www.bipm.org/documents/20126/41483022/SI-Brochure-9-EN.pdf>`_ or `CGS
<https://en.wikipedia.org/wiki/Centimetre-gram-second_system_of_units>`_::

    >>> (1.0 * u.Pa).cgs
    <Quantity 10. P / s>

The units ``mag``, ``dex``, and ``dB`` are special, being :ref:`logarithmic
units <logarithmic_units>`, for which a value is the logarithm of a physical
quantity in a given unit. These can be used with a physical unit in
parentheses to create a corresponding logarithmic quantity::

    >>> -2.5 * u.mag(u.ct / u.s)
    <Magnitude -2.5 mag(ct / s)>
    >>> from astropy import constants as c
    >>> u.Dex((c.G * u.M_sun / u.R_sun**2).cgs)  # doctest: +FLOAT_CMP
    <Dex 4.438067627303133 dex(cm / s2)>

`astropy.units` also handles :ref:`equivalencies <unit_equivalencies>`, such as
that between wavelength and frequency. To use that feature, equivalence objects
are passed to the :meth:`~astropy.units.quantity.Quantity.to` conversion
method. For instance, a conversion from wavelength to frequency does not
normally work:

    >>> (1000 * u.nm).to(u.Hz)  # doctest: +IGNORE_EXCEPTION_DETAIL
    Traceback (most recent call last):
      ...
    UnitConversionError: 'nm' (length) and 'Hz' (frequency) are not convertible

But by passing an equivalency list, in this case
:func:`~astropy.units.equivalencies.spectral`, it does:

    >>> (1000 * u.nm).to(u.Hz, equivalencies=u.spectral())  # doctest: +FLOAT_CMP
    <Quantity  2.99792458e+14 Hz>

.. EXAMPLE END

.. EXAMPLE START: Printing Quantities and Units to Strings

Quantities and units can be :ref:`printed nicely to strings
<astropy-units-format>` using the `Format String Syntax
<https://docs.python.org/3/library/string.html#format-string-syntax>`_. Format
specifiers (like ``0.03f``) in strings will be used to format the quantity
value::

    >>> q = 15.1 * u.meter / (32.0 * u.second)
    >>> q  # doctest: +FLOAT_CMP
    <Quantity 0.471875 m / s>
    >>> f"{q:0.03f}"
    '0.472 m / s'

The value and unit can also be formatted separately. Format specifiers
for units can be used to choose the unit formatter::

    >>> q = 15.1 * u.meter / (32.0 * u.second)
    >>> q  # doctest: +FLOAT_CMP
    <Quantity 0.471875 m / s>
    >>> f"{q.value:0.03f} {q.unit:FITS}"
    '0.472 m s-1'

.. EXAMPLE END

Using `astropy.units`
=====================

.. toctree::
   :maxdepth: 2

   quantity
   type_hints
   standard_units
   combining_and_defining
   decomposing_and_composing
   logarithmic_units
   structured_units
   format
   equivalencies
   physical_types
   constants_versions
   conversion

Acknowledgments
===============

This code was originally based on the `pynbody
<https://github.com/pynbody/pynbody>`__ units module written by Andrew
Pontzen, who has granted the Astropy Project permission to use the code
under a BSD license.

See Also
========

- `FITS Standard <https://fits.gsfc.nasa.gov/fits_standard.html>`_ for
  units in FITS.

- The `Units in the VO 1.0 Standard
  <http://www.ivoa.net/documents/VOUnits/>`_ for representing units in
  the VO.

- OGIP Units: A standard for storing units in `OGIP FITS files
  <https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/general/ogip_93_001/>`_.

- `Standards for astronomical catalogues: units
  <https://vizier.unistra.fr/vizier/doc/catstd-3.2.htx>`_.

- `IAU Style Manual
  <https://www.iau.org/static/publications/stylemanual1989.pdf>`_.

- `A table of astronomical unit equivalencies
  <https://www.stsci.edu/~strolger/docs/UNITS.txt>`_.

.. note that if this section gets too long, it should be moved to a separate
   doc page - see the top of performance.inc.rst for the instructions on how to do
   that
.. include:: performance.inc.rst

Reference/API
=============

.. toctree::
   :maxdepth: 2

   ref_api
