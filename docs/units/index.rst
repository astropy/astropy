.. _astropy-units:

**************************************
Units and Quantities (`astropy.units`)
**************************************

.. |quantity| replace:: :class:`~astropy.units.Quantity`

.. currentmodule:: astropy.units

Introduction
============

`astropy.units` handles defining, converting between, and performing
arithmetic with physical quantities, such as meters, seconds, Hz,
etc.  It also handles logarithmic units such as magnitude and decibel.

`astropy.units` does not know spherical geometry or sexagesimal
(hours, min, sec): if you want to deal with celestial coordinates,
see the `astropy.coordinates` package.

Getting Started
===============

Most users of the `astropy.units` package will :ref:`work with "quantities"
<quantity>`: the combination of a value and a unit.  The easiest way to create
a |quantity| is to simply multiply or divide a value by one of the built-in
units.  It works with scalars, sequences and Numpy arrays::

    >>> from astropy import units as u
    >>> 42.0 * u.meter
    <Quantity 42.0 m>
    >>> [1., 2., 3.] * u.m
    <Quantity [ 1., 2., 3.] m>
    >>> import numpy as np
    >>> np.array([1., 2., 3.]) * u.m
    <Quantity [ 1., 2., 3.] m>

You can get the unit and value from a |quantity| using the unit and
value members::

    >>> q = 42.0 * u.meter
    >>> q.value
    42.0
    >>> q.unit
    Unit("m")

From this simple building block, it's easy to start combining
quantities with different units::

    >>> 15.1 * u.meter / (32.0 * u.second)  # doctest: +FLOAT_CMP
    <Quantity 0.471875 m / s>
    >>> 3.0 * u.kilometer / (130.51 * u.meter / u.second)  # doctest: +FLOAT_CMP
    <Quantity 0.022986744310780783 km s / m>
    >>> (3.0 * u.kilometer / (130.51 * u.meter / u.second)).decompose()  # doctest: +FLOAT_CMP
    <Quantity 22.986744310780782 s>

Unit conversion is done using the
:meth:`~astropy.units.quantity.Quantity.to` method, which returns a new
|quantity| in the given unit::

    >>> x = 1.0 * u.parsec
    >>> x.to(u.km)  # doctest: +FLOAT_CMP
    <Quantity 30856775814671.914 km>

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

`astropy.units` is able to match compound units against the units it already
knows about::

    >>> (u.s ** -1).compose()  # doctest: +SKIP
    [Unit("Bq"), Unit("Hz"), Unit("3.7e+10 Ci")]

And it can convert between unit systems, such as SI or CGS:

.. doctest-skip::

    >>> (1.0 * u.Pa).cgs
    <Quantity 10.0 Ba>

The units ``mag``, ``dex`` and ``dB`` are special, being :ref:`logarithmic
units <logarithmic_units>`, for which a value is the logarithm of a physical
quantity in a given unit.  These can be used with a physical unit in
parentheses to create a corresponding logarithmic quantity::

    >>> -2.5 * u.mag(u.ct / u.s)
    <Magnitude -2.5 mag(ct / s)>
    >>> from astropy import constants as c
    >>> u.Dex((c.G * u.M_sun / u.R_sun**2).cgs)  # doctest: +FLOAT_CMP
    <Dex 4.438443765928838 dex(cm / s2)>
    
`astropy.units` also handles :ref:`equivalencies <unit_equivalencies>`, such as
that between wavelength and frequency. To use that feature, equivalence objects
are passed to the :meth:`~astropy.units.quantity.Quantity.to` conversion
method. For instance, a conversion from wavelength to frequency doesn't
normally work:

    >>> (1000 * u.nm).to(u.Hz)
    Traceback (most recent call last):
      ...
    UnitConversionError: 'nm' (length) and 'Hz' (frequency) are not convertible

but by passing an equivalency list, in this case ``spectral()``, it does:

    >>> (1000 * u.nm).to(u.Hz, equivalencies=u.spectral())
    <Quantity 299792457999999.94 Hz>

Quantities and units can be :ref:`printed nicely to strings
<astropy-units-format>` using the `Format String Syntax
<http://docs.python.org/library/string.html#format-string-syntax>`_, the
preferred string formatting syntax in recent versions of python.  Format
specifiers (like ``0.03f``) in new-style format strings will used to format the
quantity value::

    >>> q = 15.1 * u.meter / (32.0 * u.second)
    >>> q  # doctest: +FLOAT_CMP
    <Quantity 0.471875 m / s>
    >>> "{0:0.03f}".format(q)
    '0.472 m / s'

The value and unit can also be formatted separately. Format specifiers
used on units can be used to choose the unit formatter::

    >>> q = 15.1 * u.meter / (32.0 * u.second)
    >>> q  # doctest: +FLOAT_CMP
    <Quantity 0.471875 m / s>
    >>> "{0.value:0.03f} {0.unit:FITS}".format(q)
    '0.472 m s-1'

Using `astropy.units`
=====================

.. toctree::
   :maxdepth: 2

   quantity
   standard_units
   combining_and_defining
   decomposing_and_composing
   logarithmic_units
   format
   equivalencies
   conversion

See Also
========

- `FITS Standard <http://fits.gsfc.nasa.gov/fits_standard.html>`_ for
  units in FITS.

- The `Units in the VO 1.0 Standard
  <http://www.ivoa.net/Documents/VOUnits/>`_ for representing units in
  the VO.

- OGIP Units: A standard for storing units in `OGIP FITS files
  <http://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/general/ogip_93_001/>`_.

- `Standards for astronomical catalogues units
  <http://cds.u-strasbg.fr/doc/catstd-3.2.htx>`_.

- `IAU Style Manual
  <http://www.iau.org/static/publications/stylemanual1989.pdf>`_.

- `A table of astronomical unit equivalencies
  <http://astro.wku.edu/strolger/UNITS.txt>`_

Reference/API
=============

.. automodapi:: astropy.units.quantity

.. automodapi:: astropy.units

.. automodapi:: astropy.units.format

.. automodapi:: astropy.units.si

.. automodapi:: astropy.units.cgs

.. automodapi:: astropy.units.astrophys

.. automodapi:: astropy.units.function.units

.. automodapi:: astropy.units.imperial

.. automodapi:: astropy.units.cds

.. automodapi:: astropy.units.equivalencies

.. automodapi:: astropy.units.function

Acknowledgments
===============

This code is adapted from the `pynbody
<https://github.com/pynbody/pynbody>`__ units module written by Andrew
Pontzen, who has granted the Astropy project permission to use the code
under a BSD license.
