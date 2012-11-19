***********************
Units (`astropy.units`)
***********************

.. currentmodule:: astropy.units

Introduction
============

``astropy.units`` is a Python package to handle defining and converting
between physical units, and performing arithmetic with physical quantities
(numbers with associated units).

Getting Started
===============

  >>> from astropy import units as u
  >>> # Convert from parsec to meter
  >>> u.pc.to(u.m)
  3.0856776e+16
  >>> cms = u.cm / u.s
  >>> mph = u.mile / u.hour
  >>> cms.to(mph, 1)
  0.02236936292054402
  >>> cms.to(mph, [1., 1000., 5000.])
  array([  2.23693629e-02,   2.23693629e+01,   1.11846815e+02])

Units that "cancel out" become a special unit called the
"dimensionless unit":

  >>> u.m / u.m
  Unit(dimensionless)

`astropy.units` also handles equivalencies, such as that between
wavelength and frequency.  To use that feature, equivalence objects
are passed to the `~astropy.units.core.UnitBase.to` conversion method::

  # Wavelength to frequency doesn't normally work
  >>> u.nm.to(u.Hz, [1000, 2000])
  UnitsException: 'nm' (length) and 'Hz' (frequency) are not convertible
  # ...but by passing an equivalency unit (spectral()), it does...
  >>> u.nm.to(u.Hz, [1000, 2000], equivs=u.spectral())
  array([  2.99792458e+14,   1.49896229e+14])
  >>> u.nm.to(u.eV, [1000, 2000], equivs=u.spectral())
  array([ 1.23984201,  0.61992101])

Also included in the `astropy.units` package is the `Quantity` object,
which represents a numerical value with an associated unit. These objects
support arithmetic with other numbers and `Quantity` objects and preserve
units.
   >>> from astropy import units as u
   >>> 15.1*u.meter / (32.0*u.second)
   <Quantity 0.471875 m / (s)>
   >>> 3.0*u.kilometer / (130.51*u.meter/u.second)
   <Quantity 0.0229867443108 km s / (m)>
   >>> (3.0*u.kilometer / (130.51*u.meter/u.second)).simplify_units()
   <Quantity 22.9867443108 s>

Using `astropy.units`
=====================

.. toctree::
   :maxdepth: 2

   standard_units
   composing_and_defining
   conversion
   format
   equivalencies
   quantity

See Also
========

- `FITS Standard <http://fits.gsfc.nasa.gov/fits_standard.html>`_ for
  units in FITS.

- The `proposed IVOA standard
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

.. automodapi:: astropy.units.core

.. automodapi:: astropy.units.format

.. automodapi:: astropy.units.si

.. automodapi:: astropy.units.cgs

.. automodapi:: astropy.units.astrophys

.. automodapi:: astropy.units.imperial

.. automodapi:: astropy.units.equivalencies

.. automodapi:: astropy.units.quantity


Acknowledgments
===============

This code is adapted from the `pynbody
<http://code.google.com/p/pynbody/>`_ units module written by Andrew
Pontzen, who has granted the Astropy project permission to use the
code under a BSD license.
