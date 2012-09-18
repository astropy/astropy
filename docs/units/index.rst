***********************
Units (`astropy.units`)
***********************

.. currentmodule:: astropy.units

Introduction
============

``astropy.units`` is a Python package to handle defining and converting
between physical units.

Unlike some other unit-related Python packages, `astropy.units` does
not aim to provide operations on unitized values.  Instead, it just
handles the unit description that must be associated with values by
some other means.

Getting Started
===============

  >>> from astropy import units as u
  >>> # Convert from parsec to meter
  >>> u.pc.to(u.m)
  3.0856776e+16
  >>> speed_unit = u.cm / u.s
  >>> speed_unit.to(u.mile / u.hour, 1)
  0.02236936292054402
  >>> speed_unit.to([1., 1000., 5000.])
  array([  2.23693629e-02,   2.23693629e+01,   1.11846815e+02])

`astropy.units` also handles equivalencies, such as that between
wavelength and frequency.  To use that feature, equivalence objects
are passed to the `to` conversion method::

  # Wavelength to frequency doesn't normally work
  >>> u.nm.to(u.Hz, [1000, 2000])
  UnitsException: 'nm' and 'Hz' are not convertible
  # ...but by passing an equivalency unit (sp()), it does...
  >>> u.nm.to(u.Hz, [1000, 2000], u.sp())
  array([  2.99792458e+14,   1.49896229e+14])
  >>> u.nm.to(u.eV, [1000, 2000], u.sp())
  array([ 1.23984201,  0.61992101])

Using `astropy.units`
=====================

.. toctree::
   :maxdepth: 2

   standard_units
   composing_and_defining
   conversion
   format
   equivalencies

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


Acknowledgments
===============

This code is adapted from the `pynbody
<http://code.google.com/p/pynbody/>`_ units module written by Andrew
Pontzen, who has granted the Astropy project permission to use the
code under a BSD license.
