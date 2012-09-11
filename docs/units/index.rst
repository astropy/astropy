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
  >>> u.pc.to(u.m, 1)
  3.0856776e+16
  >>> speed_unit = u.cm / u.s
  >>> speed_unit.to(u.mile / u.hour, 1)
  0.02236936292054402
  >>> speed_converter = speed_unit.get_converter("mile hour^-1")
  >>> speed_converter([1., 1000., 5000.])
  array([  2.23693629e-02,   2.23693629e+01,   1.11846815e+02])

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

.. automodapi:: astropy.units.standard_units

Acknowledgments
===============

This code is adapted from the `pynbody
<http://code.google.com/p/pynbody/>`_ units module written by Andrew
Pontzen, who has granted the Astropy project permission to use the
code under a BSD license.
