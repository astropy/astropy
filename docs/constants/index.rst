*******************************
Constants (`astropy.constants`)
*******************************

.. currentmodule:: astropy.constants

Introduction
============

`astropy.constants` contains a number of physical constants useful in
Astronomy. Constants are `~astropy.units.quantity.Quantity` objects with
additional meta-data describing their provenance and uncertainties.

Getting Started
===============

To use the constants in S.I. units, you can import the constants directly from
the `astropy.constants.si` sub-package::

    >>> from astropy.constants.si import G

or, if you want to avoid having to explicitly import all the constants you
need, you can simply do:

    >>> from astropy.constants import si

and then subsequently use for example ``si.G``. Constants are fully-fleged
`~astropy.units.quantity.Quantity` objects, so you can easily convert them to
different units for example::

    >>> print si.c
      Name   = Speed of light in vacuum
      Value  = 299792458.0
      Error  = 0.0
      Units = m / (s)
      Reference = CODATA 2010

    >>> print si.c.to('km/s')
    299792.458 km / (s)

    >>> print si.c.to('pc/yr')
    0.306594845466 pc / (yr)

and you can use them in conjunction with unit and other non-constant
`~astropy.units.quantity.Quantity` objects::

    >>> F = (si.G * 3. * si.M_sun * 100 * u.kg) / (2.2 * u.au) ** 2
    >>> print F.to(u.N)
    0.367669392028 N

While it is possible to convert most constants to c.g.s using e.g.::

    >>> si.c.cgs
    <Quantity 29979245800.0 cm / (s)>

it is also possible to simply import the constants directly in the c.g.s
system::

    >>> from astropy.constants import cgs

    >>> print cgs.c
      Name   = Speed of light in vacuum
      Value  = 29979245800.0
      Error  = 0.0
      Units = cm / (s)
      Reference = CODATA 2010


Reference/API
=============

.. automodapi:: astropy.constants.si

.. automodapi:: astropy.constants.cgs
