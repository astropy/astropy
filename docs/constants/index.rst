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

    >>> print const.c
      Name   = Speed of light in vacuum
      Value  = 299792458.0
      Uncertainty  = 0.0
      Unit  = m / s
      Reference = CODATA 2010

    >>> print const.c.to('km/s')
    299792.458 km / s

    >>> print const.c.to('pc/yr')
    0.306601393788 pc / yr

and you can use them in conjunction with unit and other non-constant
`~astropy.units.Quantity` objects::

    >>> from astropy import units as u
    >>> F = (const.G * 3. * const.M_sun * 100 * u.kg) / (2.2 * u.au) ** 2
    >>> print F.to(u.N)
    0.367669392028 N

It is possible to convert most constants to cgs using e.g.::

    >>> const.c.cgs
    <Quantity 29979245800.0 cm / s>

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

Reference/API
=============

.. automodapi:: astropy.constants
