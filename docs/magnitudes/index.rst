.. _astropy_magnitudes:

*********************************
Magnitudes (`astropy.magnitudes`)
*********************************

Introduction
============

The `~astropy.magnitudes` package provides classes for representing magnitude,
as well as tools for converting between magnitude and flux.

.. warning::

    This package is currently a work-in-progress, and thus it is possible
    there will be significant API changes in later versions of Astropy.


Getting Started
===============

Define a generic magnitude object and convert it to flux:

    >>> from astropy import magnitudes as m
    >>> mag = m.Magnitude(-5)
    >>> mag
    <Magnitude -5 mag>
    >>> mag.to_flux()
    100.0


Using `astropy.magnitudes`
==========================

Conversion between magnitude and flux, using zeropoint in the unit of magnitude:

.. math::

    mag = -2.5 \; log_{10} flux + zeropoint

    flux = 10^{-0.4 \; (mag - zeropoint)}


Generic Magnitude System
------------------------

In the generic magnitude system represented by
`~astropy.magnitudes.mags.Magnitude`, the intrinsic zeropoint is 0 mag.
However, for flexibility, its
:func:`~astropy.magnitudes.mags.Magnitude.to_flux` method can take any
user-defined zeropoint (in mag) for conversion to flux.

During initialization, if a `~astropy.units.quantity.Quantity` is given and
its unit is not a magnitude, its value is assumed to be some kind of generic
flux and converted to magnitude as such.

Examples
^^^^^^^^

>>> from astropy import units as u
>>> from astropy import magnitudes as m
>>> mag = m.Magnitude([-5, u.Quantity(100)])
>>> mag
<Magnitude [-5.,-5.] mag>
>>> mag.to_flux()
array([ 100.,  100.])
>>> mag.to_flux(zeropoint=-5)
array([ 1.,  1.])
>>> mag + u.Quantity(5.0, u.mag)
<Quantity [ 0., 0.] mag>


STMAG and ABMAG
---------------

The `~astropy.magnitudes.mags.STMAG` and `~astropy.magnitudes.mags.ABMAG`
systems are primarily used in Hubble Space Telescope (HST) calibration.
They are defined such that an object with a specific constant flux distribution
at all wavelengths will have zero color at all wavelengths, as shown in the
following table:

====== =================================================================== ===============
System Constant flux distribution                                          Zeropoint (mag)
====== =================================================================== ===============
STMAG  :math:`3.63 \times 10^{-9} \; erg \; cm^{-2} \; s^{-1} \; \AA^{-1}` -21.1
ABMAG  :math:`3.63 \times 10^{-20} \; erg \; cm^{-2} \; s^{-1} \; Hz^{-1}` -48.6
====== =================================================================== ===============

Examples
^^^^^^^^

>>> from astropy import units as u
>>> from astropy import magnitudes as m

Using STMAG:

>>> mag = m.STMAG(u.Quantity(3.63e-9, u.erg / u.cm ** 2 / u.s / u.AA))
>>> mag
<STMAG 0.00023343740971526472 mag>
>>> mag.to_flux()
3.6300000000000067e-09

Using ABMAG:

>>> mag = m.ABMAG(u.Quantity(3.63e-20, u.erg / u.cm ** 2 / u.s / u.Hz))
>>> mag
<ABMAG 0.000233437409711712 mag>
>>> mag.to_flux()
3.630000000000007e-20


See Also
========

Carroll, B. W., & Ostlie, D. A. 1996, An Introduction to Modern Astrophysics (1st ed.; Reading, MA: Addison-Wesley)

`HST ACS photometric systems <http://www.stsci.edu/hst/acs/analysis/zeropoints/#keywords>`_


Reference/API
=============

.. automodapi:: astropy.magnitudes.mags
