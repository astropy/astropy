.. doctest-skip-all

.. _astropy_analytic_functions:

***************************************************
Analytic Functions (``astropy.analytic_functions``)
***************************************************

Introduction
============

The ``astropy.analytic_functions`` subpackage provides analytic functions that
are commonly used in astronomy. Unlike similar functions from other software,
these ones already understand `~astropy.units.Quantity`, i.e., they can handle
units of input and output parameters.


Getting Started
===============

>>> from astropy import units as u
>>> from astropy import analytic_functions

Calculate blackbody flux for 10000 K at 6000 Angstrom:

>>> analytic_functions.planck(6000 * u.AA, 10000 * u.K)
<Quantity 0.00018391673686797086 erg / (cm2 Hz s sr)>


Using ``astropy.analytic_functions``
====================================

.. _planck-law:

Blackbody Radiation
-------------------

Blackbody flux is calculated with Planck law
(:ref:`Rybicki & Lightman 1979 <ref-rybicki1979>`).

.. math::

    B_{\lambda}(T) = \frac{2 h c^{2} / \lambda^{5}}{exp(h c / \lambda k T) - 1}

where the unit of :math:`B_{\lambda}(T)` is
:math:`erg \; s^{-1} cm^{-2} \AA^{-1} sr^{-1}`.

:func:`~astropy.analytic_functions.planck` calculates the blackbody
flux for the given wavelength/frequency/wave number, temperature, and desired
flux unit. By default, it returns flux for :math:`B_{\nu}(T)`. If overflow or
underflow occurs, the affected flux values are set to zeroes.

Examples
^^^^^^^^

Calculate blackbody flux in :math:`B_{\lambda}(T)` for 5000 K at 100 and 10000
Angstrom:

>>> from astropy import units as u
>>> from astropy.analytic_functions import planck
>>> flam = u.erg / (u.cm ** 2 * u.s * u.AA)
>>> planck([100, 10000], 5000, flux_unit=flam)
<Quantity [0., 710190.52597893] erg / (Angstrom cm2 s sr)>


See Also
========

.. _ref-rybicki1979:

Rybicki, G. B., & Lightman, A. P. 1979, Radiative Processes in Astrophysics (New York, NY: Wiley)


Reference/API
=============

.. automodapi:: astropy.analytic_functions
