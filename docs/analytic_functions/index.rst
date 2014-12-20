.. doctest-skip-all

.. _astropy_analytic_functions:

***************************************************
Analytic Functions (``astropy.analytic_functions``)
***************************************************

Introduction
============

The ``astropy.analytic_functions`` subpackage provides analytic functions that
are commonly used in astronomy. These already understand
`~astropy.units.Quantity`, i.e., they can handle units of input and output
parameters.

In future versions of ``astropy``, many of these might be accessible as
`~astropy.modeling.core.Model`.


Getting Started
===============

>>> from astropy import units as u
>>> from astropy.analytic_functions import blackbody_lambda, blackbody_nu

Calculate blackbody flux for 10000 K at 6000 Angstrom:

>>> blackbody_lambda(6000 * u.AA, 10000 * u.K)
<Quantity 15315791.836941158 erg / (Angstrom cm2 s sr)>
>>> blackbody_nu(6000 * u.AA, 10000 * u.K)
<Quantity 0.00018391673686797075 erg / (cm2 Hz s sr)


Using ``astropy.analytic_functions``
====================================

.. _blackbody-planck-law:

Blackbody Radiation
-------------------

Blackbody flux is calculated with Planck law
(:ref:`Rybicki & Lightman 1979 <ref-rybicki1979>`):

.. math::

    B_{\lambda}(T) = \frac{2 h c^{2} / \lambda^{5}}{exp(h c / \lambda k T) - 1}

    B_{\nu}(T) = \frac{2 h \nu^{3} / c^{2}}{exp(h \nu / k T) - 1}

where the unit of :math:`B_{\lambda}(T)` is
:math:`erg \; s^{-1} cm^{-2} \AA^{-1} sr^{-1}`, and :math:`B_{\nu}(T)` is
:math:`erg \; s^{-1} cm^{-2} Hz^{-1} sr^{-1}`.
:func:`~astropy.analytic_functions.blackbody.blackbody_lambda` and
:func:`~astropy.analytic_functions.blackbody.blackbody_nu` calculate the
blackbody flux for :math:`B_{\lambda}(T)` and :math:`B_{\nu}(T)`, respectively.

.. _blackbody-examples:

Examples
^^^^^^^^

>>> import numpy as np
>>> from astropy import units as u
>>> from astropy.analytic_functions import blackbody_lambda, blackbody_nu

Calculate blackbody flux for 5000 K at 100 and 10000 Angstrom while suppressing
any Numpy warnings:

>>> wavelengths = [100, 10000] * u.AA
>>> temperature = 5000 * u.K
>>> with np.errstate(all='ignore'):
...     flux_lam = blackbody_lambda(wavelengths, temperature)
...     flux_nu = blackbody_nu(wavelengths, temperature)
>>> flux_lam
Quantity [1.27452545e-108, 7.10190526e+005] erg / (Angstrom cm2 s sr)>
>>> flux_nu
<Quantity [4.25135927e-123, 2.36894060e-005] erg / (cm2 Hz s sr)>

Plot a blackbody spectrum for 5000 K:

.. plot::

    import matplotlib.pyplot as plt
    import numpy as np
    from astropy import constants as const
    from astropy import units as u
    from astropy.analytic_functions import blackbody_lambda

    temperature = 5000 * u.K
    wavemax = (const.b_wien / temperature).to(u.AA)  # Wien's displacement law
    waveset = np.logspace(
        0, np.log10(wavemax.value + 10 * wavemax.value), num=1000) * u.AA
    with np.errstate(all='ignore'):
        flux = blackbody_lambda(waveset, temperature)

    fig, ax = plt.subplots(figsize=(8,5))
    ax.plot(waveset.value, flux.value)
    ax.axvline(wavemax.value, ls='--')
    ax.get_yaxis().get_major_formatter().set_powerlimits((0, 1))
    ax.set_xlabel(r'$\lambda$ ({0})'.format(waveset.unit))
    ax.set_ylabel(r'$B_{\lambda}(T)$')
    ax.set_title('Blackbody, T = {0}'.format(temperature))


See Also
========

.. _ref-rybicki1979:

Rybicki, G. B., & Lightman, A. P. 1979, Radiative Processes in Astrophysics (New York, NY: Wiley)


Reference/API
=============

.. automodapi:: astropy.analytic_functions.blackbody
   :no-inheritance-diagram:
