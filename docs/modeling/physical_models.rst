Physical Models
***************

These are models that are physical motivated, generally as solutions to
physical problems.  This is in contrast to those that are mathematically motivated,
generally as solutions to mathematical problems.

BlackBody
=========

.. _blackbody-planck-law:

The :class:`~astropy.modeling.physical_models.BlackBody` model provides a model
for using `Planck's Law <https://en.wikipedia.org/wiki/Planck%27s_law/>`_.
The blackbody function is

.. math::

   B_{\nu}(T) = A \frac{2 h \nu^{3} / c^{2}}{exp(h \nu / k T) - 1}

where :math:`\nu` is the frequency, :math:`T` is the temperature,
:math:`A` is the scaling factor,
:math:`h` is the Plank constant, :math:`c` is the speed of light, and
:math:`k` is the Boltzmann constant.

The two parameters of the model the scaling factor ``scale`` (A) and
the absolute temperature ``temperature`` (T).  If the ``scale`` factor does not
have units, then the result is in units of spectral radiance, specifically
ergs/(cm^2 Hz s sr).  If the ``scale`` factor is passed with spectral radiance units,
then the result is in those units (e.g., ergs/(cm^2 A s sr) or MJy/sr).
Setting the ``scale`` factor with units of ergs/(cm^2 A s sr) will give the
Planck function as :math:`B_\lambda`.
The temperature can be passed as a Quantity with any supported temperature unit.
If temperature
is passed without a unit, it is assumed to be in Kelvin.

An example plot for a blackbody with a temperature of 10000 K and a scale of 1 is
shown below.  A scale of 1 shows the Planck function with no scaling in the
default units returned by :class:`~astropy.modeling.physical_models.BlackBody`.

.. plot::
    :include-source:

    import numpy as np
    import matplotlib.pyplot as plt

    from astropy.modeling.models import BlackBody
    import astropy.units as u

    wavelengths = np.logspace(np.log10(1000), np.log10(3e4), num=1000) * u.AA

    # blackbody parameters
    temperature = 10000 * u.K

    # BlackBody1D provides the results in ergs/(cm^2 Hz s sr) when scale has no units
    bb = BlackBody(temperature=temperature, scale=10000.0)
    bb_result = bb(wavelengths)

    fig, ax = plt.subplots(ncols=1)
    ax.plot(wavelengths, bb_result, '-')

    ax.set_xscale('log')
    ax.set_xlabel(r"$\lambda$ [{}]".format(wavelengths.unit))
    ax.set_ylabel(r"$F(\lambda)$ [{}]".format(bb_result.unit))

    plt.tight_layout()
    plt.show()

The ``bolometric_flux`` member function gives the bolometric flux using
:math:`\sigma T^4/\pi` where :math:`\sigma` is the Stefan-Boltzmann constant.

The ``lambda_max`` and ``nu_max`` member functions give the wavelength and
frequency of the maximum for :math:`B_\lambda` and :math:`B_\nu`,
respectively, calculated using `Wein's Law
<https://en.wikipedia.org/wiki/Wien%27s_displacement_law>`_.

.. note::

    Prior to v4.0, the ``BlackBody1D`` and the functions ``blackbody_nu`` and ``blackbody_lambda``
    were provided.  ``BlackBody1D`` was a more limited blackbody model that was
    specific to integrated fluxes from sources.  See :doc:`blackbody_deprecated`
    and `astropy issue #9066 <https://github.com/astropy/astropy/issues/9066>`_ for details.
    The capabilities of all three
    can be obtained with :class:`~astropy.modeling.physical_models.BlackBody`.
