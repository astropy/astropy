Physical Models
***************

These are models that are physics motivated (as opposed to mathematically motivated).

BlackBody
=========

The :class:`~astropy.modeling.physical_models.BlackBody` model provides a model
for using `Plank's Law <https://en.wikipedia.org/wiki/Planck%27s_law/>`_.
The blackbody function is

.. math::

   B_{\nu}(T) = A \frac{2 h \nu^{3} / c^{2}}{exp(h \nu / k T) - 1}

where :math:`\nu` is the frequency, :math:`T` is the temperature,
:math:`A` is the scaling factor,
:math:`h` is the Plank constant, :math:`c` is the speed of light, and
:math:`k` is the Boltzmann constant.

The two parameters of the model the unitless scaling factor ``scale`` and
the absolute temperature ``temperature``.  The ``scale`` factor is unitless
as the Plank function has natural units of spectral radiance, specifically
ergs/(cm^2 Hz s sr).  The
temperature can be passed as a Quantity with any temperature unit.  If temperature
is passed without a unit, it is assumed to be in Kelvin units.

An example plot for a blackbody with a temperature of 10000 K and a scale of 1 is
shown below.

.. plot::
    :include-source:

    import numpy as np
    import matplotlib.pyplot as plt

    from astropy.modeling.models import BlackBody
    import astropy.units as u

    wavelengths = np.logspace(np.log10(1000), np.log10(3e4), num=1000) * u.AA

    # blackbody parameters
    temperature = 10000 * u.K

    # BlackBody1D provides the results in ergs/(cm^2 Hz s sr)
    bb = BlackBody(temperature=temperature, scale=10000.0)
    bb_result = bb(wavelengths)

    fig, ax = plt.subplots(ncols=1)
    ax.plot(wavelengths, bb_result, '-')

    ax.set_xscale('log')
    ax.set_xlabel(r"$\lambda$ [{}]".format(wavelengths.unit))
    ax.set_ylabel(r"$Flux$ [{}]".format(bb_result.unit))

    plt.tight_layout()
    plt.show()

The Planck function can be needed as :math:`B_\lambda` in units of
ergs/(cm^2 s A sr).  This can be easily accomplished with a
unit conversion as illustrated below.  Note that this uses the ``spectral_density``
equivalencies requiring the steradian units to be first removed, then reapplied
after the conversion.

.. code-block:: python

    import numpy as np
    import astropy.units as u
    from astropy.modeling.models import BlackBody

    wavelengths = np.logspace(np.log10(1000), np.log10(3e4), num=1000) * u.AA
    bb = BlackBody(temperature=10000.0 * u.K, scale=10000.0)
    bb_nu = bb(wavelengths)
    bb_lambda = (bb_nu * u.sr).to(
        u.erg / (u.cm * u.cm * u.Angstrom * u.s),
        equivalencies=u.spectral_density(wavelengths),
    ) / u.sr

The ``bolometric_flux()`` member function gives the bolometric flux using
:math:`\sigma T^4/\pi` where :math:`\sigma` is the Stefan-Boltzmann constant.

.. note::

    Prior to v4.0, the BlackBody1D and the functions blackbody_nu and blackbody_lambda
    were provided.  BlackBody1D was a more limited blackbody model that was
    specific to integrated fluxes from sources.  See :doc:`blackbody_deprecated` for details.
    The capabilities of all three
    can be obtained with :class:`~astropy.modeling.physical_models.BlackBody`.
