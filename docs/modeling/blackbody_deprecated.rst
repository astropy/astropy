:orphan:

Blackbody (deprecated capabilities)
===================================

Blackbody functions (deprecated)
--------------------------------

Prior to astropy v4.0, two simple functions were provided to compute
a blackbody in frequency and wavelength based units -
:func:`~astropy.modeling.blackbody.blackbody_nu`
and :func:`~astropy.modeling.blackbody.blackbody_lambda`.
These functions are deprecated as they same
results can be obtained using the
:class:`~astropy.modeling.blackbody.BlackBody1D` model which also
provides all the benefits of being a full astropy.modeling  model.

For the :func:`~astropy.modeling.blackbody.blackbody_lambda` and
:func:`~astropy.modeling.blackbody.blackbody_nu` functions, the
flux is calculated with Planck law
(:ref:`Rybicki & Lightman 1979 <ref-rybicki1979>`):

.. math::

    B_{\lambda}(T) = \frac{2 h c^{2} / \lambda^{5}}{exp(h c / \lambda k T) - 1}

    B_{\nu}(T) = \frac{2 h \nu^{3} / c^{2}}{exp(h \nu / k T) - 1}

where the unit of :math:`B_{\lambda}(T)` is
:math:`ergs \; s^{-1} cm^{-2} \mathring{A}^{-1} sr^{-1}`, and
:math:`B_{\nu}(T)` is :math:`ergs \; s^{-1} cm^{-2} Hz^{-1} sr^{-1}`.
:func:`~astropy.modeling.blackbody.blackbody_lambda` and
:func:`~astropy.modeling.blackbody.blackbody_nu` calculate the
blackbody flux for :math:`B_{\lambda}(T)` and :math:`B_{\nu}(T)`,
respectively.

Note that an array of temperatures can also be given instead of a single
temperature. In this case, the Numpy broadcasting rules apply: for instance, if
the frequency and temperature have the same shape, the output will have this
shape too, while if the frequency is a 2-d array with shape ``(n, m)`` and the
temperature is an array with shape ``(m,)``, the output will have a shape
``(n, m)``.

BlackBody1D and blackbody_nu/blackbody_lam
------------------------------------------

The equivalency between the
:class:`~astropy.modeling.blackbody.BlackBody1D`
:func:`~astropy.modeling.blackbody.blackbody_nu` and
:func:`~astropy.modeling.blackbody.blackbody_lambda`
is shown in the plot below.

.. plot::
    :include-source:

    import numpy as np
    import matplotlib.pyplot as plt

    from astropy.modeling.blackbody import blackbody_lambda, blackbody_nu
    from astropy.modeling.blackbody import BlackBody1D
    import astropy.units as u
    from astropy import constants as const

    wavelengths = np.logspace(np.log10(1000), np.log10(3e4), num=1000) * u.AA

    # blackbody parameters
    temperature = 10000 * u.K

    # provides the result in ergs/(cm^2 Hz s sr)
    spectrum_lam = blackbody_lambda(wavelengths, temperature)

    # provides the results in ergs/(cm^2 A s sr)
    spectrum_nu = blackbody_nu(wavelengths, temperature)

    # BlackBody1D provides the results in W/m^2
    #   the amplitude is specified as bolometric flux in ergs(cm^2 s)
    bolometric_flux = const.sigma_sb * temperature ** 4 / np.pi
    bolometric_flux.to(u.erg / (u.cm * u.cm * u.s))
    bb_astro = BlackBody1D(temperature=temperature, bolometric_flux=bolometric_flux)

    # plot the nu and lam results from both methods
    fig, ax = plt.subplots(ncols=2, figsize=(8.0, 4.0))

    # lambda forms
    ax[0].plot(wavelengths, spectrum_lam, label="spectrum_lam", linewidth=6, alpha=0.5)

    # the BlackBody1D has to be converted to the spectrum_nu units
    ax[0].plot(
        wavelengths,
        bb_astro(wavelengths).to(
            u.erg / (u.cm * u.cm * u.Angstrom * u.s),
            equivalencies=u.spectral_density(wavelengths),) / u.sr,
        label="BB1D_lam")

    # nu forms
    ax[1].plot(wavelengths, spectrum_nu, label="spectrum_nu", linewidth=6, alpha=0.5)

    # the BlackBody1D has to be converted to the spectrum_lam units
    ax[1].plot(
        wavelengths,
        bb_astro(wavelengths).to(u.erg / (u.cm * u.cm * u.Hz * u.s)) / u.sr,
        label="BB1D_nu")

    ax[0].set_xlabel(r"$\lambda$ [{}]".format(wavelengths.unit))
    ax[0].set_ylabel(r"$Flux$ [{}]".format(spectrum_lam.unit))
    ax[1].set_xlabel(r"$\lambda$ [{}]".format(wavelengths.unit))
    ax[1].set_ylabel(r"$Flux$ [{}]".format(spectrum_nu.unit))
    ax[0].set_xscale("log")
    ax[0].set_yscale("log")
    ax[1].set_xscale("log")
    ax[1].set_yscale("log")
    ax[0].legend(loc="best")
    ax[1].legend(loc="best")

    plt.tight_layout()
    plt.show()

See Also
^^^^^^^^

.. _ref-rybicki1979:

Rybicki, G. B., & Lightman, A. P. 1979, Radiative Processes in Astrophysics (New York, NY: Wiley)


API for deprecated model/functions
----------------------------------

.. automodapi:: astropy.modeling.blackbody
