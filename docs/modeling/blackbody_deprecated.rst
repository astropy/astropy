:orphan:

.. _deprecated-blackbody:

Blackbody Module (deprecated capabilities)
==========================================

BlackBody1D class
-----------------

The equivalency between the
:class:`~astropy.modeling.physical_models.BlackBody` and
:class:`~astropy.modeling.blackbody.BlackBody1D`
is shown in the plot below.

.. plot::
    :nofigs:
    :context: reset

    # Ignore deprecation warnings from Blackbody1D
    import warnings
    warnings.simplefilter('ignore')


.. plot::
    :include-source:
    :context:

    import numpy as np
    import matplotlib.pyplot as plt

    from astropy.modeling.blackbody import BlackBody1D
    from astropy.modeling.models import BlackBody
    import astropy.units as u

    wavelengths = np.logspace(np.log10(1000), np.log10(3e4), num=1000) * u.AA

    # blackbody parameters
    temperature = 10000 * u.K

    # BlackBody1D is normalized by default to have a bolometric flux of 1
    bb1d = BlackBody1D(temperature)

    # default results from BlackBody same units as BlackBody1D
    bb = BlackBody(temperature)

    # plot the results from both methods
    fig, ax = plt.subplots(ncols=1, figsize=(4.0, 4.0))

    ax.plot(wavelengths, bb1d(wavelengths) / u.sr, label="BlackBody1D", linewidth=6, alpha=0.5)

    # divide by bolometric flux to get equivalent of BlackBody1D
    ax.plot(wavelengths, bb(wavelengths)/bb.bolometric_flux.value, label="BlackBody")


    ax.set_xlabel(r"$\lambda$ [{}]".format(wavelengths.unit))
    ax.set_ylabel(r"$Flux$ [{}]".format(bb(wavelengths).unit))
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend(loc="best")

    plt.tight_layout()
    plt.show()

blackbody functions (deprecated)
--------------------------------

Prior to astropy v4.0, two simple functions were provided to compute
a blackbody in frequency and wavelength based units -
:func:`~astropy.modeling.blackbody.blackbody_nu`
and :func:`~astropy.modeling.blackbody.blackbody_lambda`.
These functions are deprecated as they same
results can be obtained using the
:class:`~astropy.modeling.physical_models.BlackBody` model.

Equivalent results between
and :class:`~astropy.modeling.physical_models.BlackBody` and the
func:`~astropy.modeling.blackbody.blackbody_nu` and
func:`~astropy.modeling.blackbody.blackbody_lambda` is shown in the plot below.

.. plot::
    :nofigs:
    :context: reset

    # Ignore deprecation warnings from Blackbody1D
    import warnings
    warnings.simplefilter('ignore')


.. plot::
    :include-source:
    :context:

    import numpy as np
    import matplotlib.pyplot as plt

    from astropy.modeling.blackbody import blackbody_lambda, blackbody_nu
    from astropy.modeling.models import BlackBody
    import astropy.units as u

    wavelengths = np.logspace(np.log10(1000), np.log10(3e4), num=1000) * u.AA

    # blackbody parameters
    temperature = 10000 * u.K

    # provides the results in ergs/(cm^2 Hz s sr)
    spectrum_nu = blackbody_nu(wavelengths, temperature)

    # provides the result in ergs/(cm^2 A s sr)
    spectrum_lam = blackbody_lambda(wavelengths, temperature)

    # default results from BlackBody same units as blackbody_nu
    bb_nu = BlackBody(temperature)

    # use scale to change results units for Blackbody to match blackbody_lambda
    bb_lam = BlackBody(temperature, scale=1.0 * u.erg / (u.cm ** 2 * u.AA * u.s * u.sr))

    # plot the nu and lam results from both methods
    fig, ax = plt.subplots(ncols=2, figsize=(8.0, 4.0))

    # lambda forms
    ax[0].plot(wavelengths, spectrum_lam, label="blackbody_lambda", linewidth=6, alpha=0.5)
    ax[0].plot(wavelengths, bb_lam(wavelengths), label="BlackBody")

    # nu forms
    ax[1].plot(wavelengths, spectrum_nu, label="blackbody_nu", linewidth=6, alpha=0.5)
    ax[1].plot(wavelengths, bb_nu(wavelengths), label="BlackBody")

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

API for deprecated model/functions
----------------------------------

.. automodapi:: astropy.modeling.blackbody
