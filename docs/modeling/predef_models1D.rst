.. _predef_models1D:

*********
1D Models
*********

Operations
==========

These models perform simple mathematical operations.

- :class:`~astropy.modeling.functional_models.Const1D` model returns the
  constant replicated by the number of input x values.

- :class:`~astropy.modeling.functional_models.Multiply` model multiples the
  input x values by a factor and propagates units if the factor is
  a :class:`~astropy.units.Quantity`.

- :class:`~astropy.modeling.functional_models.RedshiftScaleFactor` model
  multiples the input x values by a (1 + z) factor.

- :class:`~astropy.modeling.functional_models.Scale` model multiples by a
  factor without changing the units of the result.

- :class:`~astropy.modeling.functional_models.Shift` model adds a constant
  to the input x values.

Shapes
======

These models provide shapes, often used to model general x, y data.

- :class:`~astropy.modeling.functional_models.Linear1D` model provides a
  line parameterizied by the slope and y-intercept

- :class:`~astropy.modeling.functional_models.Sine1D` model provides a sine
  parameterized by an amplitude, frequency, and phase.

.. plot::

    import numpy as np
    import matplotlib.pyplot as plt

    from astropy.modeling.models import (Linear1D, Sine1D)

    x = np.linspace(-4.0, 6.0, num=100)

    fig, sax = plt.subplots(ncols=2, figsize=(10, 5))
    ax = sax.flatten()

    linemod = Linear1D(slope=2., intercept=1.)
    ax[0].plot(x, linemod(x), label="Linear1D")

    sinemod = Sine1D(amplitude=10., frequency=0.5, phase=0.)
    ax[1].plot(x, sinemod(x), label="Sine1D")
    ax[1].set_ylim(-11.0, 13.0)

    for k in range(2):
        ax[k].set_xlabel("x")
        ax[k].set_ylabel("y")
        ax[k].legend()

    plt.tight_layout()
    plt.show()

Profiles
========

These models provide profiles, often used for lines in spectra.

- :class:`~astropy.modeling.functional_models.Box1D` model computes a box
  function with an amplitude centered at x_0 with the specified width.

- :class:`~astropy.modeling.functional_models.Gaussian1D` model computes
  a Gaussian with an amplitude centered at x_0 with the specified width.

- :class:`~astropy.modeling.functional_models.KingProjectedAnalytic1D` model
  computes the analytic form of the a King model with an amplitude and
  core and tidal radii.

- :class:`~astropy.modeling.functional_models.Lorentz1D` model computes
  a Lorentzian with an amplitude centered at x_0 with the specified width.

- :class:`~astropy.modeling.functional_models.RickerWavelet1D` model computes
  a RickerWavelet function with an amplitude centered at x_0 with the specified width.

- :class:`~astropy.modeling.functional_models.Moffat1D` model computes a
  Moffat function with an amplitude centered at x_0 with the specified width.

- :class:`~astropy.modeling.functional_models.Sersic1D` model
  computes a Sersic model with an amplitude with an effective radius and
  the specified sersic index.

- :class:`~astropy.modeling.functional_models.Trapezoid1D` model computes a
  box with sloping sides with an amplitude centered at x_0 with the specified
  width and sides wit the specified slope.

- :class:`~astropy.modeling.functional_models.Voigt1D` model computes a
  Voigt function with an amplitude centered at x_0 with the specified
  Lorentzian and Gaussian widths.

.. plot::

    import numpy as np
    import matplotlib.pyplot as plt

    from astropy.modeling.models import (
        Box1D,
        Gaussian1D,
        RickerWavelet1D,
        Moffat1D,
        Lorentz1D,
        Sersic1D,
        Trapezoid1D,
        KingProjectedAnalytic1D,
        Voigt1D,
    )

    x = np.linspace(-4.0, 6.0, num=100)
    r = np.logspace(-1.0, 2.0, num=100)

    fig, sax = plt.subplots(nrows=3, ncols=3, figsize=(10, 10))
    ax = sax.flatten()

    mods = [
        Box1D(amplitude=10.0, x_0=1.0, width=1.0),
        Gaussian1D(amplitude=10.0, mean=1.0, stddev=1.0),
        KingProjectedAnalytic1D(amplitude=10.0, r_core=1.0, r_tide=10.0),
        Lorentz1D(amplitude=10.0, x_0=1.0, fwhm=1.0),
        RickerWavelet1D(amplitude=10.0, x_0=1.0, sigma=1.0),
        Moffat1D(amplitude=10.0, x_0=1.0, gamma=1.0, alpha=1.),
        Sersic1D(amplitude=10.0, r_eff=1.0 / 2.0, n=5),
        Trapezoid1D(amplitude=10.0, x_0=1.0, width=1.0, slope=5.0),
        Voigt1D(amplitude_L=10.0, x_0=1.0, fwhm_L=1.0, fwhm_G=1.0),
    ]

    for k, mod in enumerate(mods):
        cname = mod.__class__.__name__
        ax[k].set_title(cname)
        if cname in ["KingProjectedAnalytic1D", "Sersic1D"]:
            ax[k].plot(r, mod(r))
            ax[k].set_xscale("log")
            ax[k].set_yscale("log")
        else:
            ax[k].plot(x, mod(x))

    for k in range(len(mods)):
        ax[k].set_xlabel("x")
        ax[k].set_ylabel("y")

    # remove axis for any plots not used
    for k in range(len(mods), len(ax)):
        ax[k].axis("off")

    plt.tight_layout()
    plt.show()
