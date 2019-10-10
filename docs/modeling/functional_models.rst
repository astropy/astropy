.. _functional_models:

Functional Models
*****************

These are models that are mathematically motivated, generally as solutions to
mathematical problems.   See :doc:`physical_models` for physically motivated
models.

.. Lorentz1D, Sersic1D, Sersic2D, Voigt1D potential for moving to physical models

1D Models
---------

- :class:`~astropy.modeling.functional_models.Box1D` model computes a box
  function with an amplitude centered at x_0 with the specified width.

- :class:`~astropy.modeling.functional_models.Const1D` model returns the
  constant replicated by the number of input x values.

- :class:`~astropy.modeling.functional_models.Gaussian1D`

- :class:`~astropy.modeling.functional_models.Linear1D` model provides a
  line parameterizied by the slope and y-intercept

- :class:`~astropy.modeling.functional_models.MexicanHat1D`

- :class:`~astropy.modeling.functional_models.Moffat1D`

- :class:`~astropy.modeling.functional_models.Lorentz1D`

- :class:`~astropy.modeling.functional_models.Multiply` model multiples the
  input x values by a factor and propagates units if the factor is
  a :class:`~astropy.units.Quantity`.

- :class:`~astropy.modeling.functional_models.RedshiftScaleFactor` model
  multiples the input x values by a (1 + z) factor.

- :class:`~astropy.modeling.functional_models.Scale` model multiples by a
  factor without changing the units of the result.

- :class:`~astropy.modeling.functional_models.Sersic1D`

- :class:`~astropy.modeling.functional_models.Shift` model adds a constant
  to the input x values.

- :class:`~astropy.modeling.functional_models.Sine1D`

- :class:`~astropy.modeling.functional_models.Trapezoid1D`

- :class:`~astropy.modeling.functional_models.Voigt1D`

- :class:`~astropy.modeling.functional_models.KingProjectedAnalytic1D`

Plot showing some of the 1D models

.. plot::
    :include-source:

    import numpy as np
    import matplotlib.pyplot as plt

    from astropy.modeling.models import (Box1D, Gaussian1D, MexicanHat1D,
                                         Moffat1D, Lorentz1D, Voigt1D)
    import astropy.units as u

    x = np.linspace(-5.0, 5.0, num=100)

    boxmod = Box1D(amplitude=10., x_0=0.0, width=1.0)
    gaussmod = Gaussian1D(amplitude=10., mean=0.0, stddev=1.0)
    mexhatmod = MexicanHat1D(amplitude=10., x_0=0.0, sigma=1.0)
    moffatmod = Moffat1D(amplitude=10., x_0=0.0, gamma=1.0)
    lorentzmod = Lorentz1D(amplitude=10., x_0=0.0, fwhm=1.0)
    voigtmod = Voigt1D(amplitude_L=10., x_0=0.0, fwhm_L=1.0, fwhm_G=1.0)

    fig, ax = plt.subplots()
    ax.plot(x, boxmod(x), label="Box1D")
    ax.plot(x, gaussmod(x), label="Gaussian1D")
    ax.plot(x, mexhatmod(x), label="MexicanHat1D")
    ax.plot(x, moffatmod(x), label="Moffat1D")
    ax.plot(x, lorentzmod(x), label="Lorentz1D")
    ax.plot(x, voigtmod(x), label="Voigt1D")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.legend()
    plt.tight_layout()
    plt.show()

2D Models
---------

- :class:`~astropy.modeling.functional_models.AiryDisk2D`

- :class:`~astropy.modeling.functional_models.Box2D`

- :class:`~astropy.modeling.functional_models.Const2D` model returns the
  constant replicated by the number of input x and y values.

- :class:`~astropy.modeling.functional_models.Disk2D`

- :class:`~astropy.modeling.functional_models.Ellipse2D`

- :class:`~astropy.modeling.functional_models.Gaussian2D`

- :class:`~astropy.modeling.functional_models.MexicanHat2D`

- :class:`~astropy.modeling.functional_models.Planar2D`

- :class:`~astropy.modeling.functional_models.Sersic2D`

- :class:`~astropy.modeling.functional_models.TrapezoidDisk2D`

- :class:`~astropy.modeling.functional_models.Ring2D`
