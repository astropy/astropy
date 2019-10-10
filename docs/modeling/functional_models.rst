.. _functional_models:

Functional Models
*****************

These are models that are mathematically motivated, generally as solutions to
mathematical problems.   See :doc:`physical_models` for physically motivated
models.

.. Lorentz1D, Sersic1D, Sersic2D, Voigt1D potential for moving to physical models

1D Models
---------

- :class:`~astropy.modeling.functional_models.Box1D`

- :class:`~astropy.modeling.functional_models.Const1D` model returns the
  constant replicated by the number of input x values.

- :class:`~astropy.modeling.functional_models.Gaussian1D`

- :class:`~astropy.modeling.functional_models.Linear1D` model provides a
  line parameterizied by the slope and y-intercept

- :class:`~astropy.modeling.functional_models.MexicanHat1D`

- :class:`~astropy.modeling.functional_models.Moffat1D`

- :class:`~astropy.modeling.functional_models.Lorentz1D`

- :class:`~astropy.modeling.functional_models.Multiply` model multiples by a
  factor and propagates units if the factor is a :class:`~astropy.units.Quantity`.

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

- :class:`~astropy.modeling.functional_models.TrapazoidDisk2D`

- :class:`~astropy.modeling.functional_models.Ring2D`
