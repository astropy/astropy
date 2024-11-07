.. include:: links.inc

.. _spline_models:

****************
1D Spline Models
****************

`~astropy.modeling.spline.Spline1D` models are models which can be used
to fit a piecewise polynomial to a set of data. This means that splines
are closely tied to the method used to fit the spline to the data. Currently,
we provide three methods for fitting splines to data:

- :class:`~astropy.modeling.spline.SplineInterpolateFitter`, which
  fits an interpolating spline to the data. This means that the spline
  will exactly fit all data points.

- :class:`~astropy.modeling.spline.SplineSmoothingFitter`, which fits
  a smoothing spline to the data. This means that the number of knots
  is chosen to satisfy the "smoothing condition":

  .. math:: \sum_{i} \left(w_i * (y_i - spl(x_i))\right)^{2} \leq s

- :class:`~astropy.modeling.spline.SplineExactKnotsFitter`, which fits
  a spline to the data using an exact set of knots. This means that the
  spline will use least-squares regression using the user supplied (interior)
  knots to find the best fit spline to the data.

.. plot::
    :include-source:

    import numpy as np
    import matplotlib.pyplot as plt
    from astropy.modeling.models import Spline1D
    from astropy.modeling.fitting import (SplineInterpolateFitter,
                                          SplineSmoothingFitter,
                                          SplineExactKnotsFitter)

    rng = np.random.default_rng()
    x = np.linspace(-3, 3, 50)
    y = np.exp(-x**2) + 0.1 * rng.standard_normal(50)
    xs = np.linspace(-3, 3, 1000)
    t = [-1, 0, 1]
    spl = Spline1D()

    fitter = SplineInterpolateFitter()
    spl1 = fitter(spl, x, y)

    fitter = SplineSmoothingFitter()
    spl2 = fitter(spl, x, y, s=0.5)

    fitter = SplineExactKnotsFitter()
    spl3 = fitter(spl, x, y, t=t)

    fig, ax = plt.subplots()
    ax.plot(x, y, 'ro', label="Data")
    ax.plot(xs, spl1(xs), 'b-', label="Interpolating")
    ax.plot(xs, spl2(xs), 'g-', label="Smoothing")
    ax.plot(xs, spl3(xs), 'k-', label="Exact Knots")
    ax.legend()
    plt.show()

Note that by default, splines have `degree <astropy.modeling.spline.Spline1D.degree>` 3.
In the case of these splines, the ``degree - 1`` is the number of derivatives that
are matched by the spline across knot points. So for degree 3 splines, the value,
first, and second derivatives of the spline will match across each knot point.

.. warning::

    Splines only support integer degrees, such that ``1 <= degree <= 5``.
