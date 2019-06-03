.. _stats-bls:

***********************************
Box least squares (BLS) periodogram
***********************************

The "box least squares (BLS) periodogram" [1]_ is a statistical tool used for
detecting transiting exoplanets and eclipsing binaries in time series
photometric data.
The main interface to this implementation is the `~astropy.timeseries.BoxLeastSquares`
class.


Mathematical Background
=======================

The BLS method finds transit candidates by modeling a transit as a periodic
upside down top hat with four parameters: period, duration, depth, and a
reference time.
In this implementation, the reference time is chosen to be the mid-transit
time of the first transit in the observational baseline.
These parameters are shown in the following sketch:

.. plot::

    import numpy as np
    import matplotlib.pyplot as plt

    period = 6
    t0 = -3
    duration = 2.5
    depth = 0.1
    x = np.linspace(-5, 5, 50000)
    y = np.ones_like(x)
    y[np.abs((x-t0+0.5*period)%period-0.5*period)<0.5*duration] = 1.0 - depth
    plt.figure(figsize=(7, 4))
    plt.axvline(t0, color="k", ls="dashed", lw=0.75)
    plt.axvline(t0+period, color="k", ls="dashed", lw=0.75)
    plt.axhline(1.0-depth, color="k", ls="dashed", lw=0.75)
    plt.plot(x, y)

    kwargs = dict(
        va="center", arrowprops=dict(arrowstyle="->", lw=0.5),
        bbox={"fc": "w", "ec": "none"},
    )
    plt.annotate("period", xy=(t0+period, 1.01), xytext=(t0+0.5*period, 1.01), ha="center", **kwargs)
    plt.annotate("period", xy=(t0, 1.01), xytext=(t0+0.5*period, 1.01), ha="center", **kwargs)
    plt.annotate("duration", xy=(t0-0.5*duration, 1.0-0.5*depth), xytext=(t0, 1.0-0.5*depth), ha="center", **kwargs)
    plt.annotate("duration", xy=(t0+0.5*duration, 1.0-0.5*depth), xytext=(t0, 1.0-0.5*depth), ha="center", **kwargs)
    plt.annotate("reference time", xy=(t0, 1.0-depth-0.01), xytext=(t0+0.25*duration, 1.0-depth-0.01), ha="left", **kwargs)
    plt.annotate("depth", xy=(0.0, 1.0), xytext=(0.0, 1.0-0.5*depth), ha="center", rotation=90, **kwargs)
    plt.annotate("depth", xy=(0.0, 1.0-depth), xytext=(0.0, 1.0-0.5*depth), ha="center", rotation=90, **kwargs)


    plt.ylim(1.0 - depth - 0.02, 1.02)
    plt.xlim(-5, 5)
    plt.gca().set_yticks([])
    plt.gca().set_xticks([])
    plt.ylabel("brightness")
    plt.xlabel("time")

    # ****

Assuming that the uncertainties on the measured flux are known, independent,
and Gaussian, the maximum likelihood in-transit flux can be computed as

.. math::

    y_\mathrm{in} = \frac{\sum_\mathrm{in} y_n/{\sigma_n}^2}{\sum_\mathrm{in} 1/{\sigma_n}^2}

where :math:`y_n` are the brightness measurements, :math:`\sigma_n` are the
associated uncertainties, and both sums are computed over the in-transit data
points.
Similarly, the maximum likelihood out-of-transit flux is

.. math::

    y_\mathrm{out} = \frac{\sum_\mathrm{out} y_n/{\sigma_n}^2}{\sum_\mathrm{out} 1/{\sigma_n}^2}

where these sums are over the out-of-transit observations.
Using these results, the log likelihood of a transit model (maximized over
depth) at a given period :math:`P`, duration :math:`\tau`, and reference time
:math:`t_0` is

.. math::

    \log \mathcal{L}(P,\,\tau,\,t_0) =
    -\frac{1}{2}\,\sum_\mathrm{in}\frac{(y_n-y_\mathrm{in})^2}{{\sigma_n}^2}
    -\frac{1}{2}\,\sum_\mathrm{out}\frac{(y_n-y_\mathrm{out})^2}{{\sigma_n}^2}
    + \mathrm{constant}

This equation might be familiar because it is proportional to the "chi
squared" :math:`\chi^2` for this model and this is a direct consequence of our
assumption of Gaussian uncertainties.
This :math:`\chi^2` is called the "signal residue" by [1]_, so maximizing the
log likelihood over duration and reference time is equivalent to computing the
box least squares spectrum from [1]_.
In practice, this is achieved by finding the maximum likelihood model over a
grid in duration and reference time as specified by the ``durations`` and
``oversample`` parameters for the
`~astropy.timeseries.BoxLeastSquares.power` method.
Behind the scenes, this implementation minimizes the number of required
calculations by pre-binning the observations onto a fine grid following [1]_
and [2]_.


Basic Usage
===========

The transit periodogram takes as input time series observations where the
timestamps ``t`` and the observations ``y`` (usually brightness) are stored as
NumPy arrays or :class:`~astropy.units.Quantity`.
If known, error bars ``dy`` can also optionally be provided.
For example, to evaluate the periodogram for a simulated data set, can be
computed as follows:

>>> import numpy as np
>>> import astropy.units as u
>>> from astropy.timeseries import BoxLeastSquares
>>> np.random.seed(42)
>>> t = np.random.uniform(0, 20, 2000)
>>> y = np.ones_like(t) - 0.1*((t%3)<0.2) + 0.01*np.random.randn(len(t))
>>> model = BoxLeastSquares(t * u.day, y, dy=0.01)
>>> periodogram = model.autopower(0.2)

The output of the `.astropy.timeseries.BoxLeastSquares.autopower` method
is a `~astropy.timeseries.BoxLeastSquaresResults` object with several
useful attributes, the most useful of which are generally the ``period`` and
``power`` attributes.
This result can be plotted using matplotlib:

>>> import matplotlib.pyplot as plt                  # doctest: +SKIP
>>> plt.plot(periodogram.period, periodogram.power)  # doctest: +SKIP

.. plot::

    import numpy as np
    import astropy.units as u
    import matplotlib.pyplot as plt
    from astropy.timeseries import BoxLeastSquares

    np.random.seed(42)
    t = np.random.uniform(0, 20, 2000)
    y = np.ones_like(t) - 0.1*((t%3)<0.2) + 0.01*np.random.randn(len(t))
    model = BoxLeastSquares(t * u.day, y, dy=0.01)
    periodogram = model.autopower(0.2)

    plt.figure(figsize=(8, 4))
    plt.plot(periodogram.period, periodogram.power, "k")
    plt.xlabel("period [day]")
    plt.ylabel("power")

In this figure, you can see the peak at the correct period of 3 days.


Objectives
==========

By default, the `~astropy.timeseries.BoxLeastSquares.power` method computes the log
likelihood of the model fit and maximizes over reference time and duration.
It is also possible to use the signal-to-noise ratio with which the transit
depth is measured as an objective function.
To do this, call `~astropy.timeseries.BoxLeastSquares.power` or
`~astropy.timeseries.BoxLeastSquares.autopower` with ``objective='snr'`` as follows:

>>> model = BoxLeastSquares(t * u.day, y, dy=0.01)
>>> periodogram = model.autopower(0.2, objective="snr")

.. plot::

    import numpy as np
    import astropy.units as u
    import matplotlib.pyplot as plt
    from astropy.timeseries import BoxLeastSquares

    np.random.seed(42)
    t = np.random.uniform(0, 20, 2000)
    y = np.ones_like(t) - 0.1*((t%3)<0.2) + 0.01*np.random.randn(len(t))
    model = BoxLeastSquares(t * u.day, y, dy=0.01)
    periodogram = model.autopower(0.2, objective="snr")

    plt.figure(figsize=(8, 4))
    plt.plot(periodogram.period, periodogram.power, "k")
    plt.xlabel("period [day]")
    plt.ylabel("depth S/N")

This objective will generally produce a periodogram that is qualitatively
similar to the log likelihood spectrum, but it has been used to improve the
reliability of transit search in the presence of correlated noise.


Period Grid
===========

The transit periodogram is always computed on a grid of periods and the
results can be sensitive to the sampling.
As discussed in [1]_, the performance of the transit periodogram method is
more sensitive to the period grid than the
`~astropy.timeseries.LombScargle` periodogram.
This implementation of the transit periodogram includes a conservative
heuristic for estimating the required period grid that is used by the
`~astropy.timeseries.BoxLeastSquares.autoperiod` and
`~astropy.timeseries.BoxLeastSquares.autopower` methods and the details of
this method are given in the API documentation for
`~astropy.timeseries.BoxLeastSquares.autoperiod`.
It is also possible to provide a specific period grid as follows:

>>> model = BoxLeastSquares(t * u.day, y, dy=0.01)
>>> periods = np.linspace(2.5, 3.5, 1000) * u.day
>>> periodogram = model.power(periods, 0.2)

.. plot::

    import numpy as np
    import astropy.units as u
    import matplotlib.pyplot as plt
    from astropy.timeseries import BoxLeastSquares

    np.random.seed(42)
    t = np.random.uniform(0, 20, 2000)
    y = np.ones_like(t) - 0.1*((t%3)<0.2) + 0.01*np.random.randn(len(t))
    model = BoxLeastSquares(t * u.day, y, dy=0.01)
    periods = np.linspace(2.5, 3.5, 1000) * u.day
    periodogram = model.power(periods, 0.2)

    plt.figure(figsize=(8, 4))
    plt.plot(periodogram.period, periodogram.power, "k")
    plt.xlabel("period [day]")
    plt.ylabel("power")

However, if the period grid is too coarse, the correct period can easily be
missed.

>>> model = BoxLeastSquares(t * u.day, y, dy=0.01)
>>> periods = np.linspace(0.5, 10.5, 15) * u.day
>>> periodogram = model.power(periods, 0.2)

.. plot::

    import numpy as np
    import astropy.units as u
    import matplotlib.pyplot as plt
    from astropy.timeseries import BoxLeastSquares

    np.random.seed(42)
    t = np.random.uniform(0, 20, 2000)
    y = np.ones_like(t) - 0.1*((t%3)<0.2) + 0.01*np.random.randn(len(t))
    model = BoxLeastSquares(t * u.day, y, dy=0.01)
    periods = np.linspace(0.5, 10.5, 15) * u.day
    periodogram = model.power(periods, 0.2)

    plt.figure(figsize=(8, 4))
    plt.plot(periodogram.period, periodogram.power, "k")
    plt.xlabel("period [day]")
    plt.ylabel("power")


Peak Statistics
===============

To help in the transit vetting process and to debug problems with candidate
peaks, the `~astropy.timeseries.BoxLeastSquares.compute_stats` method can be
used to calculate several statistics of a candidate transit.
Many of these statistics are based on the VARTOOLS package described in [2]_.
This will often be used as follows to compute stats for the maximum point in
the periodogram:

>>> model = BoxLeastSquares(t * u.day, y, dy=0.01)
>>> periodogram = model.autopower(0.2)
>>> max_power = np.argmax(periodogram.power)
>>> stats = model.compute_stats(periodogram.period[max_power],
...                             periodogram.duration[max_power],
...                             periodogram.transit_time[max_power])

This calculates a dictionary with statistics about this candidate.
Each entry in this dictionary is described in the documentation for
`~astropy.timeseries.BoxLeastSquares.compute_stats`.


Literature References
=====================

.. [1] Kovacs, Zucker, & Mazeh (2002), A&A, 391, 369 (arXiv:astro-ph/0206099)
.. [2] Hartman & Bakos (2016), Astronomy & Computing, 17, 1 (arXiv:1605.06811)
