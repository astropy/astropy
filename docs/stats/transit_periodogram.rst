.. _stats-transit_periodogram:

********************
Transit Periodograms
********************

The "transit periodogram" (also known the "box least squares spectrum"
following [1]_) is a statistical tool used for detecting transiting exoplanets
and eclipsing binaries in time series photometric data.
The main interface is through the :class:`~astropy.stats.TransitPeriodogram`
class.


Mathematical Background
=======================

The transit periodogram method finds transit candidates by modeling a transit
as a periodic upside down top hat with four parameters: period, duration,
depth, and a reference time.
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

Literature References
=====================
.. [1] Kovacs, Zucker, & Mazeh (2002), A&A, 391, 369 (arXiv:astro-ph/0206099)
