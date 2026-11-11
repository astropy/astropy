.. _stats-pdm:

Phase Dispersion Minimization
*****************************

The Phase Dispersion Minimization (PDM) is a non-parametric method described by Stellingwerf (1978) [1]_ used
for period determination of time series data like light curves. This method is based on the idea of folding
the data over a range of trial periods and evaluating the dispersion of the data in some bins. The period
that minimizes the dispersion is considered the best estimate of the true period.

The :class:`~astropy.timeseries.periodograms.PhaseDispersionMinimization` class implements the Phase Dispersion Minimization method.

Basic Usage
===========
PDM is used to detect periods in unevenly spaced and supports data with gaps.

.. Important::
    The PDM method is not suitable for data with a large number of outliers, it needs enough observations at
    similar phases to estimate meaningful within-bin variances, and it is not suitable for data that shows strong
    trends because the dispersion between bins will be dominated by the trend rather than the periodic signal.

Example
-------

.. EXAMPLE START: Using the Phase Dispersion Minimization class to compute a periodogram for a noisy sine wave

>>> import numpy as np
>>> rand = np.random.default_rng(67)
>>> t = 100 * rand.random(100)
>>> y = np.sin(2 * np.pi * t) + 0.01 * np.random.randn(len(t))

100 noisy observations of a sine wave with a period of 1.0

Now we can use the :class:`~astropy.timeseries.periodograms.PhaseDispersionMinimization` class to compute the PDM periodogram.

>>> from astropy.timeseries.periodograms import PhaseDispersionMinimization
>>> pdm = PhaseDispersionMinimization(t, y)
>>> periods = np.linspace(0.1, 10, 1000)
>>> thetas = pdm.run(periods=periods)

We can plot the resulting periodogram to visualize the results.

>>> import matplotlib.pyplot as plt
>>> plt.plot(periods, thetas)
>>> plt.xlabel('Period')
>>> plt.ylabel('Theta')
>>> plt.show()

.. plot::

    import numpy as np
    import matplotlib.pyplot as plt
    from astropy.timeseries.periodograms import PhaseDispersionMinimization

    rand = np.random.default_rng(67)
    t = 100 * rand.random(100)
    y = np.sin(2 * np.pi * t) + 0.01 * np.random.randn(len(t))

    pdm = PhaseDispersionMinimization(t, y)
    periods = np.linspace(0.1, 10, 1000)
    thetas = pdm.run(periods=periods)

    plt.plot(periods, thetas)
    plt.xlabel('Period')
    plt.ylabel('Theta')
    plt.show()


.. EXAMPLE END

References
==========
.. [1] Stellingwerf, R. F., “Period determination using phase dispersion minimization.”,
    <i>The Astrophysical Journal</i>, vol. 224, IOP, pp. 953–960, 1978. doi:10.1086/156444.
