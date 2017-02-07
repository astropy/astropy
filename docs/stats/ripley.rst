.. _stats-ripley:

******************************
Ripley's K Function Estimators
******************************

Ripley's K function is used to characterize spatial point processes.
More precisely, it describes correlation in point fields.
The :class: `~astropy.stats.RipleysKEstimator` class implements some
estimators for this function which provides several methods for
edge-effects correction.

Basic Usage
===========

The actual implementation of Ripley's K function estimators lie in the method
``evaluate`` which take the following arguments ``data``, ``radii``, and,
optionally ``mode``.

The ``data`` argument is a 2D array which represents the set of observed
points (events) in the area of study. The ``radii`` argument corresponds to a
set of distances for which the estimator will be evaluated. The ``mode``
argument takes a value on the following linguistic set
``{none, translation, ohser, var-width, ripley}``; each keyword represents a
different method to perform correction due to edge-effects. See the API
documentation and references for details about these methods.

Instances of :class: `~astropy.stats.RipleysKEstimator` can also be used as
callables (which is equivalent to calling the ``evaluate`` method).

A minimal usage example is shown as follows:

.. plot::
    :include-source:
    import numpy as np
    from matplotlib import pyplot as plt
    from astropy.stats import RipleysKEstimator

    z = np.random.uniform(low=5, high=10, size=(100, 2))
    Kest = RipleysKEstimator(area=25, x_max=10, y_max=10, x_min=5, y_min=5)

    r = np.linspace(0, 2.5, 100)
    plt.plot(r, Kest.poisson(r))
    plt.plot(r, Kest(data=z, radii=r, mode='none'))
    plt.plot(r, Kest(data=z, radii=r, mode='translation'))
    plt.plot(r, Kest(data=z, radii=r, mode='ohser'))
    plt.plot(r, Kest(data=z, radii=r, mode='var-width'))
    plt.plot(r, Kest(data=z, radii=r, mode='ripley'))

References
==========
.. [1] Ripley, B.D. *The second-order analysis of stationary point processes*.
       Journal of Applied Probability. 13: 255–266, 1976.
.. [1] *Spatial descriptive statistics*.
       <https://en.wikipedia.org/wiki/Spatial_descriptive_statistics>
.. [2] Cressie, N.A.C. *Statistics for Spatial Data*, Wiley, New York.
.. [3] Stoyan, D., Stoyan, H. *Fractals, Random Shapes and Point Fields*,
       Akademie Verlag GmbH, Chichester, 1992.
