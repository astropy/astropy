.. _stats-ripley:

******************************
Ripley's K Function Estimators
******************************

Ripley's K function is used to characterize spatial point processes.
More precisely, it describes correlation in point fields.
The :class: `~astropy.stats.RipleysKEstimator` class implements some
estimators for this function which take into account edge-effects.

Basic Usage
===========

References
==========
.. [1] Ripley, B.D. *The second-order analysis of stationary point processes*.
       Journal of Applied Probability. 13: 255–266, 1976.
.. [1] Spatial descriptive statistics.
       <https://en.wikipedia.org/wiki/Spatial_descriptive_statistics>
.. [2] Cressie, N.A.C. *Statistics for Spatial Data*, Wiley, New York.
.. [3] Stoyan, D., Stoyan, H. *Fractals, Random Shapes and Point Fields*,
       Akademie Verlag GmbH, Chichester, 1992.
