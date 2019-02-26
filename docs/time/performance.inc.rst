.. note that if this is changed from the default approach of using an *include*
   (in index.rst) to a separate performance page, the header needs to be changed
   from === to ***, the filename extension needs to be changed from .inc.rst to
   .rst, and a link needs to be added in the subpackage toctree

.. _astropy-time-performance:

.. Performance Tips
.. ================
..
.. Here we provide some tips and tricks for how to optimize performance of code
.. using `astropy.time`.
Light Travel Times
------------------

Listed below are two approaches to calculating light travel times for tens of
thousands of sources in a degree patch of the sky. The first approach calculates
the travel times without iteration while the second calculates using iteration.
The second approach using iteration should be avoided as it is extraordinarily slow.

    >>> import numpy as np
    >>> import astropy.coordinates as coord
    >>> import astropy.units as u
    >>> from astropy.time import Time

    >>> ra = np.random.normal(0.0, 1.0, 50000)
    >>> dec = np.random.normal(0.0, 1.0, 50000)

    >>> coos = coord.SkyCoord(ra, dec, unit=u.deg)
    >>> observatory = coord.EarthLocation.of_site('lapalma')
    >>> %time ltts = [Time.light_travel_time(self=Time.now(), skycoord=coos, location=observatory) for coo in coos] # doctest: +SKIP

If you have an internet connection, **coord.EarthLocation.of_site('lapalma')** can be used.
However, if that does not work you, or you don't have an internet connection
The observatory below is the same Lapalma Observatory.

>>> observatory = coord.EarthLocation.from_geocentric(5327448.9957829, -1718665.73869569, 3051566.90295403, unit='m')

Approach Not Using Iteration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This approach takes advantage of the vectorized operations in numpy.

   >>> %time Time.light_travel_time(self=Time.now, skycoord=coos, location=observatory)

   **CPU times:** user 56 ms, sys: 3.21 ms, total: 59.2 ms

   **Wall time:** 58.2 ms

Approach Using Iteration
^^^^^^^^^^^^^^^^^^^^^^^^
Iterating over a SkyCoord object is not fast, due to the how a SkyCoord is constructed. Therefore
the user time is significantly longer.

   >>> observatory = coord.EarthLocation.of_site('lapalma')
   >>> %time ltts = [Time.light_travel_time(self=Time.now(),skycoord=coos, location=observatory) for coo in coos]

   **CPU times:** user 16min 45s, sys: 5.08 s, total: 16min 50s

   **Wall time:** 16min 58s

Avoid Iteration For Light Travel Time for Many Objects
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For user cases where there are thousands of times for each source, broadcasting
can be used:

Broadcasting is a feature of numpy that allows a smaller array to "broadcast"
across the larger array. If you would like to know more about broadcasting go here_.


>>> times = Time.now() + np.linspace(0, 3, 1000)*u.day
>>> coos2 = coos[:,np.newaxis]
>>> coos.shape
   (50000,)
>>> coos2.shape
   (50000, 1)
>>> %time ltts = times.light_travel_time(skycoord=coos2, location=observatory)

   **CPU times:** user 13.8 s, sys: 13 s, total: 26.8 s

   **Wall time:** 27.9 s

.. _here: https://docs.scipy.org/doc/numpy-1.15.0/user/basics.broadcasting.html