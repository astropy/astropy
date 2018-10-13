.. note that if this is changed from the default approach of using an *include* 
   (in index.rst) to a separate performance page, the header needs to be changed
   from === to ***, the filename extension needs to be changed from .inc.rst to 
   .rst, and a link needs to be added in the subpackage toctree

.. _astropy-time-performance:

Performance Tips
================

Here we provide some tips and tricks for how to optimize performance of code
using `astropy.time`.

Listed below are two approaches to calculating light travel times for tens of
thousands of sources in a degree patch of the sky.

>>> import numpy as np
>>> import astropy.coordinates as coord
>>> import astropy.units as u
>>> from astropy.time import Time

>>> ra = np.random.normal(0.0, 1.0, 50000)
>>> dec = np.random.normal(0.0, 1.0, 50000)

>>> coos = coord.SkyCoord(ra, dec, unit=u.deg)

   Lapalma

>>> observatory = coord.EarthLocation.from_geocentric(5327448.9957829, -1718665.73869569, 3051566.90295403, unit='m')

*The first approach and its time to completion below:*

>>> Time.light_travel_time(self, skycoord=coos, location=observatory)

   **CPU times:** user 56 ms, sys: 3.21 ms, total: 59.2 ms

   **Wall time:** 58.2 ms

*The second approach and its time to completion below:*

>>> ltts = [Time.light_travel_time(coo, location=coos.location) for coo in coos]


   **CPU times:** user 16min 45s, sys: 5.08 s, total: 16min 50s
   **Wall time:** 16min 58s

As you can see, the user time is significantly longer. This is because iteration
is being used.

**AVOID ITERATION FOR LIGHT TRAVEL TIME FOR MANY OBJECTS**

*For user cases where there are thousands of times for each source, broadcasting can be used:*

>>> times = Time.now() + np.linspace(0, 3, 1000)*u.day
>>> ltts = times.light_travel_time(coos[:, np.newaxis], location=observatory)

   **CPU times:** user 13.8 s, sys: 13 s, total: 26.8 s

   **Wall time:** 27.9 s