.. note that if this is changed from the default approach of using an *include*
   (in index.rst) to a separate performance page, the header needs to be changed
   from === to ***, the filename extension needs to be changed from .inc.rst to
   .rst, and a link needs to be added in the subpackage toctree

.. _astropy-coordinates-performance:

Performance Tips
================

If you are using |skycoord| for many different coordinates, you will see much
better performance if you create a single |skycoord| with arrays of coordinates
as opposed to creating individual |skycoord| objects for each individual
coordinate::

    >>> coord = SkyCoord(ra_array, dec_array, unit='deg')  # doctest: +SKIP

In addition, looping over a |skycoord| object can be slow. If you need to
transform the coordinates to a different frame, it is much faster to transform a
single |skycoord| with arrays of values as opposed to looping over the
|skycoord| and transforming them individually.

Finally, for more advanced users, note that you can use broadcasting to
transform |skycoord| objects into frames with vector properties. For example::

    >>> from astropy.coordinates import SkyCoord, EarthLocation
    >>> from astropy import coordinates as coord
    >>> from astropy.coordinates.tests.utils import randomly_sample_sphere
    >>> from astropy.time import Time
    >>> from astropy import units as u
    >>> import numpy as np

    >>> # 1000 random locations on the sky
    >>> ra, dec, _ = randomly_sample_sphere(1000)
    >>> coos = SkyCoord(ra, dec)

    >>> # 300 times over the space of 10 hours
    >>> times = Time.now() + np.linspace(-5, 5, 300)*u.hour

    >>> # note the use of broadcasting so that 300 times are broadcast against 1000 positions
    >>> lapalma = EarthLocation.from_geocentric(5327448.9957829, -1718665.73869569, 3051566.90295403, unit='m')
    >>> aa_frame = coord.AltAz(obstime=times[:, np.newaxis], location=lapalma)

    >>> # calculate alt-az of each object at each time.
    >>> aa_coos = coos.transform_to(aa_frame)  # doctest: +REMOTE_DATA
