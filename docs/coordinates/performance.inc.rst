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
transform |skycoord| objects into frames with vector properties.

Example
-------

..
  EXAMPLE START
  Performance Tips for Transforming SkyCoord Objects

To use broadcasting to transform |skycoord| objects into frames with vector
properties::

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
    >>> aa_coos = coos.transform_to(aa_frame)  # doctest: +REMOTE_DATA +IGNORE_WARNINGS

..
  EXAMPLE END

Improving Performance for Arrays of ``obstime``
===============================================

The most expensive operations when transforming between observer-dependent coordinate
frames (e.g. ``AltAz``) and sky-fixed frames (e.g. ``ICRS``) are the calculation
of the orientation and position of Earth.

If |skycoord| instances with a large ``obstime`` array are transformed,
these calculations can be sped up by factors up to 100, whilst still keeping micro-arcsecond precision,
by utilizing interpolation instead of calculating Earth orientation parameters for each individual point.

This can be achieved through the ``erfa_astrom`` state and the ``ErfaAstromInterpolator``
class like this

..
  EXAMPLE START
  Improving performance for obstime arrays

To use interpolation for the astrometric values in coordinate transformation, use::

   >>> from astropy.coordinates import SkyCoord, EarthLocation, AltAz
   >>> from astropy.coordinates.erfa_astrom import erfa_astrom, ErfaAstromInterpolator
   >>> from astropy.time import Time
   >>> from time import perf_counter
   >>> import numpy as np
   >>> import astropy.units as u


   >>> # array with 10000 obstimes
   >>> obstime = Time.now() + np.linspace(0, 6, 10000) * u.hour
   >>> location = location = EarthLocation(lon=-17.89 * u.deg, lat=28.76 * u.deg, height=2200 * u.m)
   >>> frame = AltAz(obstime=obstime, location=location)  # doctest:+REMOTE_DATA
   >>> crab = SkyCoord.from_name('Crab')  # doctest:+REMOTE_DATA

   >>> # transform with default transformation and print duration
   >>> t0 = perf_counter()
   >>> crab_altaz = crab.transform_to(frame)  # doctest:+REMOTE_DATA +IGNORE_WARNINGS
   >>> print(f'Transformation took {perf_counter() - t0:.2f} s')  # doctest:+ELLIPSIS
   Transformation took ... s

   >>> # transform with interpolating astrometric values
   >>> t0 = perf_counter()
   >>> with erfa_astrom.set(ErfaAstromInterpolator(300 * u.s)):  # doctest:+REMOTE_DATA
   ...     crab_altaz_interpolated = crab.transform_to(frame)  # doctest:+REMOTE_DATA +IGNORE_WARNINGS
   >>> print(f'Transformation took {perf_counter() - t0:.2f} s')  # doctest:+ELLIPSIS
   Transformation took ... s

   >>> err = crab_altaz.separation(crab_altaz_interpolated)  # doctest:+REMOTE_DATA
   >>> print(f'Mean error of interpolation: {err.to(u.microarcsecond).mean():.4f}')  # doctest:+ELLIPSIS +REMOTE_DATA
   Mean error of interpolation: 0.0063 uarcsec

   >>> # To set erfa_astrom for a whole session, use it without context manager:
   >>> erfa_astrom.set(ErfaAstromInterpolator(300 * u.s))

..
  EXAMPLE END
