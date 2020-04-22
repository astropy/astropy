.. include:: references.txt

.. _astropy-coordinates-satellites:

Working with Earth satellites using Astropy Coordinates
*******************************************************

Finding ``TEME`` coordinates from Two-Line-Ephemerides (TLE)
============================================================

There is currently no support in `astropy.coordinates` for computing satellite orbits
from the TLE orbital elements avaiable from a source like `Celestrak <http://celestrak.com/>`_.
Full support for handling TLE files is available in the `Skyfield <http://rhodesmill.org/skyfield/>`_
library, but if you do not wish to install that library, advice for dealing with satellite data in
astropy is below.

..
  EXAMPLE START
  Using sgp4 to get a TEME coordinate

You will need some external library to compute the position and velocity of the satellite from the
TLE orbital elements. The `sgp4 <https://pypi.org/project/sgp4/>`_ library can do this. An example
using this library to find the  `~astropy.coordinates.TEME` coordinates of a satellite is::

.. doctest-requires:: sgp4

        >>> from sgp4.api import Satrec
        >>> from sgp4.api import SGP4_ERRORS
        >>> s = '1 25544U 98067A   19343.69339541  .00001764  00000-0  38792-4 0  9991'
        >>> t = '2 25544  51.6439 211.2001 0007417  17.6667  85.6398 15.50103472202482'
        >>> satellite = Satrec.twoline2rv(s, t)

The ``satellite`` object has a method, ``satellite.sgp4``, that will try to compute the TEME position
and velocity at a given time::

.. doctest-requires:: sgp4

        >>> from astropy.time import Time
        >>> t = Time(2458827.362605, format='jd')
        >>> error_code, teme_p, teme_v = satellite.sgp4(t.jd1, t.jd2)  # in km and km/s
        >>> if error_code != 0:
        ...     raise RuntimeError(SGP4_ERRORS[error_code])

Now we have the position and velocity in kilometer and kilometers per second, we can create a
position in the `~astropy.coordinates.TEME` reference frame::

.. doctest-requires:: sgp4

        >>> from astropy.coordinates import TEME, CartesianDifferential, CartesianRepresentation
        >>> from astropy import units as u
        >>> teme_v = CartesianDifferential(teme_v*u.km/u.s)
        >>> teme_p = CartesianRepresentation(*teme_p*u.km, differentials={'s': teme_v})
        >>> teme = TEME(teme_p, obstime=t)

Note how we are careful to set the observed time of the `~astropy.coordinates.TEME` frame to
the time at which we calculated satellite position.

Transforming ``TEME`` to other coordinate systems
==================================================

Once you have satellite positions in `~astropy.coordinates.TEME` coordinates they can be easily transformed
into any `astropy.coordinates` frame.

For example, to find the overhead latitude, longitude and height of the satellite::

.. doctest-requires:: sgp4

        >>> from astropy.coordinates import ITRS
        >>> itrs = teme.transform_to(ITRS(obstime=t))  # doctest: +REMOTE_DATA
        >>> location = itrs.earth_location  # doctest: +REMOTE_DATA
        >>> location.geodetic  # doctest: +REMOTE_DATA +FLOAT_CMP
        GeodeticLocation(lon=<Longitude 160.34199789 deg>, lat=<Latitude -24.6609379 deg>, height=<Quantity 420.17927591 km>)

Or, if you want to find the altitude and azimuth of the satellite from a particular location::

.. doctest-requires:: sgp4

        >>> from astropy.coordinates import EarthLocation, AltAz
        >>> siding_spring = EarthLocation.of_site('aao')  # doctest: +REMOTE_DATA
        >>> aa = teme.transform_to(AltAz(obstime=t, location=siding_spring))  # doctest: +REMOTE_DATA
        >>> aa.alt  # doctest: +REMOTE_DATA +FLOAT_CMP
        <Latitude 10.94798427 deg>
        >>> aa.az  # doctest: +REMOTE_DATA +FLOAT_CMP
        <Longitude 59.28807348 deg>

..
  EXAMPLE END
