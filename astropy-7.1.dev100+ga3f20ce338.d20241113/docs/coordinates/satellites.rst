.. _astropy-coordinates-satellites:

Working with Earth Satellites Using Astropy Coordinates
*******************************************************
This document discusses Two-Line Element ephemerides and the True Equator, Mean Equinox frame.
For satellite ephemerides given directly in geocentric ITRS coordinates
(e.g. `ILRS ephemeris format <https://ilrs.gsfc.nasa.gov/data_and_products/formats/cpf.html>`_)
please see the example transform to `~astropy.coordinates.AltAz` below starting with
the geocentric ITRS coordinate frame.

Satellite data is normally provided in the Two-Line Element (TLE) format
(see `here <https://celestrak.org/NORAD/documentation/tle-fmt.php>`_
for a definition). These datasets are designed to be used in combination
with a theory for orbital propagation model to predict the positions
of satellites.

The history of such models is discussed in detail in
`Vallado et al (2006) <https://celestrak.org/publications/AIAA/2006-6753/AIAA-2006-6753-Rev2.pdf>`_
who also provide a reference implementation of the SGP4 orbital propagation
code, designed to be compatible with the TLE sets provided by the United
States Department of Defense, which are available from a source like
`Celestrak <http://celestrak.org/>`_.

The output coordinate frame of the SGP4 model is the True Equator, Mean Equinox
frame (TEME), which is one of the frames built-in to `astropy.coordinates`.
TEME is an Earth-centered inertial frame (i.e., it does not rotate with respect
to the stars). Several definitions exist; ``astropy`` uses the implementation described
in `Vallado et al (2006) <https://celestrak.org/publications/AIAA/2006-6753/AIAA-2006-6753-Rev2.pdf>`_.

Finding TEME Coordinates from TLE Data
======================================

There is currently no support in `astropy.coordinates` for computing satellite orbits
from TLE orbital element sets. Full support for handling TLE files is available in
the `Skyfield <https://rhodesmill.org/skyfield/>`_ library, but some advice for dealing
with satellite data in ``astropy`` is below.

.. EXAMPLE START Using sgp4 to get a TEME coordinate

You will need some external library to compute the position and velocity of the satellite from the
TLE orbital elements. The `SGP4 <https://pypi.org/project/sgp4/>`_ library can do this. An example
of using this library to find the  `~astropy.coordinates.TEME` coordinates of a satellite is:

.. doctest-requires:: sgp4

    >>> from sgp4.api import Satrec
    >>> from sgp4.api import SGP4_ERRORS
    >>> s = '1 25544U 98067A   19343.69339541  .00001764  00000-0  38792-4 0  9991'
    >>> t = '2 25544  51.6439 211.2001 0007417  17.6667  85.6398 15.50103472202482'
    >>> satellite = Satrec.twoline2rv(s, t)

The ``satellite`` object has a method, ``satellite.sgp4``, that will try to compute the TEME position
and velocity at a given time:

.. doctest-requires:: sgp4

    >>> from astropy.time import Time
    >>> t = Time(2458827.362605, format='jd')
    >>> error_code, teme_p, teme_v = satellite.sgp4(t.jd1, t.jd2)  # in km and km/s
    >>> if error_code != 0:
    ...     raise RuntimeError(SGP4_ERRORS[error_code])

Now that we have the position and velocity in kilometers and kilometers per second, we can create a
position in the `~astropy.coordinates.TEME` reference frame:

.. doctest-requires:: sgp4

    >>> from astropy.coordinates import TEME, CartesianDifferential, CartesianRepresentation
    >>> from astropy import units as u
    >>> teme_p = CartesianRepresentation(teme_p*u.km)
    >>> teme_v = CartesianDifferential(teme_v*u.km/u.s)
    >>> teme = TEME(teme_p.with_differentials(teme_v), obstime=t)

.. EXAMPLE END

Note how we are careful to set the observed time of the `~astropy.coordinates.TEME` frame to
the time at which we calculated satellite position.

Transforming TEME to Other Coordinate Systems
=============================================

Once you have satellite positions in `~astropy.coordinates.TEME` coordinates they can be transformed
into any `astropy.coordinates` frame.

For example, to find the overhead latitude, longitude, and height of the satellite:

.. EXAMPLE START Transforming TEME

.. doctest-requires:: sgp4

    >>> from astropy.coordinates import ITRS
    >>> itrs_geo = teme.transform_to(ITRS(obstime=t))
    >>> location = itrs_geo.earth_location
    >>> location.geodetic  # doctest: +FLOAT_CMP
    GeodeticLocation(lon=<Longitude 160.34199789 deg>, lat=<Latitude -24.6609379 deg>, height=<Quantity 420.17927591 km>)

.. testsetup::

    >>> from astropy.coordinates import EarthLocation
    >>> siding_spring = EarthLocation(-4680888.60272112, 2805218.44653429, -3292788.0804506, unit='m')

Or, if you want to find the altitude and azimuth of the satellite from a particular location:

.. note::
    In this example, the intermediate step of manually setting up a topocentric `~astropy.coordinates.ITRS`
    frame is necessary in order to avoid the change in stellar aberration that would occur if a direct
    transform from geocentric to topocentric coordinates using ``transform_to`` was used. Please see
    the documentation of the `~astropy.coordinates.ITRS` frame for more details.

.. doctest-requires:: sgp4

    >>> from astropy.coordinates import EarthLocation, AltAz
    >>> siding_spring = EarthLocation.of_site('aao')  # doctest: +SKIP
    >>> topo_itrs_repr = itrs_geo.cartesian.without_differentials() - siding_spring.get_itrs(t).cartesian
    >>> itrs_topo = ITRS(topo_itrs_repr, obstime = t, location=siding_spring)
    >>> aa = itrs_topo.transform_to(AltAz(obstime=t, location=siding_spring))
    >>> aa.alt  # doctest: +FLOAT_CMP
    <Latitude 10.94799670 deg>
    >>> aa.az  # doctest: +FLOAT_CMP
    <Longitude 59.28803392 deg>

For a stationary observer, velocity in the `~astropy.coordinates.ITRS` is independent of location,
so if you want to carry the velocity to the topocentric frame, you can do so as follows:

.. doctest-requires:: sgp4

    >>> itrs_geo_p = itrs_geo.cartesian.without_differentials()
    >>> itrs_geo_v = itrs_geo.cartesian.differentials['s']
    >>> topo_itrs_p = itrs_geo_p - siding_spring.get_itrs(t).cartesian
    >>> topo_itrs_repr = topo_itrs_p.with_differentials(itrs_geo_v)
    >>> itrs_topo = ITRS(topo_itrs_repr, obstime = t, location=siding_spring)

.. EXAMPLE END
