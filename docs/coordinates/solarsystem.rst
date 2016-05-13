.. include:: references.txt

.. _astropy-coordinates-solarsystem:

Solar System Ephemerides
----------------------------

`astropy.coordinates` can calculate the |SkyCoord| of some of the major
solar system objects. This functionality requires the 
`jplephem <https://pypi.python.org/pypi/jplephem>`_ package
to be installed. Coordinates are calculated using the JPL DE430 ephemeris file. 
The ephemeris file provides predictions valid for years between 1550 and 2650. 
The file is 115 MB and will be downloaded the first time, but cached after that.

Three functions are provided; :meth:`~astropy.coordinates.get_body`, 
:meth:`~astropy.coordinates.get_moon` and 
:meth:`~astropy.coordinates.get_barycentric_body_position`. The first
two functions return |SkyCoord| objects in the `~astropy.coordinates.GCRS` frame,
whilst the latter returns a `~astropy.coordinates.CartesianRepresentation` of the barycentric position
of a body (i.e in the `~astropy.coordinates.ICRS` frame).

The methods are used as follows::

    >>> from astropy.time import Time
    >>> from astropy.coordinates import get_moon, get_body
    >>> from astropy.coordinates import get_barycentric_body_position, EarthLocation
    >>> t = Time("2014-09-22 23:22")
    >>> loc = EarthLocation.of_site('greenwich')
    >>> get_moon(t, loc) # doctest: +REMOTE_DATA
    <SkyCoord (GCRS: obstime=2014-09-22 23:22:00.000, obsgeoloc=[  3.98060890e+06  -1.02475229e+02   4.96686127e+06] m, obsgeovel=[ 0.  0.  0.] m / s): (ra, dec, distance) in (deg, deg, km)
    (165.51839027, 2.32901144, 407226.55887392)>
    >>> get_body(t, 'jupiter', loc) # doctest: +REMOTE_DATA
    <SkyCoord (GCRS: obstime=2014-09-22 23:22:00.000, obsgeoloc=[  3.98060890e+06  -1.02475229e+02   4.96686127e+06] m, obsgeovel=[ 0.  0.  0.] m / s): (ra, dec, distance) in (deg, deg, km)
    (136.90234741, 17.03160607, 889196019.26282585)>
    >>> get_barycentric_body_position(t, 'moon') # doctest: +REMOTE_DATA
    <CartesianRepresentation (x, y, z) in km
    (150107535.26352832, -866789.03506676, -418963.52113854)>
       
The bodies for which positions can be calculated can be listed::

    >>> from astropy.coordinates import SOLAR_SYSTEM_BODIES
    >>> SOLAR_SYSTEM_BODIES
    ['sun', 'mercury', 'venus', 'earth-moon-barycenter', 'earth', 'moon', 'mars', 'jupiter', 'saturn', 'uranus', 'neptune', 'pluto']


