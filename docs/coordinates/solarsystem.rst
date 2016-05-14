.. include:: references.txt

.. _astropy-coordinates-solarsystem:

Solar System Ephemerides
------------------------

`astropy.coordinates` can calculate the |SkyCoord| of some of the major
solar system objects. This functionality requires the 
`jplephem <https://pypi.python.org/pypi/jplephem>`_ package
to be installed. Coordinates are calculated using the JPL DE430 ephemeris file. 
The ephemeris file provides predictions valid for years between 1550 and 2650. 
The file is 115 MB and will be downloaded the first time, but cached after that.

Three functions are provided; :meth:`~astropy.coordinates.get_body`, 
:meth:`~astropy.coordinates.get_moon` and 
:meth:`~astropy.coordinates.get_body_barycentric`. The first
two functions return |SkyCoord| objects in the `~astropy.coordinates.GCRS` frame,
whilst the latter returns a `~astropy.coordinates.CartesianRepresentation` of the barycentric position
of a body (i.e in the `~astropy.coordinates.ICRS` frame).

The methods are used as follows::

    >>> from astropy.time import Time
    >>> from astropy.coordinates import get_moon, get_body
    >>> from astropy.coordinates import get_body_barycentric, EarthLocation
    >>> t = Time("2014-09-22 23:22")
    >>> loc = EarthLocation.of_site('greenwich')
    >>> get_moon(t, loc) # doctest: +REMOTE_DATA
    <SkyCoord (GCRS: obstime=2014-09-22 23:22:00.000, obsgeoloc=[ 3949481.69039034  -550931.90976401  4961151.73716876] m, obsgeovel=[  40.17459314  288.00078055   -0.        ] m / s): (ra, dec, distance) in (deg, deg, km)
    (165.51839027, 2.32901144, 407226.55887392)>
    >>> get_body(t, 'jupiter', loc) # doctest: +REMOTE_DATA
    <SkyCoord (GCRS: obstime=2014-09-22 23:22:00.000, obsgeoloc=[ 3949481.69039034  -550931.90976401  4961151.73716876] m, obsgeovel=[  40.17459314  288.00078055   -0.        ] m / s): (ra, dec, distance) in (deg, deg, km)
    (136.90234741, 17.03160607, 889196019.26282585)>
    >>> get_body_barycentric(t, 'moon') # doctest: +REMOTE_DATA
    <CartesianRepresentation (x, y, z) in km
    (150107535.26352832, -866789.03506676, -418963.52113854)>
       
The bodies for which positions can be calculated can be listed::

    >>> from astropy.coordinates import SOLAR_SYSTEM_BODIES
    >>> SOLAR_SYSTEM_BODIES
    ('sun', 'mercury', 'venus', 'earth-moon-barycenter', 'earth', 'moon', 'mars', 'jupiter', 'saturn', 'uranus', 'neptune', 'pluto')


