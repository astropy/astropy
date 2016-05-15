.. include:: references.txt

.. _astropy-coordinates-solarsystem:

Solar System Ephemerides
------------------------

`astropy.coordinates` can calculate the |SkyCoord| of some of the major
solar system objects. Coordinates are calculated using the JPL DE430 
ephemerides. These ephemerides provide predictions valid roughly for years 
between 1550 and 2650. The file is 115 MB and will need to be downloaded 
the first time you use this functionality, but will be cached after that.

.. note::
    This functionality requires that the 
    `jplephem <https://pypi.python.org/pypi/jplephem>`_ package
    is installed. This is most easily acheived via ``pip install jplephem``,
    although whatever package management system you use might have it as well.

Three functions are provided; :meth:`~astropy.coordinates.get_body`, 
:meth:`~astropy.coordinates.get_moon` and 
:meth:`~astropy.coordinates.get_body_barycentric`. The first
two functions return |SkyCoord| objects in the `~astropy.coordinates.GCRS` frame,
whilst the latter returns a `~astropy.coordinates.CartesianRepresentation` of the barycentric position
of a body (i.e in the `~astropy.coordinates.ICRS` frame).

Here are some examples of these functions in use::

    >>> from astropy.time import Time
    >>> from astropy.coordinates import get_moon, get_body
    >>> from astropy.coordinates import get_body_barycentric, EarthLocation
    >>> t = Time("2014-09-22 23:22")
    >>> loc = EarthLocation.of_site('greenwich')
    >>> get_moon(t, loc) # doctest: +REMOTE_DATA
    <SkyCoord (GCRS: obstime=2014-09-22 23:22:00.000, obsgeoloc=[ 3949481.69039034  -550931.90976401  4961151.73716876] m, obsgeovel=[  40.17459314  288.00078055   -0.        ] m / s): (ra, dec, distance) in (deg, deg, km)
        (165.51840736, 2.32900633, 407226.68749637)>
    >>> get_body(t, 'jupiter', loc) # doctest: +REMOTE_DATA
    <SkyCoord (GCRS: obstime=2014-09-22 23:22:00.000, obsgeoloc=[ 3949481.69039034  -550931.90976401  4961151.73716876] m, obsgeovel=[  40.17459314  288.00078055   -0.        ] m / s): (ra, dec, distance) in (deg, deg, km)
        (136.90234741, 17.03160607, 889196019.26282585)>
    >>> get_body_barycentric(t, 'moon') # doctest: +REMOTE_DATA
    <CartesianRepresentation (x, y, z) in km
    (150107535.26352832, -866789.03506676, -418963.52113854)>
       
For a list of the bodies for which positions can be calculated, do::

    >>> from astropy.coordinates import SOLAR_SYSTEM_BODIES
    >>> SOLAR_SYSTEM_BODIES
    ('sun', 'mercury', 'venus', 'earth-moon-barycenter', 'earth', 'moon', 'mars', 'jupiter', 'saturn', 'uranus', 'neptune', 'pluto')

.. note ::
    While the sun is included in the these ephemerides, it is important to
    recognize that `~astropy.coordinates.get_sun` does *not* use this
    method, but instead uses a polynomial model for the location of the sun
    (as this requires no special download). So it is not safe to assume that
    ``get_body(time, 'sun')`` and ``get_sun(time)`` will give the same result.

You can also change the SPK kernel (the file used to actually locate the 
planets), although this interface should be considered preliminary (and hence 
is not yet considered part of the public API)::

    >>> from astropy import coordinates
    >>> with coordinates.solar_system.kernel_url.set('http://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de432s.bsp'):
    ...     coordinates.get_body(t,'pluto') # doctest: +REMOTE_DATA
    <SkyCoord (GCRS: obstime=2014-09-22 23:22:00.000, obsgeoloc=[ 0.  0.  0.] m, obsgeovel=[ 0.  0.  0.] m / s): (ra, dec, distance) in (deg, deg, km)
        (281.52508175, -20.60080214, 4865115955.7188015)>
