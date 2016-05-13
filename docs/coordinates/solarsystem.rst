.. include:: references.txt

.. _astropy-coordinates-solarsystem:

Solar System Ephemerides
----------------------------

`astropy.coordinates` can calculate the |SkyCoord| of some of the major
solar system objects, using the JPL DE430 ephemeris file. This ephemeris
file provides predictions valid for years between 1550 and 2650. The file is
115 MB and will be downloaded the first time, but cached after that.

Three functions are provided; :meth:`~astropy.coordinates.solar_system.get_body`, 
:meth:`~astropy.coordinates.solar_system.get_moon` and 
:meth:`~astropy.coordinates.solar_system.get_barycentric_body_position`. The first
two functions return |SkyCoord| objects in the `~astropy.coordinates.GCRS` frame,
whilst the latter returns a |CartesianRepresentation| of the barycentric position
of a body (i.e in the `~astropy.coordinates.ICRS` frame).

The methods are used as follows::

    >>> from astropy.time import Time
    >>> from astropy.coordinates import get_moon, get_body
    >>> from astropy.coordinates import get_barycentric_body_position, EarthLocation
    >>> t = Time.now()
    >>> loc = EarthLocation.of_site('greenwich')
    >>> get_moon(t, loc)
    <SkyCoord (GCRS: obstime=2016-03-29 06:35:36.927857, obsgeoloc=[-1463969.30185172 -5166673.34223433  3434985.71204565] m, obsgeovel=[ 0.  0.  0.] m / s): (ra, dec, distance) in (deg, deg, km)
    (250.94788165, -17.04585998, 400244.30166804)>
    >>> get_body(t, 'jupiter', loc)
    <SkyCoord (GCRS: obstime=2016-03-29 06:36:28.310851, obsgeoloc=[-1463969.30185172 -5166673.34223433  3434985.71204565] m, obsgeovel=[ 0.  0.  0.] m / s): (ra, dec, distance) in (deg, deg, km)
    (167.22360841, 7.07220939, 673151288.66240811)>
    >>> get_barycentric_body_position(t, 'moon')
    <CartesianRepresentation (x, y, z) in km
       (-90980973.95347136, -110026376.34185831, -47742071.3962695)>


