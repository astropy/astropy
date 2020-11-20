.. include:: references.txt

.. _astropy-coordinates-common-errors:

Common mistakes
***************

The following are some common sources of difficulty when using `~astropy.coordinates`.

Object Separation
-----------------

When calculating the separation between objects, it is important to bear in mind that
:meth:`astropy.coordinates.SkyCoord.separation` gives a different answer depending
upon the order in which is used. For example::

    >>> import numpy as np
    >>> from astropy import units as u
    >>> from astropy.coordinates import SkyCoord, GCRS
    >>> from astropy.time import Time
    >>> t = Time("2010-05-22T00:00")
    >>> moon = SkyCoord(104.29*u.deg, 23.51*u.deg, 359367.3*u.km, frame=GCRS(obstime=t))
    >>> star = SkyCoord(101.4*u.deg, 23.02*u.deg, frame='icrs')
    >>> star.separation(moon) # doctest: +FLOAT_CMP
    <Angle 139.84211884 deg>
    >>> moon.separation(star) # doctest: +FLOAT_CMP
    <Angle 2.70390995 deg>

Why do these give such different answers?

The reason is that :meth:`astropy.coordinates.SkyCoord.separation` gives the separation as measured
in the frame of the |SkyCoord| object. So ``star.separation(moon)`` gives the angular separation
in the ICRS frame. This is the separation as it would appear from the Solar System Barycenter. For a
geocentric observer, ``moon.separation(star)`` gives the correct answer, since ``moon`` is in a
geocentric frame.

AltAz calculations for Earth-based objects
------------------------------------------

One might expect that the following code snippet would produce an altitude of exactly 90 degrees::

    >>> from astropy.coordinates import EarthLocation, AltAz
    >>> from astropy.time import Time
    >>> from astropy import units as u

    >>> t = Time('J2010')
    >>> obj = EarthLocation(-1*u.deg, 52*u.deg, height=10.*u.km)
    >>> home = EarthLocation(-1*u.deg, 52*u.deg, height=0.*u.km)
    >>> altaz_frame = AltAz(obstime=t, location=home)
    >>> obj.get_itrs(t).transform_to(altaz_frame).alt # doctest: +FLOAT_CMP
    <Latitude 86.32878441 deg>

Why is the result over three degrees away from the zenith? It is only possible to understand by taking a very careful
look at what ``obj.get_itrs(t)`` returns. This call provides the ITRS position of the source ``obj``. ITRS is
a geocentric coordinate system, and includes the aberration of light. So the code above provides the ITRS position
of a source that appears directly overhead *for an observer at the geocenter*.

Due to aberration, the actual position of this source will be displaced from its apparent position, by an angle of
approximately 20.5 arcseconds. Because this source is about one Earth radius away, it's actual position is therefore
around 600 metres away from where it appears to be. This 600 metre shift, for an object 10 km above the Earth's surface
is an angular difference of just over three degrees - which is why this object does not appear overhead for a topocentric
observer - one on the surface of the Earth.

The correct way to construct a |SkyCoord| object for a source that is directly overhead a topocentric observer is
as follows::

    >>> from astropy.coordinates import EarthLocation, AltAz, ITRS, CIRS
    >>> from astropy.time import Time
    >>> from astropy import units as u

    >>> t = Time('J2010')
    >>> obj = EarthLocation(-1*u.deg, 52*u.deg, height=10.*u.km)
    >>> home = EarthLocation(-1*u.deg, 52*u.deg, height=0.*u.km)

    >>> # Now we make a ITRS vector of a straight overhead object
    >>> itrs_vec = obj.get_itrs(t).cartesian - home.get_itrs(t).cartesian

    >>> # Now we create a topocentric coordinate with this data
    >>> # Any topocentric frame will work, we use CIRS
    >>> # Start by transforming the ITRS vector to CIRS
    >>> cirs_vec = ITRS(itrs_vec, obstime=t).transform_to(CIRS(obstime=t)).cartesian
    >>> # Finally, make CIRS frame object with the correct data
    >>> cirs_topo = CIRS(cirs_vec, obstime=t, location=home)

    >>> # convert to AltAz
    >>> aa = cirs_topo.transform_to(AltAz(obstime=t, location=home))
    >>> aa.alt # doctest: +FLOAT_CMP
    <Latitude 90. deg>