Separations
-----------

The on-sky separation is easily computed with the
:meth:`astropy.coordinates.BaseCoordinateFrame.separation` or
:meth:`astropy.coordinates.SkyCoord.separation` methods,
which computes the great-circle distance (*not* the small-angle
approximation)::

    >>> from astropy import units as u
    >>> from astropy.coordinates import ICRS
    >>> c1 = ICRS('5h23m34.5s', '-69d45m22s')
    >>> c2 = ICRS('0h52m44.8s', '-72d49m43s')
    >>> sep = c1.separation(c2)
    >>> sep
    <Angle 20.746118437707... deg>


The returned object is an `~astropy.coordinates.Angle` instance, so it
is straightforward to access the angle in any of several equivalent angular
units::

    >>> sep.radian
    0.36208807374669...
    >>> sep.hour
    1.38307456251383...
    >>> sep.arcminute
    1244.76710626244...
    >>> sep.arcsecond
    74686.0263757469...

The :meth:`astropy.coordinates.BaseCoordinateFrame.separation_3d` or
:meth:`astropy.coordinates.SkyCoord.separation_3d` methods can similarly
be used to get 3D distance between two coordinates, rather than on-sky::



    >>> from astropy.coordinates import SkyCoord
    >>> c1 = SkyCoord('5h23m34.5s', '-69d45m22s', distance=70*u.kpc, frame='icrs')
    >>> c2 = SkyCoord('0h52m44.8s', '-72d49m43s', distance=80*u.kpc, frame='icrs')
    >>> sep = c1.separation_3d(c2)
    >>> sep
    <Distance 28.743988157814098 kpc>
