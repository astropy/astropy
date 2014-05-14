Creating Coordinate Objects
---------------------------

TODO: We might want to first introduce SkyCoordinate with many of the
same example as below.

Creating new coordinate objects is of course crucial to using
`~astropy.coordinates`.  The typical way to create a new coordinate object
is to directly initialize your preferred coordinate system using standard
python class creation, using the name of the class representing that
system and a number for the two angles.  For example::

    >>> from astropy.coordinates import ICRS, FK4, Galactic
    >>> import astropy.units as u
    >>> ICRS(187.70592*u.degree, 12.39112*u.degree)
    <ICRS Coordinate: ra=187.70592 deg, dec=12.39112 deg>
    >>> FK4(187.07317*u.degree, 12.66715*u.degree)
    <FK4 Coordinate: equinox=B1950.000, obstime=B1950.000, ra=187.07317 deg, dec=12.66715 deg>
    >>> Galactic(283.77763*u.degree, 74.49108*u.degree)
    <Galactic Coordinate: l=283.77763 deg, b=74.49108 deg>

Note that if you do not provide units explicitly, this will fail::

    >>> ICRS(23, 1)
    Traceback (most recent call last):
        ...
    UnitsError: No unit was given - must be some kind of angle

While the above example uses python numerical types, you can also provide
strings to create coordinates.  Strings will be interpreted using the
`~astropy.coordinates.Angle` class' parsing scheme, and has a guiding
principal of being able to interpret any *unambiguous* string specifying an
angle. For the exact rules for how each string is parsed, see the
`~astropy.coordinates.Angle` documentation.  Some examples::

    >>> ICRS("3h36m29.7888s", "-41d08m15.162342s")
    <SkyCoord (ICRS): ra=54.12412 deg, dec=-41.137545095 deg>
    >>> ICRS("14.12412 hours", "-41:08:15.162342 degrees")
    <SkyCoord (ICRS): ra=211.8618 deg, dec=-41.137545095 deg>
    >>> ICRS("14.12412", "-41:08:15.162342")
    Traceback (most recent call last):
        ...
    UnitsError: No unit specified

It's also possible to create coordinates using lists or `numpy` arrays.  The
same unit rules apply as for scalar angles.::

    >>> ICRS([187.70592, 123.45678]*u.degree, [12.39112, 9.87654]*u.degree)  # doctest: +SKIP
    <ICRS Coordinate: (ra, dec) in deg
        [(187.70592, 12.39112), (123.45677999999998, 9.87654)]>
    >>> ICRS([187.70592*u.degree, 8.23*u.hourangle], [12.39112*u.degree, 1.2*u.radian])  # doctest: +SKIP
    <SkyCoord (ICRS): (ra, dec) in deg
        [(187.70592, 12.39112), (123.44999999999999, 68.75493541569878)]>
    >>> ICRS([187.70592, 123.45678], [12.39112, 9.87654])
    Traceback (most recent call last):
        ...
    UnitsError: No unit was given - must be some kind of angle

.. warning::
    If you try to create an angle using a tuple for each angle instead of a
    list or `numpy` array, it will be interpreted as ``(hours, minutes,
    seconds)`` or ``(degrees, arcmin, arcsec)``.  So if you actually want
    multiple coordinates from a tuple, convert it to a list or array.

TODO: update for SkyCoordinate?
One final way to create coordinates is to copy them from an already
existing coordinate object::
