Creating Coordinate Objects
---------------------------

Creating new coordinate objects is of course crucial to using
`~astropy.coordinates`.  The typical way to create a new coordinate object
is to directly initialize your preferred coordinate system using standard
python class creation, using the name of the class representing that
system and a number for the two angles.  For example::

    >>> from astropy.coordinates import ICRS, FK4, Galactic
    >>> import astropy.units as u
    >>> ICRS(187.70592, 12.39112, unit=(u.degree, u.degree))
    <ICRS RA=187.70592 deg, Dec=12.39112 deg>
    >>> FK4(187.07317, 12.66715, unit=(u.degree, u.degree))
    <FK4 RA=187.07317 deg, Dec=12.66715 deg>
    >>> Galactic(283.77763, 74.49108, unit=(u.degree, u.degree))
    <Galactic l=283.77763 deg, b=74.49108 deg>

Note that if you do not provide units explicitly, this will fail::

    >>> ICRS(23, 1)
    Traceback (most recent call last):
        ...
    UnitsError: No unit was specified

While the above example uses python numerical types, you can also
provide strings to create coordinates.  If the `unit` parameter is
``(None, None)`` (the default), strings will be interpreted using the
`Angle` class' parsing scheme, and has a guiding principal of being
able to interpret any *unambiguous* string specifying an angle. For
the exact rules for how each string is parsed, see the
`~astropy.coordinates.angles.Angle` documentation.  Some examples::

    >>> ICRS("3h36m29.7888s -41d08m15.162342s", unit=(None, None))
    <ICRS RA=54.12412 deg, Dec=-41.13755 deg>
    >>> ICRS("3h36m29.7888s -41d08m15.162342s")
    <ICRS RA=54.12412 deg, Dec=-41.13755 deg>
    >>> ICRS("14.12412 hours", "-41:08:15.162342 degrees")
    <ICRS RA=211.86180 deg, Dec=-41.13755 deg>
    >>> ICRS("14.12412 -41:08:15.162342")
    Traceback (most recent call last):
        ...
    UnitsError: No unit specified

You can also directly specify the units for both to resolve
ambiguities in parsing the angle strings::

    >>> ICRS("14.12412 -41:08:15.162342", unit=(u.hour, u.degree))
    <ICRS RA=211.86180 deg, Dec=-41.13755 deg>
    >>> ICRS("54:7:26.832 -41:08:15.162342", unit=(u.degree, u.degree))
    <ICRS RA=54.12412 deg, Dec=-41.13755 deg>
    >>> ICRS('3 4 5 +6 7 8', unit=(u.hour, u.degree))
    <ICRS RA=46.02083 deg, Dec=6.11889 deg>
    >>> ICRS('3h4m5s +6d7m8s', unit=(u.hour, u.degree))
    <ICRS RA=46.02083 deg, Dec=6.11889 deg>

It's also possible to create coordinates using lists or `numpy` arrays.  The same
unit rules apply as for scalar angles.::

    >>> ICRS([187.70592, 123.45678], [12.39112, 9.87654], unit=(u.degree, u.degree))
    <ICRS RA=[ 187.70592  123.45678] deg, Dec=[ 12.39112   9.87654] deg>
    >>> ICRS([187.70592, 123.45678], [12.39112, 9.87654])
    Traceback (most recent call last):
        ...
    UnitsError: No unit was specified

.. warning::
    If you try to create an angle using a tuple for each angle instead of a list or
    `numpy` array, it will be interpreted aa ``(hours, minutes, seconds)`` or
    ``(degrees, arcmin, arcsec)``.  So if you actually want multiple coordinates from
    a tuple, convert it to a list or array.

One final way to create coordinates is to copy them from an already
existing coordinate object::

    >>> i1 = ICRS(187.70592, 12.39112, unit=(u.degree, u.degree))
    >>> i2 = ICRS(i1)
    >>> i1
    <ICRS RA=187.70592 deg, Dec=12.39112 deg>
    >>> i2
    <ICRS RA=187.70592 deg, Dec=12.39112 deg>
