Creating Coordinate Objects
---------------------------

Creating new coordinate objects is of course crucial to using
`~astropy.coordinates`.  The typical way to create a new coordinate object
is to directly initialize your preferred coordinate system using standard
python class creation, using the name of the class representing that
system and a number for the two angles.  For example::

    >>> from astropy.coordinates import ICRSCoordinates, FK4Coordinates, GalacticCoordinates
    >>> ICRSCoordinates(187.70592, 12.39112, unit=(u.degree, u.degree))
    <ICRSCoordinates RA=187.70592 deg, Dec=12.39112 deg>
    >>> FK4Coordinates(187.07317, 12.66715, unit=(u.degree, u.degree))
    <FK4Coordinates RA=187.07317 deg, Dec=12.66715 deg>
    >>> GalacticCoordinates(-76.22237, 74.49108, unit=(u.degree, u.degree))
    <GalacticCoordinates l=-76.22237 deg, b=74.49108 deg>

Note that if you do not provide units explicitly, this will fail::

    >>> ICRSCoordinates(23, 1)
    UnitsError: No unit was specified in Angle initializer; the unit parameter should be an object from the  astropy.units module (e.g. 'from astropy import units as u', then use 'u.degree').

While the above example uses python numerical types, you can also provide strings to create coordinates.
If the `unit` parameter is ``(None, None)`` (the default), strings will be interpreted using the `Angle` 
class' parsing scheme, and has a guiding principal of being able to interpret any *unambiguous* string 
specifying an angle. For the exact rules for how each string is parsed, see the 
`~astropy.coordinates.angles.Angle` documentation.  Some examples::

    >>> ICRSCoordinates("3h36m29.7888s -41d08m15.162342s", unit=(None, None))
    <ICRSCoordinates RA=54.12412 deg, Dec=-41.13755 deg>
    >>> ICRSCoordinates("3h36m29.7888s -41d08m15.162342s")
    <ICRSCoordinates RA=54.12412 deg, Dec=-41.13755 deg>
    >>> ICRSCoordinates("14.12412 hours", "-41:08:15.162342 degrees")
    <ICRSCoordinates RA=211.86180 deg, Dec=-41.13755 deg>
    >>> ICRSCoordinates("14.12412 -41:08:15.162342")
    UnitsError: Could not infer Angle units from provided string 14.12412

You can also directly specify the units for both to resolve ambiguities in parsing the angle strings::

    >>> ICRSCoordinates("14.12412 -41:08:15.162342", unit=(u.hour, u.degree))
    <ICRSCoordinates RA=211.86180 deg, Dec=-41.13755 deg>
    >>> ICRSCoordinates("54:7:26.832 -41:08:15.162342", unit=(u.degree, u.degree))
    <ICRSCoordinates RA=54.12412 deg, Dec=-41.13755 deg>
    >>> ICRSCoordinates('3 4 5 +6 7 8', unit=(u.hour, u.degree))
    <ICRSCoordinates RA=46.02083 deg, Dec=6.11889 deg>
    >>> ICRSCoordinates('3h4m5s +6d7m8s', unit=(u.hour, u.degree))
    <ICRSCoordinates RA=46.02083 deg, Dec=6.11889 deg>

This will also give you an error if you give a string with units that conflict with your desired units::

    >>> ICRSCoordinates('3d4m5s +6h7m8s', unit=(u.hour, u.degree))
    ValueError: parse_hours: Invalid input string, can't parse to HMS. (3d4m5s)