Working with Angles
-------------------

The angular components of a coordinate are represented by objects of the
`~astropy.coordinates.angles.Angle` class. These objects can be instantiated on
their own anywhere a representation of an angle is needed, and support a variety
of ways of representing the value of the angle::

    >>> from astropy.coordinates import Angle
    >>> a = Angle(1, u.radian)
    >>> a
    <astropy.coordinates.angles.Angle 57.29578 deg>
    >>> a.radians
    1
    >>> a.degrees
    57.29577951308232
    >>> a.hours
    3.819718634205488
    >>> a.hms
    (3.0, 49, 10.987083139757061)
    >>> a.dms
    (57.0, 17, 44.80624709636231)
    >>> a.format()
    '57d17m44.80625s'
    >>> a.format(sep=':')
    '57:17:44.80625'
    >>> a.format(sep=('deg','m','s'))
    '57deg17m44.80625s'
    >>> a.format(u.hour)
    '3h49m10.98708s'
    >>> a.format(u.radian)
    '1.0radian'
    >>> a.format(u.radian, decimal=True)
    '1.0'

`~astropy.corodinates.angles.Angle` objects can also have bounds.  These specify
either a limited range in which the angle is valid (if it's <360 degrees), or
the limit at which an angle is wrapped back around to 0::

    >>> Angle(90, unit=u.degree, bounds=(0,180))
    <Angle 90.00000 deg>
    >>> Angle(-270, unit=u.degree, bounds=(0,180))
    <Angle 90.00000 deg>
    >>> Angle(181, unit=u.degree, bounds=(0,180))
    BoundsError: The angle given falls outside of the specified bounds.
    >>> Angle(361, unit=u.degree, bounds=(0,360))
    <Angle 1.00000 deg>

Angles will also behave correctly for appropriate arithmetic operations::

    >>> a = Angle(1, u.radian)
    >>> a + a
    <Angle 114.59156 deg>
    >>> a - a
    <Angle 0.00000 deg>
    >>> a == a
    True
    >>> a == (a + a)
    False

Angle objects can also be used for creating coordinate objects::

    >>> ICRSCoordinates(Angle(1, u.radian), Angle(2, u.radian))
    <ICRSCoordinates RA=57.29578 deg, Dec=114.59156 deg>
    >>> ICRSCoordinates(RA(1, u.radian), Dec(2, u.radian))
    <ICRSCoordinates RA=57.29578 deg, Dec=114.59156 deg>