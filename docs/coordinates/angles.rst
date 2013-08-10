Working with Angles
-------------------

The angular components of a coordinate are represented by objects of the
`~astropy.coordinates.angles.Angle` class. These objects can be instantiated on
their own anywhere a representation of an angle is needed, and support a variety
of ways of representing the value of the angle::

    >>> from astropy import units as u
    >>> from astropy.coordinates import Angle
    >>> a = Angle(1, u.radian)
    >>> a
    <Angle 1 rad>
    >>> a.radian
    1
    >>> a.degree
    57.29577951308232
    >>> a.hour
    3.8197186342054885
    >>> a.hms
    (3.0, 49, 10.987083139757061)
    >>> a.dms
    (57.0, 17, 44.80624709636231)
    >>> a.arminute
    3437.7467707849396
    >>> a.to_string()
    '1.0rad'
    >>> a.to_string(unit=u.degree)
    '57d17m44.80625s'
    >>> a.to_string(unit=u.degree, sep=':')
    '57:17:44.80625'
    >>> a.to_string(unit=u.degree, sep=('deg', 'm', 's'))
    '57deg17m44.80625s'
    >>> a.to_string(unit=u.hour)
    '3h49m10.98708s'
    >>> a.to_string(unit=u.hour, decimal=True)
    '3.8197'

`~astropy.corodinates.angles.Angle` objects can also have bounds.
These specify either a limited range in which the angle is valid (if
it's < 360 degrees), or the limit at which an angle is wrapped back
around to 0::

    >>> Angle(90, unit=u.degree, bounds=(0, 180))
    <Angle 90.0 deg>
    >>> Angle(-270, unit=u.degree, bounds=(0, 180))
    <Angle 90.0 deg>
    >>> Angle(181, unit=u.degree, bounds=(0, 180))
    astropy.coordinates.errors.BoundsError: The angle(s) 181.0 falls
    outside of the specified bounds (0.0, 180.0)
    >>> Angle(361, unit=u.degree, bounds=(0, 360))
    <Angle 1.0 deg>

Angles will also behave correctly for appropriate arithmetic operations::

    >>> a = Angle(1, u.radian)
    >>> a + a
    <Angle 2 rad>
    >>> a - a
    <Angle 0 rad>
    >>> a == a
    True
    >>> a == (a + a)
    False

Angle objects can also be used for creating coordinate objects::

    >>> from astropy.coordinates import ICRSCoordinates
    >>> ICRSCoordinates(Angle(1, u.radian), Angle(0.5, u.radian))
    <ICRSCoordinates RA=57.29578 deg, Dec=28.64789 deg>
    >>> ICRSCoordinates(RA(1, u.radian), Dec(2, u.radian))
    <ICRSCoordinates RA=57.29578 deg, Dec=28.64789 deg>
