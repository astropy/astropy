.. |Angle| replace:: `~astropy.coordinates.angles.Angle`
.. |Longitude| replace:: `~astropy.coordinates.angles.Longitude`
.. |Latitude| replace:: `~astropy.coordinates.angles.Latitude`

Working with Angles
-------------------

The angular components of a coordinate are represented by objects of the |Angle|
class. These objects can be instantiated on their own wherever a representation of an
angle is needed.

Creation
^^^^^^^^

The creation of an |Angle| object is quite flexible and supports a wide variety of
input object types and formats.  The type of the input angle(s) can an array, scalar,
tuple, string, `~astropy.units.Quantity` or another |Angle|.  This is best illustrated with a number of
examples of valid ways to create an |Angle|::

    >>> import numpy as np
    >>> from astropy import units as u
    >>> from astropy.coordinates import Angle

    >>> Angle('10.2345d')              # String with 'd' abbreviation for degrees
    <Angle 10.2345 deg>
    >>> Angle(['10.2345d', '-20d'])    # Array of strings
    <Angle [ 10.2345,-20.    ] deg>
    >>> Angle('1:2:30.43 degrees')     # Sexagesimal degrees
    <Angle 1.041786111... deg>
    >>> Angle('1 2 0 hours')           # Sexagesimal hours
    <Angle 1.0333... hourangle>
    >>> Angle(np.arange(1, 8), unit=u.deg)  # Numpy array from 0..7 in degrees
    <Angle [1,2,3,4,5,6,7] deg>
    >>> Angle(u'1°2′3″')               # Unicode degree, arcmin and arcsec symbols
    <Angle 1.0341666... deg>
    >>> Angle('1d2m3.4s')              # Degree, arcmin, arcsec.
    <Angle 1.0342777... deg>
    >>> Angle('-1h2m3s')               # Hour, minute, second
    <Angle -1.0341666... hourangle>
    >>> Angle((-1, 2, 3), unit=u.deg)  # (degree, arcmin, arcsec)
    <Angle -1.0341666... deg>
    >>> Angle(10.2345 * u.deg)         # From a Quantity object in degrees
    <Angle 10.2345 deg>
    >>> Angle(Angle(10.2345 * u.deg))  # From another Angle object
    <Angle 10.2345 deg>


Representation
^^^^^^^^^^^^^^

The |Angle| object also supports a variety of ways of representing the value of the angle,
both as a floating point number as a string::

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
    (3.0, 49.0, 10.987083...)
    >>> a.dms
    (57.0, 17.0, 44.806247...)
    >>> a.arcminute
    3437.746770...
    >>> a.to_string()
    u'1rad'
    >>> a.to_string(unit=u.degree)
    u'57d17m44.8062s'
    >>> a.to_string(unit=u.degree, sep=':')
    u'57:17:44.8062'
    >>> a.to_string(unit=u.degree, sep=('deg', 'm', 's'))
    u'57deg17m44.8062s'
    >>> a.to_string(unit=u.hour)
    u'3h49m10.9871s'
    >>> a.to_string(unit=u.hour, decimal=True)
    u'3.81972'


Usage
^^^^^

Angles will also behave correctly for appropriate arithmetic operations::

    >>> a = Angle(1.0, u.radian)
    >>> a + 0.5 * u.radian + 2 * a
    <Angle 3.5 rad>
    >>> np.sin(a / 2)
    <Quantity 0.479425...>
    >>> a == a
    array(True, dtype=bool)
    >>> a == (a + a)
    array(False, dtype=bool)

|Angle| objects can also be used for creating coordinate objects::

    >>> from astropy.coordinates import ICRS
    >>> ICRS(Angle(1, u.radian), Angle(0.5, u.radian))
    <ICRS RA=57.29578 deg, Dec=28.64789 deg>


Wrapping and bounds
^^^^^^^^^^^^^^^^^^^

There are two utility methods that simplify working with angles that should
have bounds.  The `~astropy.coordinates.angles.Angle.wrap_at()` method allows
taking an angle or angles and wrapping to be within a single 360 degree slice.
The `~astropy.coordinates.angles.Angle.is_within_bounds()` method returns a
boolean indicating whether an angle or angles is within the specified bounds.


Longitude and Latitude objects
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

|Longitude| and |Latitude| are two specialized subclasses of the |Angle| class that are
used for all of the spherical coordinate classes.  |Longitude| is used to represent values
like right ascension, Galactic longitude, and azimuth (for ecliptic, Galactic, and Alt-Az
coordinates, respectively).  |Latitude| is used for declination, Galactic latitude, and
elevation.

Longitude
"""""""""

A |Longitude| object is distinguished from a pure |Angle| by virtue
of a ``wrap_angle`` property.  The ``wrap_angle`` specifies that all angle values
represented by the object will be in the range::

  wrap_angle - 360 * u.deg <= angle(s) < wrap_angle

The default ``wrap_angle`` is 360 deg.  Setting ``wrap_angle=180 * u.deg`` would
instead result in values between -180 and +180 deg.  Setting the ``wrap_angle``
attribute of an existing ``Longitude`` object will result in re-wrapping the
angle values in-place.  For example::

    >>> from astropy.coordinates import Longitude
    >>> a = Longitude([-20, 150, 350, 360] * u.deg)
    >>> a.degree
    array([340, 150, 350,   0])
    >>> a.wrap_angle = 180 * u.deg
    >>> a.degree
    array([-20, 150, -10,   0])

Latitude
""""""""

A Latitude object is distinguished from a pure |Angle| by virtue
of being bounded so that::

  -90.0 * u.deg <= angle(s) <= +90.0 * u.deg

Any attempt to set a value outside that range will result in a `ValueError`.
