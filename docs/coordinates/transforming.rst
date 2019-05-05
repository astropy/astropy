.. include:: references.txt

.. _astropy-coordinates-transforming:

Transforming between Systems
****************************

`astropy.coordinates` supports a rich system for transforming
coordinates from one frame to another. While common astronomy frames
are built into Astropy, the transformation infrastructure is dynamic.
This means it allows users to define new coordinate frames and their
transformations. The topic of writing your own coordinate frame or
transforms is detailed in :ref:`astropy-coordinates-design`, and this
section is focused on how to *use* transformations.

The full list of built-in coordinate frames, the included transformations,
and the frame names are shown as a (clickable) graph in the
`~astropy.coordinates` API documentation.

The recommended method of transformation is shown below::

    >>> import astropy.units as u
    >>> from astropy.coordinates import SkyCoord
    >>> gc = SkyCoord(l=0*u.degree, b=45*u.degree, frame='galactic')
    >>> gc.fk5  # doctest: +FLOAT_CMP
    <SkyCoord (FK5: equinox=J2000.000): (ra, dec) in deg
        ( 229.27251463, -1.12844288)>

While this appears to be ordinary attribute-style access, it is actually
syntactic sugar for the more general
:meth:`~astropy.coordinates.SkyCoord.transform_to` method, which can
accept either a frame name, class, or instance::

    >>> from astropy.coordinates import FK5
    >>> gc.transform_to('fk5')  # doctest: +FLOAT_CMP
    <SkyCoord (FK5: equinox=J2000.000): (ra, dec) in deg
        ( 229.27251463, -1.12844288)>
    >>> gc.transform_to(FK5)  # doctest: +FLOAT_CMP
    <SkyCoord (FK5: equinox=J2000.000): (ra, dec) in deg
        ( 229.27251463, -1.12844288)>
    >>> gc.transform_to(FK5(equinox='J1980.0'))  # doctest: +FLOAT_CMP
    <SkyCoord (FK5: equinox=J1980.000): (ra, dec) in deg
        ( 229.0146935, -1.05560349)>

As a convenience, it is also possible to use a |SkyCoord| object as the frame in
:meth:`~astropy.coordinates.SkyCoord.transform_to`. This allows for putting one
coordinate object into the frame of another::

    >>> sc = SkyCoord(ra=1.0, dec=2.0, unit='deg', frame=FK5, equinox='J1980.0')
    >>> gc.transform_to(sc)  # doctest: +FLOAT_CMP
    <SkyCoord (FK5: equinox=J1980.000): (ra, dec) in deg
        ( 229.0146935, -1.05560349)>

Additionally, some coordinate frames (including `~astropy.coordinates.FK5`,
`~astropy.coordinates.FK4`, and `~astropy.coordinates.FK4NoETerms`) support
"self transformations," meaning the *type* of frame does not change, but the
frame attributes do. Any example is precessing a coordinate from one equinox
to another in an equatorial frame. This is done by passing ``transform_to`` a
frame class with the relevant attributes, as shown below. Note that these
frames use a default equinox if you do not specify one::

    >>> fk5c = SkyCoord('02h31m49.09s', '+89d15m50.8s', frame=FK5)
    >>> fk5c.equinox
    <Time object: scale='tt' format='jyear_str' value=J2000.000>
    >>> fk5c  # doctest: +FLOAT_CMP
    <SkyCoord (FK5: equinox=J2000.000): (ra, dec) in deg
        ( 37.95454167,  89.26411111)>
    >>> fk5_2005 = FK5(equinox='J2005')  # String initializes an astropy.time.Time object
    >>> fk5c.transform_to(fk5_2005)  # doctest: +FLOAT_CMP
    <SkyCoord (FK5: equinox=J2005.000): (ra, dec) in deg
        ( 39.39317639,  89.28584422)>

You can also specify the equinox when you create a coordinate using a
`~astropy.time.Time` object::

    >>> from astropy.time import Time
    >>> fk5c = SkyCoord('02h31m49.09s', '+89d15m50.8s',
    ...                 frame=FK5(equinox=Time('J1970')))
    >>> fk5_2000 = FK5(equinox=Time(2000, format='jyear'))
    >>> fk5c.transform_to(fk5_2000)  # doctest: +FLOAT_CMP
    <SkyCoord (FK5: equinox=2000.0): (ra, dec) in deg
        ( 48.023171,  89.38672485)>

The same lower-level frame classes also have a
:meth:`~astropy.coordinates.BaseCoordinateFrame.transform_to` method
that works the same as above, but they do not support attribute-style
access. They are also subtly different in that they only use frame
attributes present in the initial or final frame, while |skycoord|
objects use any frame attributes they have for all transformation
steps. So |skycoord| can always transform from one frame to another and
back again without change, while low-level classes may lose information
and hence often do not round-trip.

.. _astropy-coordinates-transforming-ephemerides:

Transformations and Solar-System Ephemerides
============================================

Some transformations (e.g., the transformation between
`~astropy.coordinates.ICRS` and `~astropy.coordinates.GCRS`) require the use of
a Solar-system ephemeris to calculate the position and velocity of the Earth
and Sun. By default, transformations are calculated using built-in
`ERFA <https://github.com/liberfa/erfa>`_ routines, but they can also use more
precise ones using the JPL ephemerides (which are derived from dynamical
models).

To use the JPL ephemerides, use the
`~astropy.coordinates.solar_system_ephemeris` context manager, as shown below:

.. doctest-requires:: jplephem

    >>> from astropy.coordinates import solar_system_ephemeris
    >>> from astropy.coordinates import GCRS
    >>> with solar_system_ephemeris.set('jpl'): # doctest: +REMOTE_DATA +IGNORE_OUTPUT
    ...     fk5c.transform_to(GCRS(obstime=Time("J2000"))) # doctest: +REMOTE_DATA +IGNORE_OUTPUT

For locations at large distances from the Solar system, using the JPL
ephemerides will make a negligible difference on the order of micro-arcseconds.
For nearby objects, such as the Moon, the difference can be of the
order of milli-arcseconds. For more details about what ephemerides
are available, including the requirements for using JPL ephemerides, see
:ref:`astropy-coordinates-solarsystem`.
