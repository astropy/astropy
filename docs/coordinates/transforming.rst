Transforming Between Systems
----------------------------

`astropy.coordinates` supports a rich system for transforming
coordinates from one system to another.  While common astronomy frames
are  built into Astropy, the transformation infrastructure is dynamic.
This means it allows users to define new coordinate frames and their
transformations.  The topic of writing your own coordinate frame or
transforms is detailed in :ref :`astropy-coordinates-design`, and this
section is focused on how to *use* transformations.

The simplest for of transformation is simply::


    >>> import astropy.units as u
    >>> from astropy.coordinates import SkyCoord
    >>> gc = SkyCoord(l=0*u.degree, b=45*u.degree, frame='galactic')
    >>> gc.fk5
    <SkyCoord (FK5): equinox=J2000.000, ra=229.27250... deg, dec=-1.12841764... deg>

While this appears to be simple attribute-style access, it is actually
just syntactic sugar for the
:meth:`~astropy.coordinates.SkyCoord.transform_to` method, which can
accept either frame names, or `~astropy.coordinates.BaseCoordinateFrame`
classes::

    >>> from astropy.coordinates import FK5
    >>> gc.transform_to('fk5')
    <SkyCoord (FK5): equinox=J2000.000, ra=229.27250... deg, dec=-1.12841764... deg>
    >>> gc.transform_to(FK5)
    <SkyCoord (FK5): equinox=J2000.000, ra=229.27250... deg, dec=-1.12841764... deg>


The full list of supported coordinate systems and transformations is
in the `astropy.coordinates` API documentation.

Additionally, some coordinate frames support "self transformations",
meaning the *type* of frame doesn't change, but the frame attributes do.
Any example is precessing a coordinate from one equinox to another in an
equatorial system. This is done by passing `transform_to` a frame class
with the relevant attributes, as shown below. Note that these systems
have a default equinox they start with if you don't specify one::

    >>> fk5c = FK5('02h31m49.09s', '+89d15m50.8s')
    >>> fk5c.equinox
    <Time object: scale='utc' format='jyear_str' value=J2000.000>
    >>> fk5c
    <SkyCoord (FK5): equinox=J2000.000, ra=37.9545416... deg, dec=89.2641... deg>
    >>> fk5_2005 = FK5(equinox='J2005')  # internally the string becomes an astropy.time.Time object
    >>> fk5c.transform_to(fk5_2005)
    <SkyCoord (FK5): equinox=J2005.000, ra=39.3931763... deg, dec=89.2858... deg>

You can also specify the equinox when you create a coordinate using an
`~astropy.time.Time` object::

    >>> from astropy.time import Time
    >>> fk5c = FK5('02h31m49.09s', '+89d15m50.8s',
    ...            equinox=Time('J1970', scale='utc'))
    >>> fk5_2000 = FK5(equinox=Time(2000, format='jyear', scale='utc'))
    >>> fk5c.transform_to(fk5_2000)
    <SkyCoord (FK5): equinox=2000.0, ra=48.0231710002 deg, dec=89.386724854 deg>

The same lower-level frame classes also have a
:meth:`~astropy.coordinates.BaseCoordinateFrame.transform_to` method
that works the same as above, but they do not support attribute-style
access. They are also subtly different in that they only use frame
attributes present in the initial or final frame, while |skycoord|
objects use  any frame attributes they have for all transformation
steps.  So |skycoord| can always transform from one frame to another and
back again without change, while low-level classes may lose information
and hence often do not round-trip.
