.. include:: references.txt

.. _astropy-coordinates-transforming:

Transforming Between Systems
----------------------------

`astropy.coordinates` supports a rich system for transforming
coordinates from one system to another.  While common astronomy frames
are  built into Astropy, the transformation infrastructure is dynamic.
This means it allows users to define new coordinate frames and their
transformations.  The topic of writing your own coordinate frame or
transforms is detailed in :ref:`astropy-coordinates-design`, and this
section is focused on how to *use* transformations.

The simplest method of transformation is shown below::

    >>> import astropy.units as u
    >>> from astropy.coordinates import SkyCoord
    >>> gc = SkyCoord(l=0*u.degree, b=45*u.degree, frame='galactic')
    >>> gc.fk5  # doctest: +FLOAT_CMP
    <SkyCoord (FK5): equinox=J2000.000, ra=229.27250215 deg, dec=-1.12841764184 deg>

While this appears to be simple attribute-style access, it is actually
syntactic sugar for the more general
:meth:`~astropy.coordinates.SkyCoord.transform_to` method, which can
accept either a frame name, class or instance::

    >>> from astropy.coordinates import FK5
    >>> gc.transform_to('fk5')  # doctest: +FLOAT_CMP
    <SkyCoord (FK5): equinox=J2000.000, ra=229.27250215 deg, dec=-1.12841764184 deg>
    >>> gc.transform_to(FK5)  # doctest: +FLOAT_CMP
    <SkyCoord (FK5): equinox=J2000.000, ra=229.27250215 deg, dec=-1.12841764184 deg>
    >>> gc.transform_to(FK5(equinox='J1980.0'))  # doctest: +FLOAT_CMP
    <SkyCoord (FK5): equinox=J1980.000, ra=229.014681064 deg, dec=-1.05557823687 deg>

As a convenience it is also possible to use a |SkyCoord| object as the frame in
:meth:`~astropy.coordinates.SkyCoord.transform_to`.  This allows easily putting one
coordinate object into the frame of another::

    >>> sc = SkyCoord(ra=1.0, dec=2.0, unit='deg', frame=FK5, equinox='J1980.0')
    >>> gc.transform_to(sc)  # doctest: +FLOAT_CMP
    <SkyCoord (FK5): equinox=J1980.000, ra=229.014681064 deg, dec=-1.05557823687 deg>

The table below summarizes the built-in coordinate frames.  For details of
these frames and the transformations between them see the `astropy.coordinates`
API documentation and the `~astropy.coordinates.BaseCoordinateFrame` class
which forms the basis for all `astropy.coordinates` coordinate frames.

================================== ================
 Frame class                        Frame name
================================== ================
`~astropy.coordinates.ICRS`         ``icrs``
`~astropy.coordinates.FK5`          ``fk5``
`~astropy.coordinates.FK4`          ``fk4``
`~astropy.coordinates.FK4NoETerms`  ``fk4noeterms``
`~astropy.coordinates.Galactic`     ``galactic``
================================== ================

Additionally, some coordinate frames (including `~astropy.coordinates.FK5`,
`~astropy.coordinates.FK4`, and `~astropy.coordinates.FK4NoETerms`) support
"self transformations", meaning the *type* of frame doesn't change, but the
frame attributes do.  Any example is precessing a coordinate from one equinox
to another in an equatorial system. This is done by passing ``transform_to`` a
frame class with the relevant attributes, as shown below. Note that these
systems use a default equinox if you don't specify one::

    >>> fk5c = FK5('02h31m49.09s', '+89d15m50.8s')
    >>> fk5c.equinox
    <Time object: scale='utc' format='jyear_str' value=J2000.000>
    >>> fk5c  # doctest: +FLOAT_CMP
    <SkyCoord (FK5): equinox=J2000.000, ra=37.9545416667 deg, dec=89.2641111111 deg>
    >>> fk5_2005 = FK5(equinox='J2005')  # String initializes an astropy.time.Time object
    >>> fk5c.transform_to(fk5_2005)  # doctest: +FLOAT_CMP
    <SkyCoord (FK5): equinox=J2005.000, ra=39.3931763878 deg, dec=89.2858442155 deg>

You can also specify the equinox when you create a coordinate using an
`~astropy.time.Time` object::

    >>> from astropy.time import Time
    >>> fk5c = FK5('02h31m49.09s', '+89d15m50.8s',
    ...            equinox=Time('J1970', scale='utc'))
    >>> fk5_2000 = FK5(equinox=Time(2000, format='jyear', scale='utc'))
    >>> fk5c.transform_to(fk5_2000)  # doctest: +FLOAT_CMP
    <SkyCoord (FK5): equinox=2000.0, ra=48.0231710002 deg, dec=89.386724854 deg>

The same lower-level frame classes also have a
:meth:`~astropy.coordinates.BaseCoordinateFrame.transform_to` method
that works the same as above, but they do not support attribute-style
access. They are also subtly different in that they only use frame
attributes present in the initial or final frame, while |skycoord|
objects use any frame attributes they have for all transformation
steps.  So |skycoord| can always transform from one frame to another and
back again without change, while low-level classes may lose information
and hence often do not round-trip.
