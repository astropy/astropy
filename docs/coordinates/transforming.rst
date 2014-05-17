Transforming Between Systems
----------------------------

`astropy.coordinates` supports a rich system for transforming
coordinates from one system to another.  While common astronomy frames
are  built into Astropy, the transformation infrastructure is dynamic.
This means it allows users to define new coordinate frames and their
transformations.  The topic of writing your own coordinate frame or
transforms is detailed in :ref :`astropy-coordinates-design`, and this
section is focused on how to *use* transformations.

.. todo:: update for APE5


    >>> import astropy.units as u
    >>> from astropy.coordinates import SkyCoord
    >>> gc = SkyCoord(l=0*u.degree, b=45*u.degree, frame='galactic')
    >>> gc.transform_to(FK5)
    <FK5 Coordinate: equinox=J2000.000, ra=229.27250215 deg, dec=-1.12841764184 deg>
    >>> ic = ICRS(ra=0*u.degree, dec=45*u.degree)
    >>> ic.transform_to(FK5)
    <FK5 Coordinate: equinox=J2000.000, ra=1.18888896168e-05 deg, dec=45.0000025278 deg>

While this appears to be simple attribute-style access, it is actually just
syntactic sugar for the
:meth:`~astropy.coordinates.Base.transform_to` method::

.. >>> from astropy.coordinates import FK5
.. >>> gc.transform_to(FK5)
.. <FK5 RA=229.27250 deg, Dec=-1.12842 deg>
.. >>> ic.transform_to(FK5)
.. <FK5 RA=0.00001 deg, Dec=45.00000 deg>

The full list of supported coordinate systems and transformations is
in the `astropy.coordinates` API documentation below.

Additionally, some coordinate systems support precessing the
coordinate to produce a new coordinate in the same system but at a
different equinox.  Note that these systems have a default equinox
they start with if you don't specify one::

    >>> from astropy.time import Time
    >>> fk5c = FK5('02h31m49.09s', '+89d15m50.8s')
    >>> fk5c.equinox
    <Time object: scale='utc' format='jyear_str' value=J2000.000>
    >>> fk5c
    <SkyCoord (FK5): equinox=J2000.000, ra=37.9545416667 deg, dec=89.2641111111 deg>
    >>> fk5_2100 = FK5(equinox=Time(2100, format='jyear', scale='utc'))
    >>> fk5c.transform_to(fk5_2100)
    <SkyCoord (FK5): equinox=2100.0, ra=88.3239593382 deg, dec=89.5405669962 deg>

You can also specify the equinox when you create a coordinate using a
`~astropy.time.Time` object::

    >>> fk5c = FK5('02h31m49.09s', '+89d15m50.8s',
    ...            equinox=Time('J1970', scale='utc'))
    >>> fk5_2000 = FK5(equinox=Time(2000, format='jyear', scale='utc'))
    >>> fk5c.transform_to(fk5_2000)
    <SkyCoord (FK5): equinox=2000.0, ra=48.0231710002 deg, dec=89.386724854 deg>

Coordinate systems do not necessarily all support an equinox nor
precession, as it is a meaningless action for coordinate systems that
do not depend on a particular equinox.

Furthermore, coordinates typically have an ``obstime`` attribute,
intended to record the time of the observation.  Some systems
(especially FK4) require this information due to being non-inertial
frames (i.e., they rotate over time due to motions of the defining
stars).
