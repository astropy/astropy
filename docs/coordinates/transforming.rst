Transforming Between Systems
----------------------------

`astropy.coordinates` supports a rich system for transforming coordinates from
one system to another.  The key concept is that a registry of all the
transformations is used to determine which coordinates can convert to others.
When you ask for a transformation, the registry (or "transformation graph") is
searched for the shortest path from your starting coordinate to your target, and
it applies all of the transformations in that path in series.   This allows only
the simplest transformations to be defined, and the package will automatically
determine how to combine those transformations to get from one system to
another.

As described above, there are two ways of transforming coordinates.  Coordinates
that have an alias (created with
`~astropy.coordinates.transformations.coordinate_alias`) can be converted by
simply using attribute style access to any other coordinate system::

    >>> gc = GalacticCoordinates(l=0, b=45, unit=(u.degree, u.degree))
    >>> gc.fk5
    <FK5Coordinates RA=229.27250 deg, Dec=-1.12842 deg>
    >>> ic = ICRSCoordinates(ra=0, dec=45, unit=(u.degree, u.degree)))
    >>> ic.fk5
    <FK5Coordinates RA=0.00001 deg, Dec=45.00000 deg>

While this appears to be simple attribute-style access, it is actually just
syntactic sugar for the `transform_to` method::

    >>> from astropy.coordinates import FK5Coordinates
    >>> gc.transform_to(FK5Coordinates)
    <FK5Coordinates RA=229.27250 deg, Dec=-1.12842 deg>
    >>> ic.transform_to(FK5Coordinates)
    <FK5Coordinates RA=0.00001 deg, Dec=45.00000 deg>

The full list of supported coordinate systems and transformations is in the
`astropy.coordinates` API documentation below.

Additionally, some coordinate systems support precessing the coordinate to
produce a new coordinate in the same system but at a different equinox.  Note 
that these systems have a default equinox they start with if you don't specify 
one::

    >>> fk5c = FK5Coordinates('02h31m49.09s +89d15m50.8s')
    >>> fk5c.equinox
    <Time object: scale='utc' format='jyear_str' vals=J2000.000>
    >>> fk5c
    <FK5Coordinates RA=37.95454 deg, Dec=89.26411 deg>
    >>> fk5c.precess_to(Time(2100, format='jyear', scale='utc'))
    <FK5Coordinates RA=88.32396 deg, Dec=89.54057 deg>

You can also specify the equinox when you create a coordinate using an 
`astropy.time.Time` object::

    >>> from astropy.time import Time
    >>> fk5c = FK5Coordinates('02h31m49.09s +89d15m50.8s', equinox=Time('J1970', scale='utc'))
    <FK5Coordinates RA=37.95454 deg, Dec=89.26411 deg>
    >>> fk5c.precess_to(Time(2000, format='jyear', scale='utc'))
    <FK5Coordinates RA=48.02317 deg, Dec=89.38672 deg>

Coordinate systems do not necessarily all support an equinox nor precession, as it is a
meaningless action for coordinate systems that do not depend on a particular equinox.

Furthermore, coordinates typically have an `obstime` attribute, intended to record the
time of the observation.  Some systems (especially FK4) require this information due
to being non-inertial frames (i.e., they rotate over time due to motions of the
defining stars).