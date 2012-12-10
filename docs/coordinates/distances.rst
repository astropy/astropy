Distances and Cartesian Representations
---------------------------------------

Coordinates can also have line-of-sight distances.  If these are provided, a
coordinate object becomes a full-fledged point in three-dimensional space.  If
not (i.e., the `distance` attribute of the coordinate object is `None`), the
point is interpreted as lying on the (dimensionless) unit sphere.

The `~astropy.coordinates.distances.Distance` class is provided to represent a
line-of-sight distance for a coordinate.  It must include a length unit to be
valid.::

    >>> from astropy.coordinates import Distance
    >>> d = Distance(770)
    UnitsError: A unit must be provided for distance.
    >>> d = Distance(770, u.kpc)
    >>> c = ICRSCoordinates('00h42m44.3s +41d16m9s', distance=d)
    >>> c
    <ICRSCoordinates RA=10.68458 deg, Dec=41.26917 deg, Distance=7.7e+02 kpc>

If a distance is available, the coordinate can be converted into cartesian
coordinates using the `x`/`y`/`z` attributes::

    >>> c.x
    568.7128882165681
    >>> c.y
    107.3009359688103
    >>> c.z
    507.8899092486349

.. note::

    The location of the origin is different for different coordinate
    systems, but for common celestial coordinate systems it is often
    the Earth center (or for precision work, the Earth/Moon barycenter).

The cartesian coordinates can also be accessed via the
`~astropy.coordinates.distances.CartesianCoordinates` object, which has
additional capabilities like arithmetic operations::

    >>> cp = c.cartesian
    >>> cp
    <CartesianPoints (568.712888217, 107.300935969, 507.889909249) kpc>
    >>> cp.x
    568.7128882165681
    >>> cp.y
    107.3009359688103
    >>> cp.z
    507.8899092486349
    >>> cp.unit
    Unit("kpc")
    >>> cp + cp
    <CartesianPoints (1137.42577643, 214.601871938, 1015.7798185) kpc>
    >>> cp - cp
    <CartesianPoints (0.0, 0.0, 0.0) kpc>

This cartesian representation can also be used to create a new coordinate
object, either directly or through a `CartesianPoints` object::

    >>> ICRSCoordinates(x=568.7129, y=107.3009, z=507.8899, unit=u.kpc)
    <ICRSCoordinates RA=10.68458 deg, Dec=41.26917 deg, Distance=7.7e+02 kpc>
    >>> cp = CartesianPoints(x=568.7129, y=107.3009, z=507.8899, unit=u.kpc)
    >>> ICRSCoordinates(cp)
    <ICRSCoordinates RA=10.68458 deg, Dec=41.26917 deg, Distance=7.7e+02 kpc>

Finally, two coordinates with distances can be used to derive a real-space
distance (i.e., non-projected separation)::

    >>> c1 = ICRSCoordinates('5h23m34.5s -69d45m22s', distance=Distance(49, u.kpc))
    >>> c2 = ICRSCoordinates('0h52m44.8s -72d49m43s', distance=Distance(61, u.kpc))
    >>> sep3d = c1.separation_3d(c2)
    >>> sep3d
    <Distance 23.05685 kpc>
    >>> sep3d.kpc
    23.05684814695706
    >>> sep3d.Mpc
    0.02305684814695706
    >>> sep3d.au
    4755816315.663559