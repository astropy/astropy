Distances and Cartesian Representations
---------------------------------------

Coordinates can also have line-of-sight distances.  If these are provided, a
coordinate object becomes a full-fledged point in three-dimensional space.

The `~astropy.coordinates.Distance` class is provided to represent a
line-of-sight distance for a coordinate.  It must include a length unit to
be valid.::

    >>> from astropy.coordinates import Distance, ICRS, CartesianRepresentation
    >>> from astropy import units as u
    >>> d = Distance(770.)
    Traceback (most recent call last):
    ...
    UnitsError: No unit was provided for Distance
    >>> d = Distance(770., u.kpc)
    >>> d
    <Distance 770.0 kpc>
    >>> c = ICRS('00h42m44.3s', '+41d16m9s', distance=d)
    >>> c
    <ICRS Coordinate: ra=10.6845833333 deg, dec=41.2691666667 deg, distance=770.0 kpc>

Because `~astropy.coordinates.Distance` is a subclass of
`~astropy.units.Quantity`, in general a `~astropy.units.Quantity` with units
of length may be provided and it will automatically convert to a
`~astropy.coordinates.Distance`::

    >>> ICRS('00h42m44.3s', '+41d16m9s', distance=770*u.kpc)
    <ICRS Coordinate: ra=10.6845833333 deg, dec=41.2691666667 deg, distance=770.0 kpc>

By default, `~astropy.coordinates.Distance` values must be non-negative. For some cosmologies,
distance metrics may become zero or negative; negative distance values are supported
by setting the ``allow_negative`` keyword argument to ``True``::

    >>> d = Distance(-42., u.kpc)
    Traceback (most recent call last):
    ...
    ValueError: Distance must be >= 0. Set the kwarg 'allow_negative=True' to
    allow negative values.
    >>> d = Distance(-42., u.kpc, allow_negative=True)

If a ``distance`` is present, the coordinate can be converted into Cartesian
coordinates using the :attr:`~astropy.coordinates.CartesianPoints.x` /
:attr:`~astropy.coordinates.CartesianPoints.y` /
:attr:`~astropy.coordinates.CartesianPoints.z` attributes (which are
`~astropy.units.Quantity` objects)::

    >>> cart = c.represent_as(CartesianRepresentation)
    >>> cart.x
    <Quantity 568.71288821656... kpc>
    >>> cart.y
    <Quantity 107.30093596881... kpc>
    >>> cart.z
    <Quantity 507.88990924863... kpc>

If a ``distance`` is not present, the Cartesian coordinates are still
available, but the point is interpreted as lying on the (dimensionless)
unit sphere::

    >>> c2 = ICRS('00h42m44.3s', '+41d16m9s')
    >>> cart2 = c2.represent_as(CartesianRepresentation)
    >>> cart2.x
    <Quantity 0.73858816651502...>
    >>> cart2.y
    <Quantity 0.13935186489455...>
    >>> cart2.z
    <Quantity 0.65959728473848...>


.. note::

    The location of the origin is different for different coordinate
    systems, but for common celestial coordinate systems it is often
    the Earth center (or for precision work, the Earth/Moon barycenter).

This Cartesian representation can also be used to create a new coordinate
object through a `~astropy.coordinates.CartesianRepresentation` object::

    >>> cart = CartesianRepresentation(x=568.7129, y=107.3009, z=507.8899, unit=u.kpc)
    >>> ICRS(cart)
    <ICRS Coordinate: ra=10.6845796179 deg, dec=41.2691659084 deg, distance=769.99999759 kpc>

Finally, two coordinates with distances can be used to derive a real-space
distance (i.e., non-projected separation)::

    >>> c1 = ICRS('5h23m34.5s', '-69d45m22s', distance=49*u.kpc)
    >>> c2 = ICRS('0h52m44.8s', '-72d49m43s', distance=61*u.kpc)
    >>> sep3d = c1.separation_3d(c2)
    >>> sep3d
    <Distance 23.056848146957... kpc>
    >>> sep3d.kpc
    23.056848146957...
    >>> sep3d.Mpc
    0.023056848146957...
    >>> sep3d.au
    4755816315.663...
