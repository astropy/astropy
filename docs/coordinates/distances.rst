Distances and Cartesian Representations
---------------------------------------

Coordinates can also have line-of-sight distances.  If these are provided, a
coordinate object becomes a full-fledged point in three-dimensional space.

The `~astropy.coordinates.Distance` class is provided to represent a
line-of-sight distance for a coordinate.  It must include a length unit to be
valid.::

    >>> from astropy.coordinates import Distance, ICRS
    >>> from astropy import units as u
    >>> d = Distance(770.)
    Traceback (most recent call last):
    ...
    UnitsError: No unit was provided for Distance
    >>> d = Distance(770., u.kpc)
    >>> d
    <Distance 770.0 kpc>
    >>> c = ICRS('00h42m44.3s +41d16m9s', distance=d)
    >>> c
    <ICRS RA=10.68458 deg, Dec=41.26917 deg, Distance=7.7e+02 kpc>

Because `Distance` is a subclass of `~astropy.units.Quantity`, in general a
`~astropy.units.Quantity` with units of length may be provided and it will
automatically convert to a `Distance`::

    >>> ICRS('00h42m44.3s +41d16m9s', distance=770 *  u.kpc)
    <ICRS RA=10.68458 deg, Dec=41.26917 deg, Distance=7.7e+02 kpc>

If a `distance` is present, the coordinate can be converted into Cartesian
coordinates using the `x`/`y`/`z` attributes (which are
`~astropy.units.Quantity` objects)::

    >>> c.x
    <Quantity 568.71288821656... kpc>
    >>> c.y
    <Quantity 107.30093596881... kpc>
    >>> c.z
    <Quantity 507.88990924863... kpc>

If a `distance` is not present, the Cartesian coordinates are still
available, but the point is interpreted as lying on the (dimensionless)
unit sphere::

    >>> c2 = ICRS('00h42m44.3s +41d16m9s')
    >>> c2.x
    <Quantity 0.73858816651502...>
    >>> c2.y
    <Quantity 0.13935186489455...>
    >>> c2.z
    <Quantity 0.65959728473848...>


.. note::

    The location of the origin is different for different coordinate
    systems, but for common celestial coordinate systems it is often
    the Earth center (or for precision work, the Earth/Moon barycenter).

The Cartesian coordinates can also be accessed via the
`~astropy.coordinates.CartesianPoints` object, which has
additional capabilities like arithmetic operations::

    >>> cp = c.cartesian
    >>> cp
    <CartesianPoints [ 568.71288822, 107.30093597, 507.88990925] kpc>
    >>> cp.x
    <Quantity 568.71288821656... kpc>
    >>> cp.y
    <Quantity 107.30093596881... kpc>
    >>> cp.z
    <Quantity 507.88990924863... kpc>
    >>> cp.unit
    Unit("kpc")
    >>> cp + cp
    <CartesianPoints [ 1137.42577643,  214.60187194, 1015.7798185 ] kpc>
    >>> cp - cp
    <CartesianPoints [ 0., 0., 0.] kpc>

This Cartesian representation can also be used to create a new coordinate
object, either directly or through a `CartesianPoints` object::

    >>> from astropy.coordinates import CartesianPoints
    >>> ICRS(x=568.7129, y=107.3009, z=507.8899, unit=u.kpc)
    <ICRS RA=10.68458 deg, Dec=41.26917 deg, Distance=7.7e+02 kpc>
    >>> cp = CartesianPoints(x=568.7129, y=107.3009, z=507.8899, unit=u.kpc)
    >>> ICRS(cp)
    <ICRS RA=10.68458 deg, Dec=41.26917 deg, Distance=7.7e+02 kpc>

Finally, two coordinates with distances can be used to derive a real-space
distance (i.e., non-projected separation)::

    >>> c1 = ICRS('5h23m34.5s -69d45m22s', distance=Distance(49, u.kpc))
    >>> c2 = ICRS('0h52m44.8s -72d49m43s', distance=Distance(61, u.kpc))
    >>> sep3d = c1.separation_3d(c2)
    >>> sep3d
    <Distance 23.056848146957... kpc>
    >>> sep3d.kpc
    23.056848146957...
    >>> sep3d.Mpc
    0.023056848146957...
    >>> sep3d.au
    4755816315.663...
