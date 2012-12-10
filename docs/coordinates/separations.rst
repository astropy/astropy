Separations
-----------

The on-sky separation is easily computed with the `separation` method, which
computes the great-circle distance (*not* the small-angle approximation)::

    >>> c1 = ICRSCoordinates('5h23m34.5s -69d45m22s')
    >>> c2 = ICRSCoordinates('0h52m44.8s -72d49m43s')
    >>> sep = c1.separation(c2)
    >>> sep
    <AngularSeparation 20.74612 deg>

The `~astropy.coordinates.angles.AngularSeparation` object is a subclass of
`~astropy.coordinates.angles.Angle`, so it can be accessed in the same ways,
along with a few additions::

    >>> sep.radians
    0.36208807374669766
    >>> sep.hours
    1.383074562513832
    >>> sep.arcmins
    1244.7671062624488
    >>> sep.arcsecs
    74686.02637574692