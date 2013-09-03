Separations
-----------

The on-sky separation is easily computed with the `separation` method, which
computes the great-circle distance (*not* the small-angle approximation)::

    >>> from astropy.coordinates import ICRSCoordinates
    >>> c1 = ICRSCoordinates('5h23m34.5s -69d45m22s')
    >>> c2 = ICRSCoordinates('0h52m44.8s -72d49m43s')
    >>> sep = c1.separation(c2)
    >>> sep
    <Angle 20d44m46.02638s>


The `~astropy.coordinates.angles.AngularSeparation` object is a subclass of
`~astropy.coordinates.angles.Angle`, so it can be accessed in the same ways,
along with a few additions::

    >>> sep.radian
    0.36208807374669766
    >>> sep.hour
    1.383074562513832
    >>> sep.arcminute
    1244.7671062624488
    >>> sep.arcsecond
    74686.02637574692
