Separations
-----------

The on-sky separation is easily computed with the `separation` method, which
computes the great-circle distance (*not* the small-angle approximation)::

    >>> from astropy.coordinates import ICRS
    >>> c1 = ICRS('5h23m34.5s -69d45m22s')
    >>> c2 = ICRS('0h52m44.8s -72d49m43s')
    >>> sep = c1.separation(c2)
    >>> sep
    <Angle 20.746118437707477 deg>


The returned object is an `~astropy.coordinates.angles.Angle` instance, so it
is straightforward to access the angle in any of several equivalent angular
units::

    >>> sep.radian
    0.3620880737466976
    >>> sep.hour
    1.383074562513832
    >>> sep.arcminute
    1244.7671062624488
    >>> sep.arcsecond
    74686.02637574692
