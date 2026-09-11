.. _astropy-coordinates-imposing-attributes:

Imposing the value of Frame Attributes
**************************************

In certain circumstances you may want to override the value of one frame attribute with another.
Examples of when this might be useful include when you are reading a series of frames from some source and need to change an attribute for a transform or calculation, or when you are using a package such as ``reproject`` where you don't control the construction of the coordinate objects.

Take this example, where you have two coordinates where the obstime difference results in a very subtly different origin of the coordinate frame.
For this example we are constructing these coordinates explicitly::

    >>> import astropy.units as u
    >>> from astropy.coordinates import SkyCoord, impose_frame_attributes

    >>> coord1 = SkyCoord(10*u.deg, 20*u.deg, frame='hcrs', obstime='2026-01-01 00:00:00')
    >>> coord2 = SkyCoord(20*u.deg, 30*u.deg, frame='hcrs', obstime='2026-01-01 00:00:00.001')

When you try and compute the separation an error is raised because of the origin shift::

    >>> print(coord1.separation(coord2))  # doctest: +SKIP
    ...
    astropy.units.errors.UnitsError: The input HCRS coordinates do not have length units. This probably means you created coordinates with lat/lon but no distance.  Heliocentric<->ICRS transforms cannot function in this case because there is an origin shift.

Given the very small difference in observer location assuming that they are equal is likely valid to a high degree of precision::

    >>> with impose_frame_attributes(obstime=coord1.obstime):
    ...     print(coord1.separation(coord2))
    13d28m54.20360928s

The `~astropy.coordinates.impose_frame_attributes` decorator accepts keyword arguments which are any attribute on any frame class.
So in our previous example of `~astropy.coordinates.HCRS` the ``obstime`` attribute.
