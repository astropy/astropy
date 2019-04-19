In-Place Modification of Coordinates
************************************

Coordinates are generally considered to be immutable. If you want to create
another coordinate frame with different data you should use
`~astropy.coordinates.BaseCoordinateFrame.realize_frame`. This is the safest
way to change the data in a frame as it creates a new frame with that data.
Creating a new frame can be relatively slow, however, particularly for scalar
coordinates. Hence some situations may require that data be changed in-place in
an already existing frame object. This modification can be done by
modifying the values of the representation data as follows::

    >>> import astropy.units as u
    >>> from astropy.coordinates import SkyCoord
    >>> c = SkyCoord([1,2],[3,4], unit='deg')
    >>> c.data.lon[()] = [10, 20] * u.deg
    >>> c  # doctest: +FLOAT_CMP
    <SkyCoord (ICRS): (ra, dec) in deg
        [(10., 3.), (20., 4.)]>


This changes the longitude values of the frame. Unfortunately, doing this alone
introduces problems: `~astropy.coordinates.SkyCoord` and
`~astropy.coordinates.BaseCoordinateFrame` cache various kinds of information to
speed up some repeated operations. So we need to tell the cache that it should
be cleared so that it can be recalculated from the new data. This can be
achieved by doing::

    >>> c.cache.clear()

It should be noted that the only way to modify the data in a frame is by using
the ``.data`` attribute directly and not the aliases for components on the
frame (*i.e.*, the following will not work)::

    >>> c.ra[()] = 20 * u.deg

This is because a different representation object is used when accessing the
aliased component names. If you wish to inspect the mapping between frame
attributes (e.g., ``.ra`` and representation attributes ``.lon``) you can look
at the following dictionary::

    >>> c.representation_component_names
    OrderedDict([('ra', 'lon'), ('dec', 'lat'), ('distance', 'distance')])
