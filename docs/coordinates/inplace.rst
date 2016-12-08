In Place Modification of Coordinates
====================================

Coordinates are generally considered to be immutable. If you want to create
another coordinate frame with different data you should use
`~astropy.coordinates.BaseCoordinateFrame.realize_frame`, this is the safest way
to change the data in a frame as it creates a new frame with that data.

Due to the fact that this creates a new frame it is relatively slow, and in some
situations it is required that the data is replaced inside a frame with minimal
changes to the frame for maximum speed. This modification can be done by
modifiying the values of the representation data as follows::

    >>> import astropy.units as u
    >>> from astropy.coordinates import SkyCoord
    >>> c = SkyCoord([1,2],[3,4], unit='deg')
    >>> c.data.lon[()] = [10, 20] * u.deg
    >>> c
    <SkyCoord (ICRS): (ra, dec) in deg
        [( 10.,  3.), ( 20.,  4.)]>

This changes the longitude values for the frame, however
`~astropy.coordinates.SkyCoord` and `~astropy.coordinates.BaseCoordinateFrame`
cache representations of the data held within the frame as they are calculated
to speed up future operations of different representations of the same data. If
you modify the data in the representation the frame contains it is important to
clear this cache so future calculations of new representations are recomputed
using the new data. This can be achieved by doing::

    >>> del c.cache['representation']


It should be noted that the only way to modify the data in a frame is by using
the ``.data`` attribute directly and not the aliases for components on the frame
i.e.::

    >>> c.ra[()] = 20 * u.deg

will not work as a different representation object is used when acessing the
aliased component names. If you wish to inspect the mapping between frame
attributes i.e. ``.ra`` and representation attributes i.e. ``.lon`` you can look
at the::

    >>> c.representation_component_names
    OrderedDict([('ra', 'lon'), ('dec', 'lat'), ('distance', 'distance')])

dictionary.
