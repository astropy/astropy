.. _astropy-coordinates-fast-in-place:

Fast In-Place Modification of Coordinates
*****************************************

For some applications the recommended method of
:ref:`astropy-coordinates-modifying-in-place` may not be fast enough due to the
extensive validation performed in that process to ensure correctness.  Likewise,
you may find that creating another coordinate frame with different data using
`~astropy.coordinates.BaseCoordinateFrame.realize_frame` does not meet your
performance requirements.

For these high-performance situations, you can directly modify in-place the
representation data in the frame object as shown in this example::

    >>> import astropy.units as u
    >>> from astropy.coordinates import SkyCoord
    >>> sc = SkyCoord([1,2],[3,4], unit='deg')
    >>> sc.data.lon[()] = [10, 20] * u.deg
    >>> sc.data.lat[1] = 40 * u.deg

    >>> sc.cache.clear()  # IMPORTANT TO DO THIS!

    >>> sc  # doctest: +FLOAT_CMP
    <SkyCoord (ICRS): (ra, dec) in deg
        [(10., 3.), (20., 40.)]>

Notice that the ``.data`` representation object uses different names for the
components than in the coordinate object.  If you wish to inspect the
mapping between frame attributes (e.g., ``.ra``) and representation attributes
(e.g., ``.lon``) you can look at the following dictionary::

    >>> sc.representation_component_names
    {'ra': 'lon', 'dec': 'lat', 'distance': 'distance'}

.. warning::

   You *must* include the step to clear the cache as shown. Failing to do so
   will cause the object to be inconsistent and likely result in incorrect
   results. `~astropy.coordinates.SkyCoord`
   and `~astropy.coordinates.BaseCoordinateFrame` cache various kinds of
   information for performance reasons, so you need clear the cache so that
   the new representation values are used when required.

You should note that the only way to modify the data in a frame is by using
the ``.data`` attribute directly and not the aliases for components on the
frame.  For example the following will *appear* to give a correct
result but it does not actually modify the underlying representation data::

    >>> sc.ra[1] = 20 * u.deg  # THIS IS WRONG

This problem is related to the current implementation of performance-based
caching and cannot be easily resolved.
