.. _example_skycoord_to_celestial_wcs:

Create a WCS from a center and physical parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This example demonstrates how a WCS can be instantiated with a center coordinate
and field of view values:

>>> from astropy.coordinates import Angle, SkyCoord
>>> from astropy.wcs.utils import skycoord_to_celestial_wcs
>>> center = SkyCoord(266.4, -29, unit="deg")
>>> fov = Angle([0.4, 0.2], unit="deg")
>>> wcs = skycoord_to_celestial_wcs(
...         center, fov=fov,
...         pixel_angular_size=Angle([-0.005, 0.005], unit="deg")
...         )

We can also provide the shape of the image in pixels for a finer control:

>>> from astropy.coordinates import Angle, SkyCoord
>>> from astropy.wcs.utils import skycoord_to_celestial_wcs
>>> center = SkyCoord(266.4, -29, unit="deg")
>>> shape = [80, 40]
>>> wcs = skycoord_to_celestial_wcs(
...         center, pixel_shape=shape,
...         pixel_angular_size=Angle([-0.005, 0.005], unit="deg")
...         )
