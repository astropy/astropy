.. _example_quick_init_celestial:

Create a WCS from a center and field of view
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This example demonstrate how a WCS can be instantiated with a center coordinate
and field of view values:

>>> from astropy.coordinates import Angle, SkyCoord
>>> from astropy.wcs.utils import quick_init_celestial_wcs
>>> center = SkyCoord(266.4, -29, unit="deg")
>>> fov = Angle([0.4, 0.2], unit="deg")
>>> wcs = quick_init_celestial_wcs(center, fov,
...                                pixel_scale=Angle([-0.005, 0.005], unit="deg"))
