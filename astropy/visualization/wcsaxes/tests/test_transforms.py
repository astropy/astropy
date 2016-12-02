import numpy as np
from astropy.wcs import WCS

from matplotlib.transforms import Affine2D

from ..transforms import (WCSWorld2PixelTransform, WCSPixel2WorldTransform,
                          CoordinateTransform)

WCS2D = WCS(naxis=2)
WCS2D.wcs.ctype = ['x', 'y']
WCS2D.wcs.cunit = ['km', 'km']
WCS2D.wcs.crpix = [614.5, 856.5]
WCS2D.wcs.cdelt = [6.25, 6.25]
WCS2D.wcs.crval = [0., 0.]

WCS3D = WCS(naxis=3)
WCS3D.wcs.ctype = ['x', 'y', 'z']
WCS3D.wcs.cunit = ['km', 'km', 'km']
WCS3D.wcs.crpix = [614.5, 856.5, 333]
WCS3D.wcs.cdelt = [6.25, 6.25, 23]
WCS3D.wcs.crval = [0., 0., 1.]


# We add Affine2D to catch the fact that in Matplotlib, having a Composite
# transform can end up in more strict requirements for the dimensionality.


def test_2d():

    world = np.ones((10, 2))

    w1 = WCSWorld2PixelTransform(WCS2D) + Affine2D()
    pixel = w1.transform(world)
    world_2 = w1.inverted().transform(pixel)

    np.testing.assert_allclose(world, world_2)


def test_3d():

    world = np.ones((10, 3))

    w1 = WCSWorld2PixelTransform(WCS3D, slice=('y',0,'x')) + Affine2D()
    pixel = w1.transform(world)
    world_2 = w1.inverted().transform(pixel)

    np.testing.assert_allclose(world[:,0], world_2[:,0])
    np.testing.assert_allclose(world[:,2], world_2[:,2])
