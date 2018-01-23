# Licensed under a 3-clause BSD style license - see LICENSE.rst


import numpy as np

from matplotlib.transforms import Affine2D, IdentityTransform

from ....wcs import WCS

from ..transforms import WCSWorld2PixelTransform

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


def test_shorthand_inversion():
    """Test that the Matplotlib subtraction shorthand for composing and
    inverting transformations works."""
    w1 = WCS(naxis=2)
    w1.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    w1.wcs.crpix = [256.0, 256.0]
    w1.wcs.cdelt = [-0.05, 0.05]
    w1.wcs.crval = [120.0, -19.0]

    w2 = WCS(naxis=2)
    w2.wcs.ctype = ['RA---SIN', 'DEC--SIN']
    w2.wcs.crpix = [256.0, 256.0]
    w2.wcs.cdelt = [-0.05, 0.05]
    w2.wcs.crval = [235.0, +23.7]

    t1 = WCSWorld2PixelTransform(w1)
    t2 = WCSWorld2PixelTransform(w2)

    assert t1 - t2 == t1 + t2.inverted()
    assert t1 - t2 != t2.inverted() + t1
    assert t1 - t1 == IdentityTransform()


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

    w1 = WCSWorld2PixelTransform(WCS3D, slice=('y', 0, 'x')) + Affine2D()
    pixel = w1.transform(world)
    world_2 = w1.inverted().transform(pixel)

    np.testing.assert_allclose(world[:, 0], world_2[:, 0])
    np.testing.assert_allclose(world[:, 2], world_2[:, 2])
