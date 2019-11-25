# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
import warnings
from textwrap import dedent

import pytest
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.transforms import Affine2D, IdentityTransform

from astropy.io import fits
from astropy import units as u
from astropy.wcs.wcsapi import BaseLowLevelWCS, SlicedLowLevelWCS
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.units import Quantity
from astropy.tests.image_tests import IMAGE_REFERENCE_DIR
from astropy.wcs import WCS
from astropy.visualization.wcsaxes.frame import RectangularFrame, RectangularFrame1D
from astropy.visualization.wcsaxes.wcsapi import (WCSWorld2PixelTransform,
                                                  transform_coord_meta_from_wcs,
                                                  apply_slices)

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


@pytest.fixture
def wcs_4d():
    header = dedent("""\
    WCSAXES =                    4 / Number of coordinate axes
    CRPIX1  =                  0.0 / Pixel coordinate of reference point
    CRPIX2  =                  0.0 / Pixel coordinate of reference point
    CRPIX3  =                  0.0 / Pixel coordinate of reference point
    CRPIX4  =                  5.0 / Pixel coordinate of reference point
    CDELT1  =                  0.4 / [min] Coordinate increment at reference point
    CDELT2  =                2E-11 / [m] Coordinate increment at reference point
    CDELT3  =   0.0027777777777778 / [deg] Coordinate increment at reference point
    CDELT4  =   0.0013888888888889 / [deg] Coordinate increment at reference point
    CUNIT1  = 'min'                / Units of coordinate increment and value
    CUNIT2  = 'm'                  / Units of coordinate increment and value
    CUNIT3  = 'deg'                / Units of coordinate increment and value
    CUNIT4  = 'deg'                / Units of coordinate increment and value
    CTYPE1  = 'TIME'               / Coordinate type code
    CTYPE2  = 'WAVE'               / Vacuum wavelength (linear)
    CTYPE3  = 'HPLT-TAN'           / Coordinate type codegnomonic projection
    CTYPE4  = 'HPLN-TAN'           / Coordinate type codegnomonic projection
    CRVAL1  =                  0.0 / [min] Coordinate value at reference point
    CRVAL2  =                  0.0 / [m] Coordinate value at reference point
    CRVAL3  =                  0.0 / [deg] Coordinate value at reference point
    CRVAL4  =                  0.0 / [deg] Coordinate value at reference point
    LONPOLE =                180.0 / [deg] Native longitude of celestial pole
    LATPOLE =                  0.0 / [deg] Native latitude of celestial pole
    """)
    return WCS(header=fits.Header.fromstring(header, sep='\n'))


@pytest.fixture
def cube_wcs():
    data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'data'))
    cube_header = os.path.join(data_dir, 'cube_header')
    header = fits.Header.fromtextfile(cube_header)
    return WCS(header=header)


def test_shorthand_inversion():
    """
    Test that the Matplotlib subtraction shorthand for composing and inverting
    transformations works.
    """
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

    world = np.ones((10, 2))

    w1 = WCSWorld2PixelTransform(WCS3D[:, 0, :]) + Affine2D()
    pixel = w1.transform(world)
    world_2 = w1.inverted().transform(pixel)

    np.testing.assert_allclose(world[:, 0], world_2[:, 0])
    np.testing.assert_allclose(world[:, 1], world_2[:, 1])


def test_coord_type_from_ctype(cube_wcs):

    _, coord_meta = transform_coord_meta_from_wcs(cube_wcs, RectangularFrame,
                                                  slices=(50, 'y', 'x'))

    axislabel_position = coord_meta['default_axislabel_position']
    ticklabel_position = coord_meta['default_ticklabel_position']
    ticks_position = coord_meta['default_ticks_position']

    # These axes are swapped due to the pixel derivatives
    assert axislabel_position == ['l', 'r', 'b']
    assert ticklabel_position == ['l', 'r', 'b']
    assert ticks_position == ['l', 'r', 'b']

    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ['GLON-TAN', 'GLAT-TAN']
    wcs.wcs.crpix = [256.0] * 2
    wcs.wcs.cdelt = [-0.05] * 2
    wcs.wcs.crval = [50.0] * 2
    wcs.wcs.cname = ['Longitude', '']
    wcs.wcs.set()

    _, coord_meta = transform_coord_meta_from_wcs(wcs, RectangularFrame)

    assert coord_meta['type'] == ['longitude', 'latitude']
    assert coord_meta['format_unit'] == [u.deg, u.deg]
    assert coord_meta['wrap'] == [None, None]
    assert coord_meta['default_axis_label'] == ['Longitude', 'pos.galactic.lat']

    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ['HPLN-TAN', 'HPLT-TAN']
    wcs.wcs.crpix = [256.0] * 2
    wcs.wcs.cdelt = [-0.05] * 2
    wcs.wcs.crval = [50.0] * 2
    wcs.wcs.set()

    _, coord_meta = transform_coord_meta_from_wcs(wcs, RectangularFrame)

    assert coord_meta['type'] == ['longitude', 'latitude']
    assert coord_meta['format_unit'] == [u.arcsec, u.arcsec]
    assert coord_meta['wrap'] == [180., None]

    _, coord_meta = transform_coord_meta_from_wcs(wcs, RectangularFrame,
                                                  slices=('y', 'x'))

    axislabel_position = coord_meta['default_axislabel_position']
    ticklabel_position = coord_meta['default_ticklabel_position']
    ticks_position = coord_meta['default_ticks_position']

    # These axes should be swapped because of slices
    assert axislabel_position == ['l', 'b']
    assert ticklabel_position == ['l', 'b']
    assert ticks_position == ['bltr', 'bltr']

    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ['HGLN-TAN', 'HGLT-TAN']
    wcs.wcs.crpix = [256.0] * 2
    wcs.wcs.cdelt = [-0.05] * 2
    wcs.wcs.crval = [50.0] * 2
    wcs.wcs.set()

    _, coord_meta = transform_coord_meta_from_wcs(wcs, RectangularFrame)

    assert coord_meta['type'] == ['longitude', 'latitude']
    assert coord_meta['format_unit'] == [u.deg, u.deg]
    assert coord_meta['wrap'] == [180., None]

    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ['CRLN-TAN', 'CRLT-TAN']
    wcs.wcs.crpix = [256.0] * 2
    wcs.wcs.cdelt = [-0.05] * 2
    wcs.wcs.crval = [50.0] * 2
    wcs.wcs.set()

    _, coord_meta = transform_coord_meta_from_wcs(wcs, RectangularFrame)

    assert coord_meta['type'] == ['longitude', 'latitude']
    assert coord_meta['format_unit'] == [u.deg, u.deg]
    assert coord_meta['wrap'] == [360., None]

    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    wcs.wcs.crpix = [256.0] * 2
    wcs.wcs.cdelt = [-0.05] * 2
    wcs.wcs.crval = [50.0] * 2
    wcs.wcs.set()

    _, coord_meta = transform_coord_meta_from_wcs(wcs, RectangularFrame)

    assert coord_meta['type'] == ['longitude', 'latitude']
    assert coord_meta['format_unit'] == [u.hourangle, u.deg]
    assert coord_meta['wrap'] == [None, None]

    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ['spam', 'spam']
    wcs.wcs.crpix = [256.0] * 2
    wcs.wcs.cdelt = [-0.05] * 2
    wcs.wcs.crval = [50.0] * 2
    wcs.wcs.set()

    _, coord_meta = transform_coord_meta_from_wcs(wcs, RectangularFrame)

    assert coord_meta['type'] == ['scalar', 'scalar']
    assert coord_meta['format_unit'] == [u.one, u.one]
    assert coord_meta['wrap'] == [None, None]


def test_coord_type_1d_1d_wcs():
    wcs = WCS(naxis=1)
    wcs.wcs.ctype = ['WAVE']
    wcs.wcs.crpix = [256.0]
    wcs.wcs.cdelt = [-0.05]
    wcs.wcs.crval = [50.0]
    wcs.wcs.set()

    _, coord_meta = transform_coord_meta_from_wcs(wcs, RectangularFrame1D)

    assert coord_meta['type'] == ['scalar']
    assert coord_meta['format_unit'] == [u.m]
    assert coord_meta['wrap'] == [None]


def test_coord_type_1d_2d_wcs_correlated():
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ['GLON-TAN', 'GLAT-TAN']
    wcs.wcs.crpix = [256.0] * 2
    wcs.wcs.cdelt = [-0.05] * 2
    wcs.wcs.crval = [50.0] * 2
    wcs.wcs.set()

    _, coord_meta = transform_coord_meta_from_wcs(wcs, RectangularFrame1D, slices=('x', 0))

    assert coord_meta['type'] == ['longitude', 'latitude']
    assert coord_meta['format_unit'] == [u.deg, u.deg]
    assert coord_meta['wrap'] == [None, None]
    assert coord_meta['visible'] == [True, True]


def test_coord_type_1d_2d_wcs_uncorrelated():
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ['WAVE', 'UTC']
    wcs.wcs.crpix = [256.0] * 2
    wcs.wcs.cdelt = [-0.05] * 2
    wcs.wcs.crval = [50.0] * 2
    wcs.wcs.cunit = ['nm', 's']
    wcs.wcs.set()

    _, coord_meta = transform_coord_meta_from_wcs(wcs, RectangularFrame1D, slices=('x', 0))

    assert coord_meta['type'] == ['scalar', 'scalar']
    assert coord_meta['format_unit'] == [u.m, u.s]
    assert coord_meta['wrap'] == [None, None]
    assert coord_meta['visible'] == [True, False]


def test_coord_meta_4d(wcs_4d):

    _, coord_meta = transform_coord_meta_from_wcs(wcs_4d, RectangularFrame, slices=(0, 0, 'x', 'y'))

    axislabel_position = coord_meta['default_axislabel_position']
    ticklabel_position = coord_meta['default_ticklabel_position']
    ticks_position = coord_meta['default_ticks_position']

    assert axislabel_position == ['', '', 'b', 'l']
    assert ticklabel_position == ['', '', 'b', 'l']
    assert ticks_position == ['', '', 'bltr', 'bltr']


def test_coord_meta_4d_line_plot(wcs_4d):

    _, coord_meta = transform_coord_meta_from_wcs(wcs_4d, RectangularFrame1D, slices=(0, 0, 0, 'x'))

    axislabel_position = coord_meta['default_axislabel_position']
    ticklabel_position = coord_meta['default_ticklabel_position']
    ticks_position = coord_meta['default_ticks_position']

    # These axes are swapped due to the pixel derivatives
    assert axislabel_position == ['', '', 't', 'b']
    assert ticklabel_position == ['', '', 't', 'b']
    assert ticks_position == ['', '', 't', 'b']


@pytest.fixture
def sub_wcs(wcs_4d, wcs_slice):
    return SlicedLowLevelWCS(wcs_4d, wcs_slice)


@pytest.mark.parametrize(("wcs_slice", "wcsaxes_slices", "world_map", "ndim"),
                         [
                             (np.s_[...], [0,0,'x','y'], (2, 3), 2),
                             (np.s_[...], [0,'x',0,'y'], (1, 2, 3), 3),
                             (np.s_[...], ['x',0,0,'y'], (0, 2, 3), 3),
                             (np.s_[...], ['x','y',0,0], (0, 1), 2),
                             (np.s_[:,:,0,:], [0, 'x', 'y'], (1, 2), 2),
                             (np.s_[:,:,0,:], ['x', 0, 'y'], (0, 1, 2), 3),
                             (np.s_[:,:,0,:], ['x', 'y', 0], (0, 1, 2), 3),
                             (np.s_[:,0,:,:], ['x', 'y', 0], (0, 1), 2),
                         ])
def test_apply_slices(sub_wcs, wcs_slice, wcsaxes_slices, world_map, ndim):
    transform_wcs, _, out_world_map = apply_slices(sub_wcs, wcsaxes_slices)
    assert transform_wcs.world_n_dim == ndim

    assert out_world_map == world_map


# parametrize here to pass to the fixture
@pytest.mark.parametrize("wcs_slice", [np.s_[:,:,0,:]])
def test_sliced_ND_input(sub_wcs, wcs_slice):
    slices_wcsaxes = [0, 'x', 'y']

    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=FutureWarning)
        _, coord_meta = transform_coord_meta_from_wcs(sub_wcs, RectangularFrame, slices=slices_wcsaxes)

    assert all(len(x) == 3 for x in coord_meta.values())

    coord_meta['name'] = ['time', 'custom:pos.helioprojective.lat', 'custom:pos.helioprojective.lon']
    coord_meta['type'] = ['scalar', 'latitude', 'longitude']
    coord_meta['wrap'] = [None, None, 180.0]
    coord_meta['unit'] = [u.Unit("min"), u.Unit("deg"), u.Unit("deg")]
    coord_meta['visible'] = [False, True, True]
    coord_meta['format_unit'] = [u.Unit("min"), u.Unit("arcsec"), u.Unit("arcsec")]
    coord_meta['default_axislabel_position'] = ['', 'b', 't']
    coord_meta['default_ticklabel_position'] = ['', 'b', 't']
    coord_meta['default_ticks_position'] = ['', 'btlr', 'btlr']

    # Validate the axes initialize correctly
    plt.subplot(projection=sub_wcs, slices=slices_wcsaxes)
    plt.close('all')


class LowLevelWCS5D(BaseLowLevelWCS):

    @property
    def pixel_n_dim(self):
        return 2

    @property
    def world_n_dim(self):
        return 5

    @property
    def world_axis_physical_types(self):
        return ['em.freq', 'time', 'pos.eq.ra', 'pos.eq.dec', 'phys.polarization.stokes']

    @property
    def world_axis_units(self):
        return ['Hz', 'day', 'deg', 'deg', '']

    def pixel_to_world_values(self, *pixel_arrays):
        pixel_arrays = (list(pixel_arrays) * 3)[:-1]  # make list have 5 elements
        return [np.asarray(pix) * scale for pix, scale in zip(pixel_arrays, [10, 0.2, 0.4, 0.39, 2])]

    def array_index_to_world_values(self, *index_arrays):
        return self.pixel_to_world_values(index_arrays[::-1])[::-1]

    def world_to_pixel_values(self, *world_arrays):
        world_arrays = world_arrays[:2]  # make list have 2 elements
        return [np.asarray(world) / scale for world, scale in zip(world_arrays, [10, 0.2])]

    def world_to_array_index_values(self, *world_arrays):
        return np.round(self.world_to_array_index_values(world_arrays[::-1])[::-1]).astype(int)

    @property
    def world_axis_object_components(self):
        return [('freq', 0, 'value'),
                ('time', 0, 'mjd'),
                ('celestial', 0, 'spherical.lon.degree'),
                ('celestial', 1, 'spherical.lat.degree'),
                ('stokes', 0, 'value')]

    @property
    def world_axis_object_classes(self):
        return {'celestial': (SkyCoord, (), {'unit': 'deg'}),
                'time': (Time, (), {'format': 'mjd'}),
                'freq': (Quantity, (), {'unit': 'Hz'}),
                'stokes': (Quantity, (), {'unit': 'one'})}


class TestWCSAPI:

    def teardown_method(self, method):
        plt.close('all')

    @pytest.mark.remote_data(source='astropy')
    @pytest.mark.mpl_image_compare(baseline_dir=IMAGE_REFERENCE_DIR,
                                   tolerance=0, style={})
    def test_wcsapi_5d(self):
        # Test for plotting image and also setting values of ticks
        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], projection=LowLevelWCS5D())
        ax.set_xlim(-0.5, 148.5)
        ax.set_ylim(-0.5, 148.5)
        return fig
