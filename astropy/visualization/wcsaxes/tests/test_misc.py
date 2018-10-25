# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
import warnings

import pytest
import numpy as np
import matplotlib.pyplot as plt

from .... import units as u
from ....wcs import WCS
from ....io import fits
from ....coordinates import SkyCoord
from ....tests.helper import catch_warnings
from ....tests.image_tests import ignore_matplotlibrc

from ..core import WCSAxes
from ..utils import get_coord_meta
from ..transforms import CurvedTransform

DATA = os.path.abspath(os.path.join(os.path.dirname(__file__), 'data'))


@ignore_matplotlibrc
def test_grid_regression():
    # Regression test for a bug that meant that if the rc parameter
    # axes.grid was set to True, WCSAxes would crash upon initalization.
    plt.rc('axes', grid=True)
    fig = plt.figure(figsize=(3, 3))
    WCSAxes(fig, [0.1, 0.1, 0.8, 0.8])


@ignore_matplotlibrc
def test_format_coord_regression(tmpdir):
    # Regression test for a bug that meant that if format_coord was called by
    # Matplotlib before the axes were drawn, an error occurred.
    fig = plt.figure(figsize=(3, 3))
    ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8])
    fig.add_axes(ax)
    assert ax.format_coord(10, 10) == ""
    assert ax.coords[0].format_coord(10) == ""
    assert ax.coords[1].format_coord(10) == ""
    fig.savefig(tmpdir.join('nothing').strpath)
    assert ax.format_coord(10, 10) == "10.0 10.0 (world)"
    assert ax.coords[0].format_coord(10) == "10.0"
    assert ax.coords[1].format_coord(10) == "10.0"


TARGET_HEADER = fits.Header.fromstring("""
NAXIS   =                    2
NAXIS1  =                  200
NAXIS2  =                  100
CTYPE1  = 'RA---MOL'
CRPIX1  =                  500
CRVAL1  =                180.0
CDELT1  =                 -0.4
CUNIT1  = 'deg     '
CTYPE2  = 'DEC--MOL'
CRPIX2  =                  400
CRVAL2  =                  0.0
CDELT2  =                  0.4
CUNIT2  = 'deg     '
COORDSYS= 'icrs    '
""", sep='\n')


@ignore_matplotlibrc
def test_no_numpy_warnings(tmpdir):

    # Make sure that no warnings are raised if some pixels are outside WCS
    # (since this is normal)

    ax = plt.subplot(1, 1, 1, projection=WCS(TARGET_HEADER))
    ax.imshow(np.zeros((100, 200)))
    ax.coords.grid(color='white')

    with catch_warnings(RuntimeWarning) as ws:
        plt.savefig(tmpdir.join('test.png').strpath)

    # For debugging
    for w in ws:
        print(w)

    assert len(ws) == 0


@ignore_matplotlibrc
def test_invalid_frame_overlay():

    # Make sure a nice error is returned if a frame doesn't exist
    ax = plt.subplot(1, 1, 1, projection=WCS(TARGET_HEADER))
    with pytest.raises(ValueError) as exc:
        ax.get_coords_overlay('banana')
    assert exc.value.args[0] == 'Unknown frame: banana'

    with pytest.raises(ValueError) as exc:
        get_coord_meta('banana')
    assert exc.value.args[0] == 'Unknown frame: banana'


@ignore_matplotlibrc
def test_plot_coord_transform():

    twoMASS_k_header = os.path.join(DATA, '2MASS_k_header')
    twoMASS_k_header = fits.Header.fromtextfile(twoMASS_k_header)
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_axes([0.15, 0.15, 0.8, 0.8],
                      projection=WCS(twoMASS_k_header),
                      aspect='equal')
    ax.set_xlim(-0.5, 720.5)
    ax.set_ylim(-0.5, 720.5)

    c = SkyCoord(359.76045223*u.deg, 0.26876217*u.deg)
    with pytest.raises(TypeError):
        ax.plot_coord(c, 'o', transform=ax.get_transform('galactic'))


@ignore_matplotlibrc
def test_set_label_properties():

    # Regression test to make sure that arguments passed to
    # set_xlabel/set_ylabel are passed to the underlying coordinate helpers

    ax = plt.subplot(1, 1, 1, projection=WCS(TARGET_HEADER))

    ax.set_xlabel('Test x label', labelpad=2, color='red')
    ax.set_ylabel('Test y label', labelpad=3, color='green')

    assert ax.coords[0].axislabels.get_text() == 'Test x label'
    assert ax.coords[0].axislabels.get_minpad('b') == 2
    assert ax.coords[0].axislabels.get_color() == 'red'

    assert ax.coords[1].axislabels.get_text() == 'Test y label'
    assert ax.coords[1].axislabels.get_minpad('l') == 3
    assert ax.coords[1].axislabels.get_color() == 'green'


GAL_HEADER = fits.Header.fromstring("""
SIMPLE  =                    T / conforms to FITS standard
BITPIX  =                  -32 / array data type
NAXIS   =                    3 / number of array dimensions
NAXIS1  =                   31
NAXIS2  =                 2881
NAXIS3  =                  480
EXTEND  =                    T
CTYPE1  = 'DISTMOD '
CRVAL1  =                  3.5
CDELT1  =                  0.5
CRPIX1  =                  1.0
CTYPE2  = 'GLON-CAR'
CRVAL2  =                180.0
CDELT2  =               -0.125
CRPIX2  =                  1.0
CTYPE3  = 'GLAT-CAR'
CRVAL3  =                  0.0
CDELT3  =                0.125
CRPIX3  =                241.0
""", sep='\n')


@ignore_matplotlibrc
def test_slicing_warnings(tmpdir):

    # Regression test to make sure that no warnings are emitted by the tick
    # locator for the sliced axis when slicing a cube.

    # Scalar case

    wcs3d = WCS(naxis=3)
    wcs3d.wcs.ctype = ['x', 'y', 'z']
    wcs3d.wcs.cunit = ['deg', 'deg', 'km/s']
    wcs3d.wcs.crpix = [614.5, 856.5, 333]
    wcs3d.wcs.cdelt = [6.25, 6.25, 23]
    wcs3d.wcs.crval = [0., 0., 1.]

    with warnings.catch_warnings(record=True) as warning_lines:
        warnings.resetwarnings()
        plt.subplot(1, 1, 1, projection=wcs3d, slices=('x', 'y', 1))
        plt.savefig(tmpdir.join('test.png').strpath)

    # For easy debugging if there are indeed warnings
    for warning in warning_lines:
        print(warning)

    assert len(warning_lines) == 0

    # Angle case

    wcs3d = WCS(GAL_HEADER)

    with warnings.catch_warnings(record=True) as warning_lines:
        warnings.resetwarnings()
        plt.subplot(1, 1, 1, projection=wcs3d, slices=('x', 'y', 2))
        plt.savefig(tmpdir.join('test.png').strpath)

    # For easy debugging if there are indeed warnings
    for warning in warning_lines:
        print(warning)

    assert len(warning_lines) == 0


def test_plt_xlabel_ylabel(tmpdir):

    # Regression test for a bug that happened when using plt.xlabel
    # and plt.ylabel with Matplotlib 3.0

    plt.subplot(projection=WCS())
    plt.xlabel('Galactic Longitude')
    plt.ylabel('Galactic Latitude')
    plt.savefig(tmpdir.join('test.png').strpath)


def test_grid_type_contours_transform(tmpdir):

    # Regression test for a bug that caused grid_type='contours' to not work
    # with custom transforms

    class CustomTransform(CurvedTransform):

        # We deliberately don't define the inverse, and has_inverse should
        # default to False.

        def transform(self, values):
            return values * 1.3

    transform = CustomTransform()
    coord_meta = {'type': ('scalar', 'scalar'),
                  'unit': (u.m, u.s),
                  'wrap': (None, None),
                  'name': ('x', 'y')}

    fig = plt.figure()
    ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8],
                 transform=transform, coord_meta=coord_meta)
    fig.add_axes(ax)
    ax.grid(grid_type='contours')
    fig.savefig(tmpdir.join('test.png').strpath)


def test_plt_imshow_origin():

    # Regression test for a bug that caused origin to be set to upper when
    # plt.imshow was called.

    ax = plt.subplot(projection=WCS())
    plt.imshow(np.ones((2, 2)))
    assert ax.get_xlim() == (-0.5, 1.5)
    assert ax.get_ylim() == (-0.5, 1.5)


def test_ax_imshow_origin():

    # Regression test for a bug that caused origin to be set to upper when
    # ax.imshow was called with no origin

    ax = plt.subplot(projection=WCS())
    ax.imshow(np.ones((2, 2)))
    assert ax.get_xlim() == (-0.5, 1.5)
    assert ax.get_ylim() == (-0.5, 1.5)


def test_grid_contour_large_spacing(tmpdir):

    # Regression test for a bug that caused a crash when grid was called and
    # didn't produce grid lines (due e.g. to too large spacing) and was then
    # called again.

    filename = tmpdir.join('test.png').strpath

    ax = plt.subplot(projection=WCS())
    ax.set_xlim(-0.5, 1.5)
    ax.set_ylim(-0.5, 1.5)
    ax.coords[0].set_ticks(values=[] * u.one)

    ax.coords[0].grid(grid_type='contours')
    plt.savefig(filename)

    ax.coords[0].grid(grid_type='contours')
    plt.savefig(filename)
