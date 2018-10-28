# Licensed under a 3-clause BSD style license - see LICENSE.rst
import os

import pytest
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle
from matplotlib import rc_context

from .... import units as u
from ....io import fits
from ....wcs import WCS
from ....coordinates import SkyCoord

from ..patches import SphericalCircle
from .. import WCSAxes
from . import datasets
from ....tests.image_tests import IMAGE_REFERENCE_DIR
from ..frame import EllipticalFrame


class BaseImageTests:

    @classmethod
    def setup_class(cls):

        cls._data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'data'))

        msx_header = os.path.join(cls._data_dir, 'msx_header')
        cls.msx_header = fits.Header.fromtextfile(msx_header)

        rosat_header = os.path.join(cls._data_dir, 'rosat_header')
        cls.rosat_header = fits.Header.fromtextfile(rosat_header)

        twoMASS_k_header = os.path.join(cls._data_dir, '2MASS_k_header')
        cls.twoMASS_k_header = fits.Header.fromtextfile(twoMASS_k_header)

        cube_header = os.path.join(cls._data_dir, 'cube_header')
        cls.cube_header = fits.Header.fromtextfile(cube_header)

        slice_header = os.path.join(cls._data_dir, 'slice_header')
        cls.slice_header = fits.Header.fromtextfile(slice_header)


class TestBasic(BaseImageTests):

    @pytest.mark.remote_data(source='astropy')
    @pytest.mark.mpl_image_compare(baseline_dir=IMAGE_REFERENCE_DIR,
                                   tolerance=0, style={})
    def test_image_plot(self):
        # Test for plotting image and also setting values of ticks
        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], projection=WCS(self.msx_header), aspect='equal')
        ax.set_xlim(-0.5, 148.5)
        ax.set_ylim(-0.5, 148.5)
        ax.coords[0].set_ticks([-0.30, 0., 0.20] * u.degree, size=5, width=1)
        return fig

    @pytest.mark.remote_data(source='astropy')
    @pytest.mark.mpl_image_compare(baseline_dir=IMAGE_REFERENCE_DIR,
                                   tolerance=1.5, style={})
    @pytest.mark.parametrize('axisbelow', [True, False, 'line'])
    def test_axisbelow(self, axisbelow):
        # Test that tick marks, labels, and gridlines are drawn with the
        # correct zorder controlled by the axisbelow property.
        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], projection=WCS(self.msx_header), aspect='equal')
        ax.set_axisbelow(axisbelow)
        ax.set_xlim(-0.5, 148.5)
        ax.set_ylim(-0.5, 148.5)
        ax.coords[0].set_ticks([-0.30, 0., 0.20] * u.degree, size=5, width=1)
        ax.grid()

        # Add an image (default zorder=0).
        ax.imshow(np.zeros((64, 64)))

        # Add a patch (default zorder=1).
        r = Rectangle((30., 50.), 60., 50., facecolor='green', edgecolor='red')
        ax.add_patch(r)

        # Add a line (default zorder=2).
        ax.plot([32, 128], [32, 128], linewidth=10)

        return fig

    @pytest.mark.remote_data(source='astropy')
    @pytest.mark.mpl_image_compare(baseline_dir=IMAGE_REFERENCE_DIR,
                                   tolerance=0, style={})
    def test_contour_overlay(self):
        # Test for overlaying contours on images
        hdu_msx = datasets.fetch_msx_hdu()
        wcs_msx = WCS(self.msx_header)

        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_axes([0.15, 0.15, 0.8, 0.8],
                          projection=WCS(self.twoMASS_k_header),
                          aspect='equal')
        ax.set_xlim(-0.5, 720.5)
        ax.set_ylim(-0.5, 720.5)

        # Overplot contour
        ax.contour(hdu_msx.data, transform=ax.get_transform(wcs_msx),
                   colors='orange', levels=[2.5e-5, 5e-5, 1.e-4])
        ax.coords[0].set_ticks(size=5, width=1)
        ax.coords[1].set_ticks(size=5, width=1)
        ax.set_xlim(0., 720.)
        ax.set_ylim(0., 720.)

        # In previous versions, all angle axes defaulted to being displayed in
        # degrees. We now automatically show RA axes in hour angle units, but
        # for backward-compatibility with previous reference images we
        # explicitly use degrees here.
        ax.coords[0].set_format_unit(u.degree)

        return fig

    @pytest.mark.remote_data(source='astropy')
    @pytest.mark.mpl_image_compare(baseline_dir=IMAGE_REFERENCE_DIR,
                                   tolerance=0, style={})
    def test_contourf_overlay(self):
        # Test for overlaying contours on images
        hdu_msx = datasets.fetch_msx_hdu()
        wcs_msx = WCS(self.msx_header)

        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_axes([0.15, 0.15, 0.8, 0.8],
                          projection=WCS(self.twoMASS_k_header),
                          aspect='equal')
        ax.set_xlim(-0.5, 720.5)
        ax.set_ylim(-0.5, 720.5)

        # Overplot contour
        ax.contourf(hdu_msx.data, transform=ax.get_transform(wcs_msx),
                    levels=[2.5e-5, 5e-5, 1.e-4])
        ax.coords[0].set_ticks(size=5, width=1)
        ax.coords[1].set_ticks(size=5, width=1)
        ax.set_xlim(0., 720.)
        ax.set_ylim(0., 720.)

        # In previous versions, all angle axes defaulted to being displayed in
        # degrees. We now automatically show RA axes in hour angle units, but
        # for backward-compatibility with previous reference images we
        # explicitly use degrees here.
        ax.coords[0].set_format_unit(u.degree)

        return fig

    @pytest.mark.remote_data(source='astropy')
    @pytest.mark.mpl_image_compare(baseline_dir=IMAGE_REFERENCE_DIR,
                                   tolerance=0, style={})
    def test_overlay_features_image(self):

        # Test for overlaying grid, changing format of ticks, setting spacing
        # and number of ticks

        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_axes([0.25, 0.25, 0.65, 0.65],
                          projection=WCS(self.msx_header), aspect='equal')

        # Change the format of the ticks
        ax.coords[0].set_major_formatter('dd:mm:ss')
        ax.coords[1].set_major_formatter('dd:mm:ss.ssss')

        # Overlay grid on image
        ax.grid(color='red', alpha=1.0, lw=1, linestyle='dashed')

        # Set the spacing of ticks on the 'glon' axis to 4 arcsec
        ax.coords['glon'].set_ticks(spacing=4 * u.arcsec, size=5, width=1)

        # Set the number of ticks on the 'glat' axis to 9
        ax.coords['glat'].set_ticks(number=9, size=5, width=1)

        # Set labels on axes
        ax.coords['glon'].set_axislabel('Galactic Longitude', minpad=1.6)
        ax.coords['glat'].set_axislabel('Galactic Latitude', minpad=-0.75)

        # Change the frame linewidth and color
        ax.coords.frame.set_color('red')
        ax.coords.frame.set_linewidth(2)

        assert ax.coords.frame.get_color() == 'red'
        assert ax.coords.frame.get_linewidth() == 2

        return fig

    @pytest.mark.remote_data(source='astropy')
    @pytest.mark.mpl_image_compare(baseline_dir=IMAGE_REFERENCE_DIR,
                                   tolerance=0, style={})
    def test_curvilinear_grid_patches_image(self):

        # Overlay curvilinear grid and patches on image

        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8],
                          projection=WCS(self.rosat_header), aspect='equal')

        ax.set_xlim(-0.5, 479.5)
        ax.set_ylim(-0.5, 239.5)

        ax.grid(color='black', alpha=1.0, lw=1, linestyle='dashed')

        p = Circle((300, 100), radius=40, ec='yellow', fc='none')
        ax.add_patch(p)

        p = Circle((30., 20.), radius=20., ec='orange', fc='none',
                   transform=ax.get_transform('world'))
        ax.add_patch(p)

        p = Circle((60., 50.), radius=20., ec='red', fc='none',
                   transform=ax.get_transform('fk5'))
        ax.add_patch(p)

        p = Circle((40., 60.), radius=20., ec='green', fc='none',
                   transform=ax.get_transform('galactic'))
        ax.add_patch(p)

        return fig

    @pytest.mark.remote_data(source='astropy')
    @pytest.mark.mpl_image_compare(baseline_dir=IMAGE_REFERENCE_DIR,
                                   tolerance=0, style={})
    def test_cube_slice_image(self):

        # Test for cube slicing

        fig = plt.figure()
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8],
                          projection=WCS(self.cube_header),
                          slices=(50, 'y', 'x'), aspect='equal')

        ax.set_xlim(-0.5, 52.5)
        ax.set_ylim(-0.5, 106.5)

        ax.coords[2].set_axislabel('Velocity m/s')

        ax.coords[1].set_ticks(spacing=0.2 * u.deg, width=1)
        ax.coords[2].set_ticks(spacing=400 * u.m / u.s, width=1)

        ax.coords[1].set_ticklabel(exclude_overlapping=True)
        ax.coords[2].set_ticklabel(exclude_overlapping=True)

        ax.coords[1].grid(grid_type='contours', color='red', linestyle='solid')
        ax.coords[2].grid(grid_type='contours', color='red', linestyle='solid')

        return fig

    @pytest.mark.remote_data(source='astropy')
    @pytest.mark.mpl_image_compare(baseline_dir=IMAGE_REFERENCE_DIR,
                                   tolerance=0, style={})
    def test_cube_slice_image_lonlat(self):

        # Test for cube slicing. Here we test with longitude and latitude since
        # there is some longitude-specific code in _update_grid_contour.

        fig = plt.figure()
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8],
                          projection=WCS(self.cube_header),
                          slices=('x', 'y', 50), aspect='equal')

        ax.set_xlim(-0.5, 106.5)
        ax.set_ylim(-0.5, 106.5)

        ax.coords[0].grid(grid_type='contours', color='blue', linestyle='solid')
        ax.coords[1].grid(grid_type='contours', color='red', linestyle='solid')

        # In previous versions, all angle axes defaulted to being displayed in
        # degrees. We now automatically show RA axes in hour angle units, but
        # for backward-compatibility with previous reference images we
        # explicitly use degrees here.
        ax.coords[0].set_format_unit(u.degree)

        return fig

    @pytest.mark.remote_data(source='astropy')
    @pytest.mark.mpl_image_compare(baseline_dir=IMAGE_REFERENCE_DIR,
                                   tolerance=0, style={})
    def test_plot_coord(self):
        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_axes([0.15, 0.15, 0.8, 0.8],
                          projection=WCS(self.twoMASS_k_header),
                          aspect='equal')
        ax.set_xlim(-0.5, 720.5)
        ax.set_ylim(-0.5, 720.5)

        c = SkyCoord(266 * u.deg, -29 * u.deg)
        ax.plot_coord(c, 'o')

        # In previous versions, all angle axes defaulted to being displayed in
        # degrees. We now automatically show RA axes in hour angle units, but
        # for backward-compatibility with previous reference images we
        # explicitly use degrees here.
        ax.coords[0].set_format_unit(u.degree)

        return fig

    @pytest.mark.remote_data(source='astropy')
    @pytest.mark.mpl_image_compare(baseline_dir=IMAGE_REFERENCE_DIR,
                                   tolerance=0, style={})
    def test_plot_line(self):
        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_axes([0.15, 0.15, 0.8, 0.8],
                          projection=WCS(self.twoMASS_k_header),
                          aspect='equal')
        ax.set_xlim(-0.5, 720.5)
        ax.set_ylim(-0.5, 720.5)

        c = SkyCoord([266, 266.8] * u.deg, [-29, -28.9] * u.deg)
        ax.plot_coord(c)

        # In previous versions, all angle axes defaulted to being displayed in
        # degrees. We now automatically show RA axes in hour angle units, but
        # for backward-compatibility with previous reference images we
        # explicitly use degrees here.
        ax.coords[0].set_format_unit(u.degree)

        return fig

    @pytest.mark.remote_data(source='astropy')
    @pytest.mark.mpl_image_compare(baseline_dir=IMAGE_REFERENCE_DIR,
                                   tolerance=0, style={})
    def test_changed_axis_units(self):
        # Test to see if changing the units of axis works
        fig = plt.figure()
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8],
                          projection=WCS(self.cube_header),
                          slices=(50, 'y', 'x'), aspect='equal')
        ax.set_xlim(-0.5, 52.5)
        ax.set_ylim(-0.5, 106.5)
        ax.coords[2].set_major_formatter('x.xx')
        ax.coords[2].set_format_unit(u.km / u.s)
        ax.coords[2].set_axislabel('Velocity km/s')
        ax.coords[1].set_ticks(width=1)
        ax.coords[2].set_ticks(width=1)
        ax.coords[1].set_ticklabel(exclude_overlapping=True)
        ax.coords[2].set_ticklabel(exclude_overlapping=True)

        return fig

    @pytest.mark.remote_data(source='astropy')
    @pytest.mark.mpl_image_compare(baseline_dir=IMAGE_REFERENCE_DIR,
                                   tolerance=0, style={})
    def test_minor_ticks(self):
        # Test for drawing minor ticks
        fig = plt.figure()
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8],
                          projection=WCS(self.cube_header),
                          slices=(50, 'y', 'x'), aspect='equal')
        ax.set_xlim(-0.5, 52.5)
        ax.set_ylim(-0.5, 106.5)
        ax.coords[2].set_ticklabel(exclude_overlapping=True)
        ax.coords[1].set_ticklabel(exclude_overlapping=True)
        ax.coords[2].display_minor_ticks(True)
        ax.coords[1].display_minor_ticks(True)
        ax.coords[2].set_minor_frequency(3)
        ax.coords[1].set_minor_frequency(10)

        return fig

    @pytest.mark.remote_data(source='astropy')
    @pytest.mark.mpl_image_compare(baseline_dir=IMAGE_REFERENCE_DIR,
                                   tolerance=0, style={})
    def test_ticks_labels(self):
        fig = plt.figure(figsize=(6, 6))
        ax = WCSAxes(fig, [0.1, 0.1, 0.7, 0.7], wcs=None)
        fig.add_axes(ax)
        ax.set_xlim(-0.5, 2)
        ax.set_ylim(-0.5, 2)
        ax.coords[0].set_ticks(size=10, color='blue', alpha=0.2, width=1)
        ax.coords[1].set_ticks(size=20, color='red', alpha=0.9, width=1)
        ax.coords[0].set_ticks_position('all')
        ax.coords[1].set_ticks_position('all')
        ax.coords[0].set_axislabel('X-axis', size=20)
        ax.coords[1].set_axislabel('Y-axis', color='green', size=25,
                                   weight='regular', style='normal',
                                   family='cmtt10')
        ax.coords[0].set_axislabel_position('t')
        ax.coords[1].set_axislabel_position('r')
        ax.coords[0].set_ticklabel(color='purple', size=15, alpha=1,
                                   weight='light', style='normal',
                                   family='cmss10')
        ax.coords[1].set_ticklabel(color='black', size=18, alpha=0.9,
                                   weight='bold', family='cmr10')
        ax.coords[0].set_ticklabel_position('all')
        ax.coords[1].set_ticklabel_position('r')

        return fig

    @pytest.mark.remote_data(source='astropy')
    @pytest.mark.mpl_image_compare(baseline_dir=IMAGE_REFERENCE_DIR,
                                   tolerance=0, style={})
    def test_rcparams(self):

        # Test custom rcParams

        with rc_context({

                'axes.labelcolor': 'purple',
                'axes.labelsize': 14,
                'axes.labelweight': 'bold',

                'axes.linewidth': 3,
                'axes.facecolor': '0.5',
                'axes.edgecolor': 'green',

                'xtick.color': 'red',
                'xtick.labelsize': 8,
                'xtick.direction': 'in',

                'xtick.minor.visible': True,
                'xtick.minor.size': 5,

                'xtick.major.size': 20,
                'xtick.major.width': 3,
                'xtick.major.pad': 10,

                'grid.color': 'blue',
                'grid.linestyle': ':',
                'grid.linewidth': 1,
                'grid.alpha': 0.5}):

            fig = plt.figure(figsize=(6, 6))
            ax = WCSAxes(fig, [0.15, 0.1, 0.7, 0.7], wcs=None)
            fig.add_axes(ax)
            ax.set_xlim(-0.5, 2)
            ax.set_ylim(-0.5, 2)
            ax.grid()
            ax.set_xlabel('X label')
            ax.set_ylabel('Y label')
            ax.coords[0].set_ticklabel(exclude_overlapping=True)
            ax.coords[1].set_ticklabel(exclude_overlapping=True)
            return fig

    @pytest.mark.remote_data(source='astropy')
    @pytest.mark.mpl_image_compare(baseline_dir=IMAGE_REFERENCE_DIR,
                                   tolerance=0, style={})
    def test_tick_angles(self):
        # Test that tick marks point in the correct direction, even when the
        # axes limits extend only over a few FITS pixels. Addresses #45, #46.
        w = WCS()
        w.wcs.ctype = ['RA---TAN', 'DEC--TAN']
        w.wcs.crval = [90, 70]
        w.wcs.cdelt = [16, 16]
        w.wcs.crpix = [1, 1]
        w.wcs.radesys = 'ICRS'
        w.wcs.equinox = 2000.0
        fig = plt.figure(figsize=(3, 3))
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], projection=w)
        ax.set_xlim(1, -1)
        ax.set_ylim(-1, 1)
        ax.grid(color='gray', alpha=0.5, linestyle='solid')
        ax.coords['ra'].set_ticks(color='red', size=20)
        ax.coords['dec'].set_ticks(color='red', size=20)
        # In previous versions, all angle axes defaulted to being displayed in
        # degrees. We now automatically show RA axes in hour angle units, but
        # for backward-compatibility with previous reference images we
        # explicitly use degrees here.
        ax.coords[0].set_format_unit(u.degree)
        return fig

    @pytest.mark.remote_data(source='astropy')
    @pytest.mark.mpl_image_compare(baseline_dir=IMAGE_REFERENCE_DIR,
                                   tolerance=0, style={})
    def test_tick_angles_non_square_axes(self):
        # Test that tick marks point in the correct direction, even when the
        # axes limits extend only over a few FITS pixels, and the axes are
        # non-square.
        w = WCS()
        w.wcs.ctype = ['RA---TAN', 'DEC--TAN']
        w.wcs.crval = [90, 70]
        w.wcs.cdelt = [16, 16]
        w.wcs.crpix = [1, 1]
        w.wcs.radesys = 'ICRS'
        w.wcs.equinox = 2000.0
        fig = plt.figure(figsize=(6, 3))
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], projection=w)
        ax.set_xlim(1, -1)
        ax.set_ylim(-1, 1)
        ax.grid(color='gray', alpha=0.5, linestyle='solid')
        ax.coords['ra'].set_ticks(color='red', size=20)
        ax.coords['dec'].set_ticks(color='red', size=20)
        # In previous versions, all angle axes defaulted to being displayed in
        # degrees. We now automatically show RA axes in hour angle units, but
        # for backward-compatibility with previous reference images we
        # explicitly use degrees here.
        ax.coords[0].set_format_unit(u.degree)
        return fig

    @pytest.mark.remote_data(source='astropy')
    @pytest.mark.mpl_image_compare(baseline_dir=IMAGE_REFERENCE_DIR,
                                   tolerance=0, style={})
    def test_set_coord_type(self):
        # Test for setting coord_type
        fig = plt.figure(figsize=(3, 3))
        ax = fig.add_axes([0.2, 0.2, 0.6, 0.6],
                          projection=WCS(self.msx_header),
                          aspect='equal')
        ax.set_xlim(-0.5, 148.5)
        ax.set_ylim(-0.5, 148.5)
        ax.coords[0].set_coord_type('scalar')
        ax.coords[1].set_coord_type('scalar')
        ax.coords[0].set_major_formatter('x.xxx')
        ax.coords[1].set_major_formatter('x.xxx')
        ax.coords[0].set_ticklabel(exclude_overlapping=True)
        ax.coords[1].set_ticklabel(exclude_overlapping=True)
        return fig

    @pytest.mark.remote_data(source='astropy')
    @pytest.mark.mpl_image_compare(baseline_dir=IMAGE_REFERENCE_DIR,
                                   tolerance=0, style={})
    def test_ticks_regression(self):
        # Regression test for a bug that caused ticks aligned exactly with a
        # sampled frame point to not appear. This also checks that tick labels
        # don't get added more than once, and that no error occurs when e.g.
        # the top part of the frame is all at the same coordinate as one of the
        # potential ticks (which causes the tick angle calculation to return
        # NaN).
        wcs = WCS(self.slice_header)
        fig = plt.figure(figsize=(3, 3))
        ax = fig.add_axes([0.25, 0.25, 0.5, 0.5],
                          projection=wcs, aspect='auto')
        limits = wcs.wcs_world2pix([0, 0], [35e3, 80e3], 0)[1]
        ax.set_ylim(*limits)
        ax.coords[0].set_ticks(spacing=0.002 * u.deg)
        ax.coords[1].set_ticks(spacing=5 * u.km / u.s)
        ax.coords[0].set_ticklabel(alpha=0.5)  # to see multiple labels
        ax.coords[1].set_ticklabel(alpha=0.5)
        ax.coords[0].set_ticklabel_position('all')
        ax.coords[1].set_ticklabel_position('all')
        return fig

    @pytest.mark.remote_data(source='astropy')
    @pytest.mark.mpl_image_compare(baseline_dir=IMAGE_REFERENCE_DIR,
                                   savefig_kwargs={'bbox_inches': 'tight'},
                                   tolerance=0, style={})
    def test_axislabels_regression(self):
        # Regression test for a bug that meant that if tick labels were made
        # invisible with ``set_visible(False)``, they were still added to the
        # list of bounding boxes for tick labels, but with default values of 0
        # to 1, which caused issues.
        wcs = WCS(self.msx_header)
        fig = plt.figure(figsize=(3, 3))
        ax = fig.add_axes([0.25, 0.25, 0.5, 0.5], projection=wcs, aspect='auto')
        ax.coords[0].set_axislabel("Label 1")
        ax.coords[1].set_axislabel("Label 2")
        ax.coords[1].set_axislabel_visibility_rule('always')
        ax.coords[1].ticklabels.set_visible(False)
        return fig

    @pytest.mark.remote_data(source='astropy')
    @pytest.mark.mpl_image_compare(baseline_dir=IMAGE_REFERENCE_DIR,
                                   savefig_kwargs={'bbox_inches': 'tight'},
                                   tolerance=0, style={})
    def test_noncelestial_angular(self, tmpdir):
        # Regression test for a bug that meant that when passing a WCS that had
        # angular axes and using set_coord_type to set the coordinates to
        # longitude/latitude, but where the WCS wasn't recognized as celestial,
        # the WCS units are not converted to deg, so we can't assume that
        # transform will always return degrees.

        wcs = WCS(naxis=2)

        wcs.wcs.ctype = ['solar-x', 'solar-y']
        wcs.wcs.cunit = ['arcsec', 'arcsec']

        fig = plt.figure(figsize=(3, 3))
        ax = fig.add_subplot(1, 1, 1, projection=wcs)

        ax.imshow(np.zeros([1024, 1024]), origin='lower')

        ax.coords[0].set_coord_type('longitude', coord_wrap=180)
        ax.coords[1].set_coord_type('latitude')

        ax.coords[0].set_major_formatter('s.s')
        ax.coords[1].set_major_formatter('s.s')

        ax.coords[0].set_format_unit(u.arcsec, show_decimal_unit=False)
        ax.coords[1].set_format_unit(u.arcsec, show_decimal_unit=False)

        ax.grid(color='white', ls='solid')

        # Force drawing (needed for format_coord)
        fig.savefig(tmpdir.join('nothing').strpath)

        assert ax.format_coord(512, 512) == '513.0 513.0 (world)'

        return fig

    @pytest.mark.remote_data(source='astropy')
    @pytest.mark.mpl_image_compare(baseline_dir=IMAGE_REFERENCE_DIR,
                                   savefig_kwargs={'bbox_inches': 'tight'},
                                   tolerance=0, style={})
    def test_patches_distortion(self, tmpdir):

        # Check how patches get distorted (and make sure that scatter markers
        # and SphericalCircle don't)

        wcs = WCS(self.msx_header)
        fig = plt.figure(figsize=(3, 3))
        ax = fig.add_axes([0.25, 0.25, 0.5, 0.5], projection=wcs, aspect='equal')

        # Pixel coordinates
        r = Rectangle((30., 50.), 60., 50., edgecolor='green', facecolor='none')
        ax.add_patch(r)

        # FK5 coordinates
        r = Rectangle((266.4, -28.9), 0.3, 0.3, edgecolor='cyan', facecolor='none',
                      transform=ax.get_transform('fk5'))
        ax.add_patch(r)

        # FK5 coordinates
        c = Circle((266.4, -29.1), 0.15, edgecolor='magenta', facecolor='none',
                   transform=ax.get_transform('fk5'))
        ax.add_patch(c)

        # Pixel coordinates
        ax.scatter([40, 100, 130], [30, 130, 60], s=100, edgecolor='red', facecolor=(1, 0, 0, 0.5))

        # World coordinates (should not be distorted)
        ax.scatter(266.78238, -28.769255, transform=ax.get_transform('fk5'), s=300,
                   edgecolor='red', facecolor='none')

        # World coordinates (should not be distorted)
        r = SphericalCircle((266.4 * u.deg, -29.1 * u.deg), 0.15 * u.degree,
                            edgecolor='purple', facecolor='none',
                            transform=ax.get_transform('fk5'))
        ax.add_patch(r)

        ax.coords[0].set_ticklabel_visible(False)
        ax.coords[1].set_ticklabel_visible(False)

        return fig

    @pytest.mark.remote_data(source='astropy')
    @pytest.mark.mpl_image_compare(baseline_dir=IMAGE_REFERENCE_DIR,
                                   tolerance=0, style={})
    def test_elliptical_frame(self):

        # Regression test for a bug (astropy/astropy#6063) that caused labels to
        # be incorrectly simplified.

        wcs = WCS(self.msx_header)
        fig = plt.figure(figsize=(5, 3))
        ax = fig.add_axes([0.2, 0.2, 0.6, 0.6], projection=wcs, frame_class=EllipticalFrame)
        return fig

    @pytest.mark.remote_data(source='astropy')
    @pytest.mark.mpl_image_compare(baseline_dir=IMAGE_REFERENCE_DIR,
                                   tolerance=0, style={})
    def test_hms_labels(self):
        # This tests the apparance of the hms superscripts in tick labels
        fig = plt.figure(figsize=(3, 3))
        ax = fig.add_axes([0.3, 0.2, 0.65, 0.6],
                          projection=WCS(self.twoMASS_k_header),
                          aspect='equal')
        ax.set_xlim(-0.5, 0.5)
        ax.set_ylim(-0.5, 0.5)
        ax.coords[0].set_ticks(spacing=0.2 * 15 * u.arcsec)
        return fig

    @pytest.mark.remote_data(source='astropy')
    @pytest.mark.mpl_image_compare(baseline_dir=IMAGE_REFERENCE_DIR,
                                   tolerance=0, style={'text.usetex': True})
    def test_latex_labels(self):
        fig = plt.figure(figsize=(3, 3))
        ax = fig.add_axes([0.3, 0.2, 0.65, 0.6],
                          projection=WCS(self.twoMASS_k_header),
                          aspect='equal')
        ax.set_xlim(-0.5, 0.5)
        ax.set_ylim(-0.5, 0.5)
        ax.coords[0].set_ticks(spacing=0.2 * 15 * u.arcsec)
        return fig

    @pytest.mark.remote_data(source='astropy')
    @pytest.mark.mpl_image_compare(baseline_dir=IMAGE_REFERENCE_DIR,
                                   tolerance=0, style={})
    def test_tick_params(self):

        # This is a test to make sure that tick_params works correctly. We try
        # and test as much as possible with a single reference image.

        wcs = WCS()
        wcs.wcs.ctype = ['lon', 'lat']

        fig = plt.figure(figsize=(6, 6))

        # The first subplot tests:
        # - that plt.tick_params works
        # - that by default both axes are changed
        # - changing the tick direction and appearance, the label appearance and padding
        ax = fig.add_subplot(2, 2, 1, projection=wcs)
        plt.tick_params(direction='in', length=20, width=5, pad=6, labelsize=6,
                        color='red', labelcolor='blue')

        # The second subplot tests:
        # - that specifying grid parameters doesn't actually cause the grid to
        #   be shown (as expected)
        # - that axis= can be given integer coordinates or their string name
        # - that the tick positioning works (bottom/left/top/right)
        # Make sure that we can pass things that can index coords
        ax = fig.add_subplot(2, 2, 2, projection=wcs)
        plt.tick_params(axis=0, direction='in', length=20, width=5, pad=4, labelsize=6,
                        color='red', labelcolor='blue', bottom=True, grid_color='purple')
        plt.tick_params(axis='lat', direction='out', labelsize=8,
                        color='blue', labelcolor='purple', left=True, right=True,
                        grid_color='red')

        # The third subplot tests:
        # - that ax.tick_params works
        # - that the grid has the correct settings once shown explicitly
        # - that we can use axis='x' and axis='y'
        ax = fig.add_subplot(2, 2, 3, projection=wcs)
        ax.tick_params(axis='x', direction='in', length=20, width=5, pad=20, labelsize=6,
                       color='red', labelcolor='blue', bottom=True,
                       grid_color='purple')
        ax.tick_params(axis='y', direction='out', labelsize=8,
                       color='blue', labelcolor='purple', left=True, right=True,
                       grid_color='red')
        plt.grid()

        # The final subplot tests:
        # - that we can use tick_params on a specific coordinate
        # - that the label positioning can be customized
        # - that the colors argument works
        # - that which='minor' works
        ax = fig.add_subplot(2, 2, 4, projection=wcs)
        ax.coords[0].tick_params(length=4, pad=2, colors='orange', labelbottom=True,
                                 labeltop=True, labelsize=10)
        ax.coords[1].display_minor_ticks(True)
        ax.coords[1].tick_params(which='minor', length=6)

        return fig
