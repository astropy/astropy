# Licensed under a 3-clause BSD style license - see LICENSE.rst
import os
import shutil
import tempfile
import matplotlib.pyplot as plt
from matplotlib.testing.compare import compare_images
from matplotlib.patches import Circle
from ..rc_utils import rc_context

from astropy import units as u
from astropy.io import fits
from astropy.tests.helper import pytest
from astropy.tests.helper import remote_data
from .. import datasets, WCS, WCSAxes


class BaseImageTests(object):

    @classmethod
    def setup_class(cls):
        cls._moduledir = os.path.dirname(__file__)
        cls._data_dir = os.path.abspath(os.path.join(cls._moduledir, 'data'))
        cls._baseline_images_dir = os.path.abspath(os.path.join(cls._moduledir, 'baseline_images'))

        cls._tolerance = 1.5

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

    # method to create baseline or test images
    def generate_or_test(self, generate, figure, image, bbox_inches=None):
        if generate is None:
            result_dir = tempfile.mkdtemp()
            test_image = os.path.abspath(os.path.join(result_dir, image))

            # distutils will put the baseline images in non-accessible places,
            # copy to our tmpdir to be sure to keep them in case of failure
            orig_baseline_image = os.path.abspath(os.path.join(self._baseline_images_dir, image))
            baseline_image = os.path.abspath(os.path.join(result_dir, 'baseline-'+image))
            shutil.copyfile(orig_baseline_image, baseline_image)

            figure.savefig(test_image, bbox_inches=bbox_inches)

            if not os.path.exists(baseline_image):
                raise Exception("""Image file not found for comparision test
                                Generated Image:
                                \t{test}
                                This is expected for new tests.""".format(
                                    test=test_image))

            msg = compare_images(baseline_image, test_image, tol=self._tolerance)

            if msg is None:
                shutil.rmtree(result_dir)
            else:
                raise Exception(msg)
        else:
            figure.savefig(os.path.abspath(os.path.join(generate, image)), bbox_inches=bbox_inches)
            pytest.skip("Skipping test, since generating data")

class TestBasic(BaseImageTests):

    # Test for plotting image and also setting values of ticks
    def test_image_plot(self, generate):
        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], projection=WCS(self.msx_header), aspect='equal')
        ax.set_xlim(-0.5, 148.5)
        ax.set_ylim(-0.5, 148.5)
        ax.coords[0].set_ticks([-0.30, 0., 0.20] * u.degree, size=5, width=1)

        self.generate_or_test(generate, fig, 'image_plot.png')

    # Test for overlaying contours on images
    @remote_data
    def test_contour_overlay(self, generate):

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

        self.generate_or_test(generate, fig, 'contour_overlay.png')

    # Test for overlaying grid, changing format of ticks, setting spacing and number of ticks
    def test_overlay_features_image(self, generate):

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

        self.generate_or_test(generate, fig, 'overlay_features_image.png')

    # Overlay curvilinear grid and patches on image
    def test_curvilinear_grid_patches_image(self, generate):

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

        self.generate_or_test(generate, fig, 'curvlinear_grid_patches_image.png')

    # Test for cube slicing
    def test_cube_slice_image(self, generate):

        fig = plt.figure()
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8],
                          projection=WCS(self.cube_header),
                          slices=(50, 'y', 'x'), aspect='equal')

        ax.set_xlim(-0.5, 52.5)
        ax.set_ylim(-0.5, 106.5)

        ax.coords[2].set_axislabel('Velocity m/s')

        ax.coords[1].set_ticks(spacing=0.2 * u.deg, width=1,
                               exclude_overlapping=True)
        ax.coords[2].set_ticks(spacing=400 * u.m / u.s, width=1,
                               exclude_overlapping=True)

        ax.coords[1].grid(grid_type='contours', color='red', linestyle='solid')
        ax.coords[2].grid(grid_type='contours', color='red', linestyle='solid')

        self.generate_or_test(generate, fig, 'cube_slice_image.png')

    # Test for cube slicing. Here we test with longitude and latitude since
    # there is some longitude-specific code in _update_grid_contour.
    def test_cube_slice_image_lonlat(self, generate):

        fig = plt.figure()
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8],
                          projection=WCS(self.cube_header),
                          slices=('x', 'y', 50), aspect='equal')

        ax.set_xlim(-0.5, 106.5)
        ax.set_ylim(-0.5, 106.5)

        ax.coords[0].grid(grid_type='contours', color='blue', linestyle='solid')
        ax.coords[1].grid(grid_type='contours', color='red', linestyle='solid')

        self.generate_or_test(generate, fig, 'cube_slice_image_lonlat.png')

    # Test to see if changing the units of axis works
    def test_changed_axis_units(self, generate):
        fig = plt.figure()
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8],
                          projection=WCS(self.cube_header),
                          slices=(50, 'y', 'x'), aspect='equal')
        ax.set_xlim(-0.5, 52.5)
        ax.set_ylim(-0.5, 106.5)
        ax.coords[2].set_major_formatter('x.xx')
        ax.coords[2].set_format_unit(u.km / u.s)
        ax.coords[2].set_axislabel('Velocity km/s')
        ax.coords[1].set_ticks(width=1, exclude_overlapping=True)
        ax.coords[2].set_ticks(width=1, exclude_overlapping=True)

        self.generate_or_test(generate, fig, 'changed_axis_units.png')

    # Test to for drawing minor ticks
    def test_minor_ticks(self, generate):
        fig = plt.figure()
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8],
                          projection=WCS(self.cube_header),
                          slices=(50, 'y', 'x'), aspect='equal')
        ax.set_xlim(-0.5, 52.5)
        ax.set_ylim(-0.5, 106.5)
        ax.coords[2].set_ticks(exclude_overlapping=True)
        ax.coords[1].set_ticks(exclude_overlapping=True)
        ax.coords[2].display_minor_ticks(True)
        ax.coords[1].display_minor_ticks(True)
        ax.coords[2].set_minor_frequency(3)
        ax.coords[1].set_minor_frequency(10)

        self.generate_or_test(generate, fig, 'minor_ticks_image.png')

    def test_ticks_labels(self, generate):
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
                                   family='monospace')
        ax.coords[0].set_axislabel_position('t')
        ax.coords[1].set_axislabel_position('r')
        ax.coords[0].set_ticklabel(color='purple', size=15, alpha=1,
                                   weight='light', style='normal',
                                   family='sans-serif')
        ax.coords[1].set_ticklabel(color='black', size=18, alpha=0.9,
                                   weight='bold', family='serif')
        ax.coords[0].set_ticklabel_position('all')
        ax.coords[1].set_ticklabel_position('r')

        self.generate_or_test(generate, fig, 'ticks_labels.png',
                              bbox_inches='tight')

    # Test default style (matplotlib.rcParams) for ticks and gridlines
    def test_rcparams(self, generate):
        with rc_context({
                'xtick.color': 'red',
                'xtick.major.size': 20,
                'xtick.major.width': 2,
                'grid.color': 'blue',
                'grid.linestyle': ':',
                'grid.linewidth': 1,
                'grid.alpha': 0.5}):
            fig = plt.figure(figsize=(6, 6))
            ax = WCSAxes(fig, [0.1, 0.1, 0.7, 0.7], wcs=None)
            fig.add_axes(ax)
            ax.set_xlim(-0.5, 2)
            ax.set_ylim(-0.5, 2)
            ax.grid()
            ax.coords[0].set_ticks(exclude_overlapping=True)
            ax.coords[1].set_ticks(exclude_overlapping=True)
            self.generate_or_test(generate, fig, 'rcparams.png')

    # Test that tick marks point in the correct direction, even when the
    # axes limits extend only over a few FITS pixels. Addresses #45, #46.
    def test_tick_angles(self, generate):
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
        self.generate_or_test(generate, fig, 'tick_angles.png')

    # Test that tick marks point in the correct direction, even when the
    # axes limits extend only over a few FITS pixels, and the axes are
    # non-square.
    def test_tick_angles_non_square_axes(self, generate):
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
        self.generate_or_test(generate, fig, 'tick_angles_non_square_axes.png')

    # Test for setting coord_type
    def test_set_coord_type(self, generate):
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
        ax.coords[0].set_ticks(exclude_overlapping=True)
        ax.coords[1].set_ticks(exclude_overlapping=True)

        self.generate_or_test(generate, fig, 'set_coord_type.png')

    def test_ticks_regression(self, generate):
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
        limits = wcs.wcs_world2pix([0, 0], [35e3, 80e3],0)[1]
        ax.set_ylim(*limits)
        ax.coords[0].set_ticks(spacing=0.002 * u.deg)
        ax.coords[1].set_ticks(spacing=5 * u.km / u.s)
        ax.coords[0].set_ticklabel(alpha=0.5)  # to see multiple labels
        ax.coords[1].set_ticklabel(alpha=0.5)
        ax.coords[0].set_ticklabel_position('all')
        ax.coords[1].set_ticklabel_position('all')
        self.generate_or_test(generate, fig, 'test_ticks_regression_1.png')

    def test_axislabels_regression(self, generate):
        # Regression test for a bug that meant that if tick labels were made
        # invisible with ``set_visible(False)``, they were still added to the
        # list of bounding boxes for tick labels, but with default values of 0
        # to 1, which caused issues.
        wcs = WCS(self.msx_header)
        fig = plt.figure(figsize=(3, 3))
        ax = fig.add_axes([0.25, 0.25, 0.5, 0.5], projection=wcs, aspect='auto')
        ax.coords[0].set_axislabel("Label 1")
        ax.coords[1].set_axislabel("Label 2")
        ax.coords[1].ticklabels.set_visible(False)
        self.generate_or_test(generate, fig, 'test_axislabels_regression.png',
                              bbox_inches='tight')
