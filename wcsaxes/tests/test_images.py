import pytest
import os
import numpy as np
from astropy import units as u
import matplotlib.pyplot as plt
from matplotlib.testing.compare import compare_images
from matplotlib.patches import Circle
from astropy.wcs import WCS
from astropy.io import fits
from wcsaxes import WCSAxes
from matplotlib import cbook
from astropy.tests.helper import pytest
from astropy.utils.data import download_file
from astropy.tests.helper import remote_data
import wcsaxes.datasets as datasets


class TestImages(object):

    @classmethod
    def setup_class(cls):
        cls._filedir = os.path.abspath(__file__)
        cls._basedir = os.path.split(cls._filedir)[0]
        cls._baseline_images_dir = os.path.join(cls._basedir, 'baseline_images')
        cls._result_dir = os.path.abspath(os.path.join(cls._basedir, 'test_result_images'))

        if not os.path.exists(cls._result_dir):
            cbook.mkdirs(cls._result_dir)

        if not os.path.exists(cls._baseline_images_dir):
            cbook.mkdirs(cls._baseline_images_dir)

        cls._tolerance = 1

        cls._image1 = datasets.msx()
        cls._image2 = datasets.rosat()
        cls._image3 = datasets.twoMASS_k()
        cls._data_cube = datasets.l1448_co()

    # method to create baseline or test images
    def generate_or_test(self, generate, figure, image, test_image=None, baseline_image=None):
        baseline_image = os.path.abspath(os.path.join(self._baseline_images_dir, image))
        test_image = os.path.abspath(os.path.join(self._result_dir, image))
        if generate:
            figure.savefig(baseline_image)
            pytest.skip("Skipping test, since generating data")
        else:
            figure.savefig(test_image)
            msg = compare_images(baseline_image, test_image, tol=self._tolerance)
            assert msg is None

    # Test for plotting image and also setting values of ticks
    @remote_data
    def test_image_plot(self, generate):
        hdu = fits.open(self._image1)[0]
        fig = plt.figure(figsize=(6, 6))
        ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8], wcs=WCS(hdu.header))
        fig.add_axes(ax)
        ax.imshow(hdu.data, vmin=-1e-5, vmax=1e-4, origin='lower')
        ax.coords[0].set_ticks([-0.30, 0., 0.20] * u.degree)

        self.generate_or_test(generate, fig, 'image_plot.png')

    # Test for overlaying contours on images
    @remote_data
    def test_contour_overlay(self, generate):
        hdu = fits.open(self._image3)[0]

        hdu_msx = fits.open(self._image1)[0]
        wcs_msx = WCS(hdu_msx.header)

        fig = plt.figure(figsize=(6, 6))
        ax = WCSAxes(fig, [0.15, 0.15, 0.8, 0.8], wcs=WCS(hdu.header))
        fig.add_axes(ax)
        ax.imshow(hdu.data, vmin=-100, vmax=3000, origin='lower')
        # Overplot contour
        ax.contour(hdu_msx.data, transform=ax.get_transform(wcs_msx), colors='orange', levels=[2.5e-5, 5e-5, 1.e-4])
        ax.set_xlim(0., 720.)
        ax.set_ylim(0., 720.)

        self.generate_or_test(generate, fig, 'contour_overlay.png')

    # Test for overlaying grid, changing format of ticks, setting spacing and number of ticks
    @remote_data
    def test_overlay_features_image(self, generate):
        hdu = fits.open(self._image1)[0]
        fig = plt.figure(figsize=(6, 6))
        ax = WCSAxes(fig, [0.25, 0.25, 0.65, 0.65], wcs=WCS(hdu.header))
        fig.add_axes(ax)
        # Change the format of the ticks
        ax.coords[0].set_major_formatter('dd:mm:ss')
        ax.coords[1].set_major_formatter('dd:mm:ss.ssss')

        # Overlay grid on image
        ax.grid(color='red', alpha=1.0, lw=1, linestyle='dashed')

        # Set the spacing of ticks on the 'glon' axis to 4 arcsec
        ax.coords['glon'].set_ticks(spacing=4 * u.arcsec)
        # Set the number of ticks on the 'glat' axis to 9
        ax.coords['glat'].set_ticks(number=9)
        # Set labels on axes
        ax.coords['glon'].set_axislabel('Galactic Longitude')
        ax.coords['glat'].set_axislabel('Galactic Latitude')

        self.generate_or_test(generate, fig, 'overlay_features_image.png')

    # Overlay curvilinear grid and patches on image
    @remote_data
    def test_curvilinear_grid_patches_image(self, generate):
        hdu = fits.open(self._image2)[0]
        fig = plt.figure(figsize=(8, 8))
        ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8], wcs=WCS(hdu.header))
        fig.add_axes(ax)
        ax.imshow(hdu.data, vmax=1000, origin='lower')
        ax.grid(color='white', alpha=1.0, lw=1, linestyle='dashed')
        p = Circle((300, 100), radius=40, ec='yellow', fc='none')
        ax.add_patch(p)
        p = Circle((30., 20.), radius=20., ec='orange', fc='none', transform=ax.get_transform('world'))
        ax.add_patch(p)
        p = Circle((60., 50.), radius=20., ec='red', fc='none', transform=ax.get_transform('fk5'))
        ax.add_patch(p)

        self.generate_or_test(generate, fig, 'curvlinear_grid_patches_image.png')

    @remote_data
    def test_cube_slice_image(self, generate):
        image = fits.getdata(self._data_cube)
        w = WCS(self._data_cube)
        fig = plt.figure()
        ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8], w, slices=(50, 'y', 'x'))
        fig.add_axes(ax)
        ax.imshow(image[:, :, 100].transpose(), cmap=plt.cm.gist_heat)
        ax.coords[2].set_axislabel('Velocity m/s')

        self.generate_or_test(generate, fig, 'cube_slice_image.png')

    # Test for axes and ticks sizes, labels etc
    def test_ticks_labels(self, generate):
        fig = plt.figure(figsize=(6, 6))
        ax = WCSAxes(fig, [0.1, 0.1, 0.7, 0.7], wcs=None)
        fig.add_axes(ax)
        ax.set_xlim(-0.5, 2)
        ax.set_ylim(-0.5, 2)
        ax.coords[0].set_ticks(size=10, color='blue', alpha=0.2)
        ax.coords[1].set_ticks(size=20, color='red', alpha=0.9)
        ax.coords[0].set_ticks_position('all')
        ax.coords[1].set_ticks_position('all')
        ax.coords[0].set_axislabel('X-axis', size=20)
        ax.coords[1].set_axislabel('Y-axis', color='green', size=25, weight='regular', style='normal', family='monospace')
        ax.coords[0].set_axislabel_position('t')
        ax.coords[1].set_axislabel_position('r')
        ax.coords[0].set_ticklabel(color='purple', size=15, alpha=1, weight='light', style='normal', family='sans-serif')
        ax.coords[1].set_ticklabel(color='black', size=18, alpha=0.9, weight='bold', family='serif')
        ax.coords[0].set_ticklabel_position('all')
        ax.coords[1].set_ticklabel_position('r')

        self.generate_or_test(generate, fig, 'ticks_labels.png')
