# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest
import numpy as np
import matplotlib.pyplot as plt

from ....wcs import WCS

from .. import WCSAxes
from ..frame import BaseFrame

from ....tests.image_tests import IMAGE_REFERENCE_DIR
from .test_images import BaseImageTests


class HexagonalFrame(BaseFrame):

    spine_names = 'abcdef'

    def update_spines(self):

        xmin, xmax = self.parent_axes.get_xlim()
        ymin, ymax = self.parent_axes.get_ylim()

        ymid = 0.5 * (ymin + ymax)
        xmid1 = (xmin + xmax) / 4.
        xmid2 = (xmin + xmax) * 3. / 4.

        self['a'].data = np.array(([xmid1, ymin], [xmid2, ymin]))
        self['b'].data = np.array(([xmid2, ymin], [xmax, ymid]))
        self['c'].data = np.array(([xmax, ymid], [xmid2, ymax]))
        self['d'].data = np.array(([xmid2, ymax], [xmid1, ymax]))
        self['e'].data = np.array(([xmid1, ymax], [xmin, ymid]))
        self['f'].data = np.array(([xmin, ymid], [xmid1, ymin]))


class TestFrame(BaseImageTests):

    @pytest.mark.remote_data(source='astropy')
    @pytest.mark.mpl_image_compare(baseline_dir=IMAGE_REFERENCE_DIR,
                                   filename='custom_frame.png',
                                   tolerance=0, style={})
    def test_custom_frame(self):

        wcs = WCS(self.msx_header)

        fig = plt.figure(figsize=(4, 4))

        ax = WCSAxes(fig, [0.15, 0.15, 0.7, 0.7],
                     wcs=wcs,
                     frame_class=HexagonalFrame)
        fig.add_axes(ax)

        ax.coords.grid(color='white')

        im = ax.imshow(np.ones((149, 149)), vmin=0., vmax=2.,
                       origin='lower', cmap=plt.cm.gist_heat)

        minpad = {}
        minpad['a'] = minpad['d'] = 1
        minpad['b'] = minpad['c'] = minpad['e'] = minpad['f'] = 2.75

        ax.coords['glon'].set_axislabel("Longitude", minpad=minpad)
        ax.coords['glon'].set_axislabel_position('ad')

        ax.coords['glat'].set_axislabel("Latitude", minpad=minpad)
        ax.coords['glat'].set_axislabel_position('bcef')

        ax.coords['glon'].set_ticklabel_position('ad')
        ax.coords['glat'].set_ticklabel_position('bcef')

        # Set limits so that no labels overlap
        ax.set_xlim(5.5, 100.5)
        ax.set_ylim(5.5, 110.5)

        # Clip the image to the frame
        im.set_clip_path(ax.coords.frame.patch)

        return fig

    @pytest.mark.remote_data(source='astropy')
    @pytest.mark.mpl_image_compare(baseline_dir=IMAGE_REFERENCE_DIR,
                                   filename='update_clip_path_rectangular.png',
                                   tolerance=0, style={})
    def test_update_clip_path_rectangular(self, tmpdir):

        fig = plt.figure()
        ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8], aspect='equal')

        fig.add_axes(ax)

        ax.set_xlim(0., 2.)
        ax.set_ylim(0., 2.)

        # Force drawing, which freezes the clip path returned by WCSAxes
        fig.savefig(tmpdir.join('nothing').strpath)

        ax.imshow(np.zeros((12, 4)))

        ax.set_xlim(-0.5, 3.5)
        ax.set_ylim(-0.5, 11.5)

        return fig

    @pytest.mark.remote_data(source='astropy')
    @pytest.mark.mpl_image_compare(baseline_dir=IMAGE_REFERENCE_DIR,
                                   filename='update_clip_path_nonrectangular.png',
                                   tolerance=0, style={})
    def test_update_clip_path_nonrectangular(self, tmpdir):

        fig = plt.figure()
        ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8], aspect='equal',
                     frame_class=HexagonalFrame)

        fig.add_axes(ax)

        ax.set_xlim(0., 2.)
        ax.set_ylim(0., 2.)

        # Force drawing, which freezes the clip path returned by WCSAxes
        fig.savefig(tmpdir.join('nothing').strpath)

        ax.imshow(np.zeros((12, 4)))

        ax.set_xlim(-0.5, 3.5)
        ax.set_ylim(-0.5, 11.5)

        return fig

    @pytest.mark.remote_data(source='astropy')
    @pytest.mark.mpl_image_compare(baseline_dir=IMAGE_REFERENCE_DIR,
                                   filename='update_clip_path_change_wcs.png',
                                   tolerance=0, style={})
    def test_update_clip_path_change_wcs(self, tmpdir):

        # When WCS is changed, a new frame is created, so we need to make sure
        # that the path is carried over to the new frame.

        fig = plt.figure()
        ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8], aspect='equal')

        fig.add_axes(ax)

        ax.set_xlim(0., 2.)
        ax.set_ylim(0., 2.)

        # Force drawing, which freezes the clip path returned by WCSAxes
        fig.savefig(tmpdir.join('nothing').strpath)

        ax.reset_wcs()

        ax.imshow(np.zeros((12, 4)))

        ax.set_xlim(-0.5, 3.5)
        ax.set_ylim(-0.5, 11.5)

        return fig

    def test_copy_frame_properties_change_wcs(self):

        # When WCS is changed, a new frame is created, so we need to make sure
        # that the color and linewidth are transferred over

        fig = plt.figure()
        ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8])
        fig.add_axes(ax)
        ax.coords.frame.set_linewidth(5)
        ax.coords.frame.set_color('purple')
        ax.reset_wcs()
        assert ax.coords.frame.get_linewidth() == 5
        assert ax.coords.frame.get_color() == 'purple'
