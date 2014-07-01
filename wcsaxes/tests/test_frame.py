import numpy as np
import matplotlib.pyplot as plt

from astropy.wcs import WCS

from .. import WCSAxes
from .. import datasets
from ..frame import BaseFrame

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

    def test_custom_frame(self, generate):

        wcs = WCS(self.msx_header)

        fig = plt.figure(figsize=(4,4))

        ax = WCSAxes(fig, [0.15, 0.15, 0.7, 0.7],
                     wcs=wcs,
                     frame_class=HexagonalFrame)
        fig.add_axes(ax)

        ax.coords.grid(color='white')

        im = ax.imshow(np.ones((149, 149)), vmin=0., vmax=2., cmap=plt.cm.gist_heat,
                  origin='lower')

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

        self.generate_or_test(generate, fig, 'custom_frame.png')
