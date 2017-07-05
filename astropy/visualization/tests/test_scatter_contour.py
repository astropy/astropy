# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

try:
    import matplotlib.pyplot as plt
    HAS_PLT = True
except ImportError:
    HAS_PLT = False

import pytest
import numpy as np

from ..scatter_contour import scatter_contour

class TestScatterContour(object):

    def setup(self):
        rng = np.random.RandomState(42)

        # need many points for contour - arbitrary locations and scales
        self.x, self.y = rng.normal([-0.4, 15.,], [12., 0.8], size=(1024,2)).T

    @pytest.mark.skipif('not HAS_PLT')
    def test_scatter_contour_basic(self):

        # lower the threshold so it succeeds
        scatter_contour(self.x, self.y, threshold=10)

        # or, could make bins larger
        scatter_contour(self.x, self.y, threshold=100,
                        histogram2d_kwargs=dict(bins=5))

        # fails when threshold is too high, no valid contours to draw
        with pytest.raises(ValueError):
            scatter_contour(self.x, self.y, threshold=100)

    @pytest.mark.skipif('not HAS_PLT')
    def test_hist_specify_ax(self):

        fig, axes = plt.subplots(2)
        _,contours = scatter_contour(self.x, self.y, threshold=10, ax=axes[0])
        assert contours.ax is axes[0]

        _,contours = scatter_contour(self.x, self.y, threshold=10, ax=axes[1])
        assert contours.ax is axes[1]

    @pytest.mark.skipif('not HAS_PLT')
    def test_hist_kwargs(self):

        # test passing arguments
        scatter_contour(self.x, self.y, threshold=10, levels=4)
        scatter_contour(self.x, self.y, threshold=10,
                        levels=np.linspace(8, 100, 8))

        scatter_contour(self.x, self.y, threshold=10, filled_contour=False)

        scatter_contour(self.x, self.y, threshold=10, levels=4)
        scatter_contour(self.x, self.y, threshold=10, levels=3, log_counts=True)

        scatter_contour(self.x, self.y, threshold=2,
                        plot_kwargs=dict(marker='o'),
                        contour_kwargs=dict(cmap='Greys'),
                        histogram2d_kwargs=dict(bins=(16,16)))
