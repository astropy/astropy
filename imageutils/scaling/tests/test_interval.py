import numpy as np

from ..interval import *


class TestInterval(object):

    def setup_class(self):
        self.data = np.linspace(-20., 60., 100)

    def test_manual(self):
        interval = ManualInterval(-10., +15.)
        vmin, vmax = interval.get_limits(self.data)
        np.testing.assert_allclose(vmin, -10.)
        np.testing.assert_allclose(vmax, +15.)

    def test_minmax(self):
        interval = MinMaxInterval()
        vmin, vmax = interval.get_limits(self.data)
        np.testing.assert_allclose(vmin, -20.)
        np.testing.assert_allclose(vmax, +60.)

    def test_percentile(self):
        interval = PercentileInterval(62.2)
        vmin, vmax = interval.get_limits(self.data)
        np.testing.assert_allclose(vmin, -4.88)
        np.testing.assert_allclose(vmax, 44.88)

    def test_asymmetric_percentile(self):
        interval = AsymmetricPercentileInterval(10.5, 70.5)
        vmin, vmax = interval.get_limits(self.data)
        np.testing.assert_allclose(vmin, -11.6)
        np.testing.assert_allclose(vmax, 36.4)
