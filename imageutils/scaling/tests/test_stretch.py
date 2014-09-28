import numpy as np

from ..stretch import *


class TestInterval(object):

    def setup_class(self):
        self.data = np.array([0.00, 0.25, 0.50, 0.75, 1.00])

    def test_linear(self):

        stretch = LinearStretch()
        np.testing.assert_allclose(stretch(self.data),
                                   self.data)

        np.testing.assert_allclose(stretch.inverted()(stretch(self.data)),
                                   self.data)

    def test_sqrt(self):

        stretch = SqrtStretch()
        np.testing.assert_allclose(stretch(self.data),
                                   np.array([0., 0.5, 0.70710678, 0.8660254, 1.]))

        np.testing.assert_allclose(stretch.inverted()(stretch(self.data)),
                                   self.data)

    def test_squared(self):

        stretch = SquaredStretch()
        np.testing.assert_allclose(stretch(self.data),
                                   np.array([0., 0.5, 0.70710678, 0.8660254, 1.]))

        np.testing.assert_allclose(stretch.inverted()(stretch(self.data)),
                                   self.data)

    def test_power(self):

        stretch = PowerStretch(0.5)
        np.testing.assert_allclose(stretch(self.data),
                                   np.array([0., 0.5, 0.70710678, 0.8660254, 1.]))

        np.testing.assert_allclose(stretch.inverted()(stretch(self.data)),
                                   self.data)

    def test_power_dist(self):

        stretch = PowerDistStretch()
        np.testing.assert_allclose(stretch(self.data),
                                   np.array([0., 0.00462341, 0.03062278, 0.17682794, 0.999]),
                                   rtol=1e-3)

        np.testing.assert_allclose(stretch.inverted()(stretch(self.data)),
                                   self.data, rtol=1e-3)

    def test_log(self):

        stretch = LogStretch()
        np.testing.assert_allclose(stretch(self.data),
                                   np.array([ 0., 0.79989124, 0.89994591, 0.95854665, 1.]),
                                   rtol=1e-3)

        np.testing.assert_allclose(stretch.inverted()(stretch(self.data)),
                                   self.data, rtol=1e-2)

    def test_asinh(self):

        stretch = AsinhStretch()
        np.testing.assert_allclose(stretch(self.data),
                                   np.array([0., 0.54907705, 0.77081278, 0.9041551, 0.99940765]),
                                   rtol=1e-3)

        np.testing.assert_allclose(stretch.inverted()(stretch(self.data)),
                                   self.data, rtol=1e-3)

    def test_sinh(self):

        stretch = SinhStretch()
        np.testing.assert_allclose(stretch(self.data),
                                   np.array([0., 0.08223167, 0.21292795, 0.46911683, 1.]))

        np.testing.assert_allclose(stretch.inverted()(stretch(self.data)),
                                   self.data, rtol=1e-3)

    def test_histeq(self):

        stretch = HistEqStretch(self.data)
        np.testing.assert_allclose(stretch(self.data),
                                   self.data)

        np.testing.assert_allclose(stretch.inverted()(stretch(self.data)),
                                   self.data)

        stretch = HistEqStretch(self.data[::-1])
        np.testing.assert_allclose(stretch(self.data),
                                   self.data)

        np.testing.assert_allclose(stretch.inverted()(stretch(self.data)),
                                   self.data)

        stretch = HistEqStretch(self.data ** 0.5)
        np.testing.assert_allclose(stretch(self.data),
                                   np.array([0., 0.125, 0.25, 0.5674767, 1.]))

        np.testing.assert_allclose(stretch.inverted()(stretch(self.data)),
                                   self.data)

    def test_contrast_bias(self):
        stretch = ContrastBiasStretch(contrast=2., bias=0.4)
        np.testing.assert_allclose(stretch(self.data),
                                   np.array([0., 0.2, 0.7, 1., 1.]))

        # This stretch does not invert cleanly because of clipping

    def test_inverted(self):
        stretch_1 = SqrtStretch().inverted()
        stretch_2 = PowerStretch(2)
        np.testing.assert_allclose(stretch_1(self.data),
                                   stretch_2(self.data))

    def test_chaining(self):
        stretch_1 = SqrtStretch() + SqrtStretch()
        stretch_2 = PowerStretch(0.25)
        np.testing.assert_allclose(stretch_1(self.data),
                                   stretch_2(self.data))
