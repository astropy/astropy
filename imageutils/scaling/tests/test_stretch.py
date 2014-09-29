import numpy as np

from ..stretch import *


class TestStretch(object):

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
                                   np.array([0., 0.0625, 0.25, 0.5625, 1.]))

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
                                   np.array([0., 0.004628, 0.030653, 0.177005, 1.]),
                                   atol=1.e-6)

        np.testing.assert_allclose(stretch.inverted()(stretch(self.data)),
                                   self.data)

    def test_log(self):

        stretch = LogStretch()
        np.testing.assert_allclose(stretch(self.data),
                                   np.array([0., 0.799776, 0.899816, 0.958408, 1.]),
                                   atol=1.e-6)

        np.testing.assert_allclose(stretch.inverted()(stretch(self.data)),
                                   self.data)

    def test_asinh(self):

        stretch = AsinhStretch()
        np.testing.assert_allclose(stretch(self.data),
                                   np.array([0., 0.549402, 0.77127, 0.904691, 1.]),
                                   atol=1.e-6)

        np.testing.assert_allclose(stretch.inverted()(stretch(self.data)),
                                   self.data)

    def test_sinh(self):

        stretch = SinhStretch()
        np.testing.assert_allclose(stretch(self.data),
                                   np.array([0., 0.082085, 0.212548, 0.46828, 1.]),
                                   atol=1.e-6)

        np.testing.assert_allclose(stretch.inverted()(stretch(self.data)),
                                   self.data)

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
