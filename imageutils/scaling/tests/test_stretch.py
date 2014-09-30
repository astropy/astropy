import numpy as np

from astropy.tests.helper import pytest

from ..stretch import *

DATA = np.array([0.00, 0.25, 0.50, 0.75, 1.00])

RESULTS = {}
RESULTS[LinearStretch()] = np.array([0.00, 0.25, 0.50, 0.75, 1.00])
RESULTS[SqrtStretch()] = np.array([0., 0.5, 0.70710678, 0.8660254, 1.])
RESULTS[SquaredStretch()] = np.array([0., 0.0625, 0.25, 0.5625, 1.])
RESULTS[PowerStretch(0.5)] = np.array([0., 0.5, 0.70710678, 0.8660254, 1.])
RESULTS[PowerDistStretch()] = np.array([0., 0.004628, 0.030653, 0.177005, 1.])
RESULTS[LogStretch()] = np.array([0., 0.799776, 0.899816, 0.958408, 1.])
RESULTS[AsinhStretch()] = np.array([0., 0.549402, 0.77127, 0.904691, 1.])
RESULTS[SinhStretch()] = np.array([0., 0.082085, 0.212548, 0.46828, 1.])
RESULTS[ContrastBiasStretch(contrast=2., bias=0.4)] = np.array([-0.3, 0.2, 0.7, 1.2, 1.7])
RESULTS[HistEqStretch(DATA)] = DATA
RESULTS[HistEqStretch(DATA[::-1])] = DATA
RESULTS[HistEqStretch(DATA ** 0.5)] = np.array([0., 0.125, 0.25, 0.5674767, 1.])


class TestStretch(object):

    def setup_class(self):
        DATA = np.array([0.00, 0.25, 0.50, 0.75, 1.00])


    @pytest.mark.parametrize('stretch', RESULTS.keys())
    def test_no_clip(self, stretch):

        np.testing.assert_allclose(stretch(DATA),
                                           RESULTS[stretch], atol=1.e-6)

    @pytest.mark.parametrize('stretch', RESULTS.keys())
    def test_clip(self, stretch):

        np.testing.assert_allclose(stretch(DATA, clip=True),
                                   np.clip(RESULTS[stretch], 0., 1), atol=1.e-6)

    @pytest.mark.parametrize('stretch', RESULTS.keys())
    def test_inplace(self, stretch):

        data_in = DATA.copy()
        result = np.zeros(DATA.shape)
        stretch(data_in, out=result)
        np.testing.assert_allclose(result, RESULTS[stretch], atol=1.e-6)
        np.testing.assert_allclose(data_in, DATA)

    @pytest.mark.parametrize('stretch', RESULTS.keys())
    def test_round_trip(self, stretch):

        np.testing.assert_allclose(stretch.inverted()(stretch(DATA)),
                                   DATA)

    @pytest.mark.parametrize('stretch', RESULTS.keys())
    def test_inplace_roundtrip(self, stretch):

        result = np.zeros(DATA.shape)
        stretch(DATA, out=result)
        stretch.inverted()(result, out=result)
        np.testing.assert_allclose(result, DATA)

    @pytest.mark.parametrize('stretch', RESULTS.keys())
    def test_double_inverse(self, stretch):
        np.testing.assert_allclose(stretch.inverted().inverted()(DATA), stretch(DATA), atol=1.e-6)

    def test_inverted(self):
        stretch_1 = SqrtStretch().inverted()
        stretch_2 = PowerStretch(2)
        np.testing.assert_allclose(stretch_1(DATA),
                                   stretch_2(DATA))

    def test_chaining(self):

        stretch_1 = SqrtStretch() + SqrtStretch()
        stretch_2 = PowerStretch(0.25)
        stretch_3 = PowerStretch(4.)

        np.testing.assert_allclose(stretch_1(DATA),
                                   stretch_2(DATA))

        np.testing.assert_allclose(stretch_1.inverted()(DATA),
                                   stretch_3(DATA))
