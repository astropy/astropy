# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_equal

from astropy.utils.exceptions import AstropyDeprecationWarning
from astropy.visualization.stretch import (
    AsinhStretch,
    ContrastBiasStretch,
    HistEqStretch,
    InvertedHistEqStretch,
    InvertedLogStretch,
    InvertedPowerDistStretch,
    LinearStretch,
    LogStretch,
    PowerDistStretch,
    PowerStretch,
    SinhStretch,
    SqrtStretch,
    SquaredStretch,
)

DATA = np.array([0.00, 0.25, 0.50, 0.75, 1.00])

RESULTS = {}
RESULTS[LinearStretch()] = np.array([0.00, 0.25, 0.50, 0.75, 1.00])
RESULTS[LinearStretch(intercept=0.5) + LinearStretch(slope=0.5)] = np.array(
    [0.5, 0.625, 0.75, 0.875, 1.0]
)
RESULTS[SqrtStretch()] = np.array([0.0, 0.5, 0.70710678, 0.8660254, 1.0])
RESULTS[SquaredStretch()] = np.array([0.0, 0.0625, 0.25, 0.5625, 1.0])
RESULTS[PowerStretch(0.5)] = np.array([0.0, 0.5, 0.70710678, 0.8660254, 1.0])
RESULTS[PowerDistStretch()] = np.array([0.0, 0.004628, 0.030653, 0.177005, 1.0])
RESULTS[LogStretch()] = np.array([0.0, 0.799776, 0.899816, 0.958408, 1.0])
RESULTS[AsinhStretch()] = np.array([0.0, 0.549402, 0.77127, 0.904691, 1.0])
RESULTS[SinhStretch()] = np.array([0.0, 0.082085, 0.212548, 0.46828, 1.0])
RESULTS[ContrastBiasStretch(contrast=2.0, bias=0.4)] = np.array(
    [-0.3, 0.2, 0.7, 1.2, 1.7]
)
RESULTS[HistEqStretch(DATA)] = DATA
RESULTS[HistEqStretch(DATA[::-1])] = DATA
RESULTS[HistEqStretch(DATA**0.5)] = np.array([0.0, 0.125, 0.25, 0.5674767, 1.0])


class TestStretch:
    @pytest.mark.parametrize("stretch", RESULTS.keys())
    def test_no_clip(self, stretch):
        np.testing.assert_allclose(
            stretch(DATA, clip=False), RESULTS[stretch], atol=1.0e-6
        )

    @pytest.mark.parametrize("ndim", [2, 3])
    @pytest.mark.parametrize("stretch", RESULTS.keys())
    def test_clip_ndimensional(self, stretch, ndim):
        new_shape = DATA.shape + (1,) * ndim

        np.testing.assert_allclose(
            stretch(DATA.reshape(new_shape), clip=True).ravel(),
            np.clip(RESULTS[stretch], 0.0, 1),
            atol=1.0e-6,
        )

    @pytest.mark.parametrize("stretch", RESULTS.keys())
    def test_clip(self, stretch):
        np.testing.assert_allclose(
            stretch(DATA, clip=True), np.clip(RESULTS[stretch], 0.0, 1), atol=1.0e-6
        )

    @pytest.mark.parametrize("stretch", RESULTS.keys())
    def test_inplace(self, stretch):
        data_in = DATA.copy()
        result = np.zeros(DATA.shape)
        stretch(data_in, out=result, clip=False)
        np.testing.assert_allclose(result, RESULTS[stretch], atol=1.0e-6)
        np.testing.assert_allclose(data_in, DATA)

    @pytest.mark.parametrize("stretch", RESULTS.keys())
    def test_round_trip(self, stretch):
        np.testing.assert_allclose(
            stretch.inverse(stretch(DATA, clip=False), clip=False), DATA
        )

    @pytest.mark.parametrize("stretch", RESULTS.keys())
    def test_inplace_roundtrip(self, stretch):
        result = np.zeros(DATA.shape)
        stretch(DATA, out=result, clip=False)
        stretch.inverse(result, out=result, clip=False)
        np.testing.assert_allclose(result, DATA)

    @pytest.mark.parametrize("stretch", RESULTS.keys())
    def test_double_inverse(self, stretch):
        np.testing.assert_allclose(
            stretch.inverse.inverse(DATA), stretch(DATA), atol=1.0e-6
        )

    def test_inverted(self):
        stretch_1 = SqrtStretch().inverse
        stretch_2 = PowerStretch(2)
        np.testing.assert_allclose(stretch_1(DATA), stretch_2(DATA))

    def test_chaining(self):
        stretch_1 = SqrtStretch() + SqrtStretch()
        stretch_2 = PowerStretch(0.25)
        stretch_3 = PowerStretch(4.0)

        np.testing.assert_allclose(stretch_1(DATA), stretch_2(DATA))

        np.testing.assert_allclose(stretch_1.inverse(DATA), stretch_3(DATA))


def test_clip_invalid():
    stretch = SqrtStretch()

    values = stretch([-1.0, 0.0, 0.5, 1.0, 1.5])
    np.testing.assert_allclose(values, [0.0, 0.0, 0.70710678, 1.0, 1.0])

    values = stretch([-1.0, 0.0, 0.5, 1.0, 1.5], clip=False)
    np.testing.assert_allclose(values, [np.nan, 0.0, 0.70710678, 1.0, 1.2247448])


@pytest.mark.parametrize("a", [-2.0, -1, 1.0])
def test_invalid_powerdist_a(a):
    match = "a must be >= 0, but cannot be set to 1"
    with pytest.raises(ValueError, match=match):
        PowerDistStretch(a=a)
    with pytest.raises(ValueError, match=match):
        InvertedPowerDistStretch(a=a)


@pytest.mark.parametrize("a", [-2.0, -1, 0.0])
def test_invalid_power_log_a(a):
    match = "a must be > 0"
    with pytest.raises(ValueError, match=match):
        PowerStretch(a=a)
    with pytest.raises(ValueError, match=match):
        LogStretch(a=a)
    with pytest.raises(ValueError, match=match):
        InvertedLogStretch(a=a)


@pytest.mark.parametrize("a", [-2.0, -1, 0.0])
def test_invalid_sinh_a(a):
    match = "a must be > 0"
    with pytest.raises(ValueError, match=match):
        AsinhStretch(a=a)
    with pytest.raises(ValueError, match=match):
        SinhStretch(a=a)


def test_sinh_a():
    a = 0.9
    a_inv = 1.0 / np.arcsinh(1.0 / a)
    z = AsinhStretch(a=a)
    assert_allclose(z.inverse.a, a_inv)


def test_histeqstretch_invalid():
    data = np.array([-np.inf, 0.00, 0.25, 0.50, 0.75, 1.00, np.inf])
    result = np.array([0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0])
    assert_equal(HistEqStretch(data)(data), result)
    assert_equal(InvertedHistEqStretch(data)(data), result)


def test_deprecated_attrs():
    match = "The power attribute is deprecated"
    with pytest.warns(AstropyDeprecationWarning, match=match):
        stretch = PowerStretch(a=0.5)
        assert stretch.power == stretch.a

    match = "The exp attribute is deprecated"
    with pytest.warns(AstropyDeprecationWarning, match=match):
        stretch = PowerDistStretch(a=0.5)
        assert stretch.exp == stretch.a
    with pytest.warns(AstropyDeprecationWarning, match=match):
        stretch = InvertedPowerDistStretch(a=0.5)
        assert stretch.exp == stretch.a

    with pytest.warns(AstropyDeprecationWarning, match=match):
        stretch = LogStretch(a=0.5)
        assert stretch.exp == stretch.a
    with pytest.warns(AstropyDeprecationWarning, match=match):
        stretch = InvertedLogStretch(a=0.5)
        assert stretch.exp == stretch.a
