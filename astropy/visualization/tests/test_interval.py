# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
import pytest
from numpy.testing import assert_allclose

from astropy.utils import NumpyRNGContext
from astropy.utils.masked import Masked
from astropy.visualization.interval import (
    AsymmetricPercentileInterval,
    ManualInterval,
    MinMaxInterval,
    PercentileInterval,
    ZScaleInterval,
)


class TestInterval:
    data = np.linspace(-20.0, 60.0, 100)

    def test_manual(self):
        interval = ManualInterval(-10.0, +15.0)
        vmin, vmax = interval.get_limits(self.data)
        assert_allclose(vmin, -10.0)
        assert_allclose(vmax, +15.0)

    def test_manual_defaults(self):
        interval = ManualInterval(vmin=-10.0)
        vmin, vmax = interval.get_limits(self.data)
        assert_allclose(vmin, -10.0)
        assert_allclose(vmax, np.max(self.data))

        interval = ManualInterval(vmax=15.0)
        vmin, vmax = interval.get_limits(self.data)
        assert_allclose(vmin, np.min(self.data))
        assert_allclose(vmax, 15.0)

    def test_manual_zero_limit(self):
        # Regression test for a bug that caused ManualInterval to compute the
        # limit (min or max) if it was set to zero.
        interval = ManualInterval(vmin=0, vmax=0)
        vmin, vmax = interval.get_limits(self.data)
        assert_allclose(vmin, 0)
        assert_allclose(vmax, 0)

    def test_manual_defaults_with_nan(self):
        interval = ManualInterval()
        data = np.copy(self.data)
        data[0] = np.nan
        vmin, vmax = interval.get_limits(self.data)
        assert_allclose(vmin, -20)
        assert_allclose(vmax, +60)

    def test_minmax(self):
        interval = MinMaxInterval()
        vmin, vmax = interval.get_limits(self.data)
        assert_allclose(vmin, -20.0)
        assert_allclose(vmax, +60.0)

    def test_percentile(self):
        interval = PercentileInterval(62.2)
        vmin, vmax = interval.get_limits(self.data)
        assert_allclose(vmin, -4.88)
        assert_allclose(vmax, 44.88)

    def test_asymmetric_percentile(self):
        interval = AsymmetricPercentileInterval(10.5, 70.5)
        vmin, vmax = interval.get_limits(self.data)
        assert_allclose(vmin, -11.6)
        assert_allclose(vmax, 36.4)

    def test_asymmetric_percentile_nsamples(self):
        with NumpyRNGContext(12345):
            interval = AsymmetricPercentileInterval(10.5, 70.5, n_samples=20)
            vmin, vmax = interval.get_limits(self.data)
        assert_allclose(vmin, -14.367676767676768)
        assert_allclose(vmax, 40.266666666666666)


class TestIntervalList(TestInterval):
    # Make sure intervals work with lists
    data = np.linspace(-20.0, 60.0, 100).tolist()


class TestInterval2D(TestInterval):
    # Make sure intervals work with 2d arrays
    data = np.linspace(-20.0, 60.0, 100).reshape(100, 1)


class TestIntervalMaskedArray(TestInterval):
    # Make sure intervals work with MaskedArray
    data = np.concatenate((np.linspace(-20.0, 60.0, 100), np.full(100, 1e6)))
    data = np.ma.MaskedArray(data, data > 1000)


class TestIntervalMaskedNDArray(TestInterval):
    # Make sure intervals work with MaskedArray
    data = np.concatenate((np.linspace(-20.0, 60.0, 100), np.full(100, 1e6)))
    data = Masked(data, data > 1000)


def test_zscale():
    np.random.seed(42)
    data = np.random.randn(100, 100) * 5 + 10
    interval = ZScaleInterval()
    vmin, vmax = interval.get_limits(data)
    assert_allclose(vmin, -9.6, atol=0.1)
    assert_allclose(vmax, 25.4, atol=0.1)

    data = list(range(1000)) + [np.nan]
    interval = ZScaleInterval()
    vmin, vmax = interval.get_limits(data)
    assert_allclose(vmin, 0, atol=0.1)
    assert_allclose(vmax, 999, atol=0.1)

    data = list(range(100))
    interval = ZScaleInterval()
    vmin, vmax = interval.get_limits(data)
    assert_allclose(vmin, 0, atol=0.1)
    assert_allclose(vmax, 99, atol=0.1)


def test_zscale_npoints():
    """
    Regression test to ensure ZScaleInterval returns the minimum and
    maximum of the data if the number of data points is less than
    ``min_pixels``.
    """

    data = np.arange(4).reshape((2, 2))
    interval = ZScaleInterval(min_npixels=5)
    vmin, vmax = interval.get_limits(data)
    assert vmin == 0
    assert vmax == 3


def test_integers():
    # Need to make sure integers get cast to float
    interval = MinMaxInterval()
    values = interval([1, 3, 4, 5, 6])
    assert_allclose(values, [0.0, 0.4, 0.6, 0.8, 1.0])

    # Don't accept integer array in output
    out = np.zeros(5, dtype=int)
    with pytest.raises(
        TypeError, match=r"Can only do in-place scaling for floating-point arrays"
    ):
        values = interval([1, 3, 4, 5, 6], out=out)

    # But integer input and floating point output is fine
    out = np.zeros(5, dtype=float)
    interval([1, 3, 4, 5, 6], out=out)
    assert_allclose(out, [0.0, 0.4, 0.6, 0.8, 1.0])


def test_constant_data():
    """Test intervals with constant data (avoiding divide-by-zero)."""
    shape = (10, 10)
    data = np.ones(shape)
    interval = MinMaxInterval()
    limits = interval.get_limits(data)
    values = interval(data)
    assert_allclose(limits, (1.0, 1.0))
    assert_allclose(values, np.zeros(shape))
