# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest
import numpy as np

from ...utils import NumpyRNGContext

from ..interval import (ManualInterval, MinMaxInterval, PercentileInterval,
                        AsymmetricPercentileInterval, ZScaleInterval)


class TestInterval:

    data = np.linspace(-20., 60., 100)

    def test_manual(self):
        interval = ManualInterval(-10., +15.)
        vmin, vmax = interval.get_limits(self.data)
        np.testing.assert_allclose(vmin, -10.)
        np.testing.assert_allclose(vmax, +15.)

    def test_manual_defaults(self):

        interval = ManualInterval(vmin=-10.)
        vmin, vmax = interval.get_limits(self.data)
        np.testing.assert_allclose(vmin, -10.)
        np.testing.assert_allclose(vmax, np.max(self.data))

        interval = ManualInterval(vmax=15.)
        vmin, vmax = interval.get_limits(self.data)
        np.testing.assert_allclose(vmin, np.min(self.data))
        np.testing.assert_allclose(vmax, 15.)

    def test_manual_zero_limit(self):
        # Regression test for a bug that caused ManualInterval to compute the
        # limit (min or max) if it was set to zero.
        interval = ManualInterval(vmin=0, vmax=0)
        vmin, vmax = interval.get_limits(self.data)
        np.testing.assert_allclose(vmin, 0)
        np.testing.assert_allclose(vmax, 0)

    def test_manual_defaults_with_nan(self):
        interval = ManualInterval()
        data = np.copy(self.data)
        data[0] = np.nan
        vmin, vmax = interval.get_limits(self.data)
        np.testing.assert_allclose(vmin, -20)
        np.testing.assert_allclose(vmax, +60)

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

    def test_asymmetric_percentile_nsamples(self):
        with NumpyRNGContext(12345):
            interval = AsymmetricPercentileInterval(10.5, 70.5, n_samples=20)
            vmin, vmax = interval.get_limits(self.data)
        np.testing.assert_allclose(vmin, -14.367676767676768)
        np.testing.assert_allclose(vmax, 40.266666666666666)


class TestIntervalList(TestInterval):

    # Make sure intervals work with lists
    data = np.linspace(-20., 60., 100).tolist()


class TestInterval2D(TestInterval):

    # Make sure intervals work with 2d arrays
    data = np.linspace(-20., 60., 100).reshape(100, 1)


def test_zscale():
    np.random.seed(42)
    data = np.random.randn(100, 100) * 5 + 10
    interval = ZScaleInterval()
    vmin, vmax = interval.get_limits(data)
    np.testing.assert_allclose(vmin, -9.6, atol=0.1)
    np.testing.assert_allclose(vmax, 25.4, atol=0.1)

    data = list(range(1000)) + [np.nan]
    interval = ZScaleInterval()
    vmin, vmax = interval.get_limits(data)
    np.testing.assert_allclose(vmin, 0, atol=0.1)
    np.testing.assert_allclose(vmax, 999, atol=0.1)

    data = list(range(100))
    interval = ZScaleInterval()
    vmin, vmax = interval.get_limits(data)
    np.testing.assert_allclose(vmin, 0, atol=0.1)
    np.testing.assert_allclose(vmax, 99, atol=0.1)


def test_integers():
    # Need to make sure integers get cast to float
    interval = MinMaxInterval()
    values = interval([1, 3, 4, 5, 6])
    np.testing.assert_allclose(values, [0., 0.4, 0.6, 0.8, 1.0])

    # Don't accept integer array in output
    out = np.zeros(5, dtype=int)
    with pytest.raises(TypeError) as exc:
        values = interval([1, 3, 4, 5, 6], out=out)
    assert exc.value.args[0] == ("Can only do in-place scaling for "
                                 "floating-point arrays")

    # But integer input and floating point output is fine
    out = np.zeros(5, dtype=float)
    interval([1, 3, 4, 5, 6], out=out)
    np.testing.assert_allclose(out, [0., 0.4, 0.6, 0.8, 1.0])


def test_constant_data():
    """Test intervals with constant data (avoiding divide-by-zero)."""
    shape = (10, 10)
    data = np.ones(shape)
    interval = MinMaxInterval()
    limits = interval.get_limits(data)
    values = interval(data)
    np.testing.assert_allclose(limits, (1., 1.))
    np.testing.assert_allclose(values, np.zeros(shape))
