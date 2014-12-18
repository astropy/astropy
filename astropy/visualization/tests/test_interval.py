# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np

from ...tests.helper import pytest

from ..interval import (ManualInterval, MinMaxInterval, PercentileInterval,
                        AsymmetricPercentileInterval)


class TestInterval(object):

    data = np.linspace(-20., 60., 100)

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


class TestIntervalList(TestInterval):

    # Make sure intervals work with lists

    data = np.linspace(-20., 60., 100).tolist()


def test_integers():

    # Need to make sure integers get cast to float
    interval = MinMaxInterval()
    values = interval([1, 3, 4, 5, 6])
    np.testing.assert_allclose(values, [0., 0.4, 0.6, 0.8, 1.0])

    # Don't accept integer array in output
    out = np.zeros(5, dtype=int)
    with pytest.raises(TypeError) as exc:
        values = interval([1, 3, 4, 5, 6], out=out)
    assert exc.value.args[0] == "Can only do in-place scaling for floating-point arrays"

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
