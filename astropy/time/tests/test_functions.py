# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest
import numpy as np

from astropy.time import Time, TimeDelta
from astropy.units.quantity_helper.function_helpers import ARRAY_FUNCTION_ENABLED


class TestFunctionsTime:

    def setup_class(cls):
        cls.t = Time(50000, np.arange(8).reshape(4, 2), format='mjd',
                     scale='tai')

    def check(self, func, cls=None, scale=None, format=None, *args, **kwargs):
        if cls is None:
            cls = self.t.__class__
        if scale is None:
            scale = self.t.scale
        if format is None:
            format = self.t.format
        out = func(self.t, *args, **kwargs)
        jd1 = func(self.t.jd1, *args, **kwargs)
        jd2 = func(self.t.jd2, *args, **kwargs)
        expected = cls(jd1, jd2, format=format, scale=scale)
        if isinstance(out, np.ndarray):
            expected = np.array(expected)

        assert np.all(out == expected)

    @pytest.mark.parametrize('axis', (0, 1))
    def test_diff(self, axis):
        self.check(np.diff, axis=axis, cls=TimeDelta, format='jd')


class TestFunctionsTimeDelta(TestFunctionsTime):

    def setup_class(cls):
        cls.t = TimeDelta(np.arange(8).reshape(4, 2), format='jd',
                          scale='tai')

    @pytest.mark.parametrize('axis', (0, 1, None))
    @pytest.mark.parametrize('func', (np.sum, np.mean, np.median))
    def test_sum_like(self, func, axis):
        self.check(func, axis=axis)


@pytest.mark.xfail(not ARRAY_FUNCTION_ENABLED,
                   reason="Needs __array_function__ support")
@pytest.mark.parametrize('attribute', ['shape', 'ndim', 'size'])
@pytest.mark.parametrize('t', [
    Time('2001-02-03T04:05:06'),
    Time(50000, np.arange(8).reshape(4, 2), format='mjd', scale='tai'),
    TimeDelta(100, format='jd')])
def test_shape_attribute_functions(t, attribute):
    # Regression test for
    # https://github.com/astropy/astropy/issues/8610#issuecomment-736855217
    function = getattr(np, attribute)
    result = function(t)
    assert result == getattr(t, attribute)
