# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Test numpy functions overridden with __array_function__

See test_methods.py for tests of numpy functions that call methods or
that are overridden by ShapedLikeNDArray.
"""

import numpy as np
import pytest

from astropy.time import Time, TimeDelta

from .test_methods import assert_time_all_equal


@pytest.mark.parametrize(
    "a",
    [
        Time("J2010"),
        Time(50000.0, [1.0, 2.0], format="mjd"),
        Time("2002-03-04T15:16:17.18", precision=5),
    ],
)
@pytest.mark.parametrize("shape", [None, (2, 2)])
def test_zeros_like(a, shape):
    """Test np.zeros_like __array_function__ implementation."""
    res = np.zeros_like(a, shape=shape)
    assert np.all(res.jd == 2451544.5)
    assert res.shape == (shape or a.shape)
    assert res.scale == a.scale
    assert res.format == a.format
    assert res.precision == a.precision


class TestConcatenate:
    @classmethod
    def setup_class(cls):
        cls.a = Time(50000.0, np.arange(6.0).reshape(2, 3), format="mjd")
        cls.b = Time(["2010-11-12", "2011-10-09"]).reshape(2, 1)

    def check(self, func, *args, **kwargs):
        # Check assumes no different locations, etc.
        t_list = kwargs.pop("t_list", [self.a, self.a])
        out = func(t_list, *args, **kwargs)
        exp_jd1 = func([t.jd1 for t in t_list], *args, **kwargs)
        exp_jd2 = func([t.jd2 for t in t_list], *args, **kwargs)
        exp = Time(exp_jd1, exp_jd2, format="jd")
        assert_time_all_equal(out, exp)

    def test_concatenate(self):
        self.check(np.concatenate)
        self.check(np.concatenate, axis=1)
        self.check(np.concatenate, t_list=[self.a, self.b], axis=1)

    def test_hstack(self):
        self.check(np.hstack)

    def test_vstack(self):
        self.check(np.vstack)

    def test_dstack(self):
        self.check(np.dstack)


def test_linspace():
    """Test `np.linspace` `__array_func__` implementation for scalar and arrays."""
    t1 = Time(["2021-01-01 00:00:00", "2021-01-02 00:00:00"])
    t2 = Time(["2021-01-01 01:00:00", "2021-12-28 00:00:00"])
    atol = 2 * np.finfo(float).eps * abs(t1 - t2).max()

    ts = np.linspace(t1[0], t2[0], 3)
    assert ts[0].isclose(Time("2021-01-01 00:00:00"), atol=atol)
    assert ts[1].isclose(Time("2021-01-01 00:30:00"), atol=atol)
    assert ts[2].isclose(Time("2021-01-01 01:00:00"), atol=atol)

    ts = np.linspace(t1, t2[0], 2, endpoint=False)
    assert ts.shape == (2, 2)
    assert all(
        ts[0].isclose(Time(["2021-01-01 00:00:00", "2021-01-02 00:00:00"]), atol=atol)
    )
    assert all(
        ts[1].isclose(Time(["2021-01-01 00:30:00", "2021-01-01 12:30:00"]), atol=atol)
    )

    ts = np.linspace(t1, t2, 7)
    assert ts.shape == (7, 2)
    assert all(
        ts[0].isclose(Time(["2021-01-01 00:00:00", "2021-01-02 00:00:00"]), atol=atol)
    )
    assert all(
        ts[1].isclose(Time(["2021-01-01 00:10:00", "2021-03-03 00:00:00"]), atol=atol)
    )
    assert all(
        ts[5].isclose(Time(["2021-01-01 00:50:00", "2021-10-29 00:00:00"]), atol=atol)
    )
    assert all(
        ts[6].isclose(Time(["2021-01-01 01:00:00", "2021-12-28 00:00:00"]), atol=atol)
    )


def test_linspace_steps():
    """Test `np.linspace` `retstep` option."""
    t1 = Time(["2021-01-01 00:00:00", "2021-01-01 12:00:00"])
    t2 = Time("2021-01-02 00:00:00")
    atol = 2 * np.finfo(float).eps * abs(t1 - t2).max()

    ts, st = np.linspace(t1, t2, 7, retstep=True)
    assert ts.shape == (7, 2)
    assert st.shape == (2,)
    assert all(ts[1].isclose(ts[0] + st, atol=atol))
    assert all(ts[6].isclose(ts[0] + 6 * st, atol=atol))
    assert all(st.isclose(TimeDelta([14400, 7200], format="sec"), atol=atol))


def test_linspace_fmts():
    """Test `np.linspace` `__array_func__` implementation for start/endpoints
    from different formats/systems.
    """
    t1 = Time(["2020-01-01 00:00:00", "2020-01-02 00:00:00"])
    t2 = Time(2458850, format="jd")
    t3 = Time(1578009600, format="unix")
    atol = 2 * np.finfo(float).eps * abs(t1 - Time([t2, t3])).max()

    ts = np.linspace(t1, t2, 3)
    assert ts.shape == (3, 2)
    assert all(
        ts[0].isclose(Time(["2020-01-01 00:00:00", "2020-01-02 00:00:00"]), atol=atol)
    )
    assert all(
        ts[1].isclose(Time(["2020-01-01 06:00:00", "2020-01-01 18:00:00"]), atol=atol)
    )
    assert all(
        ts[2].isclose(Time(["2020-01-01 12:00:00", "2020-01-01 12:00:00"]), atol=atol)
    )

    ts = np.linspace(t1, Time([t2, t3]), 3)
    assert ts.shape == (3, 2)
    assert all(
        ts[0].isclose(Time(["2020-01-01 00:00:00", "2020-01-02 00:00:00"]), atol=atol)
    )
    assert all(
        ts[1].isclose(Time(["2020-01-01 06:00:00", "2020-01-02 12:00:00"]), atol=atol)
    )
    assert all(
        ts[2].isclose(Time(["2020-01-01 12:00:00", "2020-01-03 00:00:00"]), atol=atol)
    )
