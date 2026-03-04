"""
Test ``allclose`` and ``isclose``.
``allclose`` was ``quantity_allclose`` in ``astropy.tests.helper``.
"""

import numpy as np
import pytest

from astropy import units as u
from astropy.units.quantity_helper import function_helpers as fh


@pytest.mark.parametrize(
    ("a", "b"),
    [
        ([1, 2], [1, 2]),
        ([1, 2] * u.m, [100, 200] * u.cm),
        (1 * u.s, 1000 * u.ms),
    ],
)
def test_allclose_isclose_default(a, b):
    assert u.allclose(a, b)
    assert np.all(u.isclose(a, b))


def test_allclose_isclose():
    a = [1, 2] * u.m
    b = [101, 201] * u.cm
    delta = 2 * u.cm
    assert u.allclose(a, b, atol=delta)
    assert np.all(u.isclose(a, b, atol=delta))

    c = [90, 200] * u.cm
    assert not u.allclose(a, c)
    assert not np.all(u.isclose(a, c))


def test_unwrap_arange_args_stop_none_raises_runtimeerror_not_asserts():
    """
    Regression: production assert replaced with explicit exception.
    Under `python -O`, assert statements are silently removed,
    making this validation disappear.
    """
    start, stop, step = fh.unwrap_arange_args(start_or_stop=1, stop_=None, step_=2)
    assert start is None
    assert stop == 1
    assert step == 2

    with pytest.raises(RuntimeError, match="unwrap_arange_args"):
        fh.unwrap_arange_args(start_or_stop=None, stop_=None, step_=1)


def test_wrap_arange_args_stop_none_raises_runtimeerror_not_asserts():
    """
    Regression: production assert replaced with explicit exception.
    Under `python -O`, assert statements are silently removed,
    making this validation disappear.
    """
    args, kwargs = fh.wrap_arange_args(
        start=None,
        stop=5 * u.m,
        step=1,
        expected_out_unit=fh.UNIT_FROM_LIKE_ARG,
    )
    assert args == (np.float64(5.0),)
    assert kwargs == {}

    with pytest.raises(RuntimeError, match="wrap_arange_args"):
        fh.wrap_arange_args(
            start=None,
            stop=None,
            step=1,
            expected_out_unit=fh.UNIT_FROM_LIKE_ARG,
        )


def test_wrap_arange_args_expected_out_unit_mismatch_raises_runtimeerror_not_asserts(
    monkeypatch,
):
    """
    Regression: production assert replaced with explicit exception.
    Under `python -O`, assert statements are silently removed,
    making this validation disappear.
    """

    def fake_quantities2arrays(*args):
        return args, u.m

    monkeypatch.setattr(fh, "_quantities2arrays", fake_quantities2arrays)

    with pytest.raises(RuntimeError, match="Unit mismatch after conversion"):
        fh.wrap_arange_args(
            start=None,
            stop=5 * u.m,
            step=1,
            expected_out_unit=u.s,
        )


def test_wrap_arange_args_stop_unit_mismatch_raises_runtimeerror_not_asserts(
    monkeypatch,
):
    """
    Regression: production assert replaced with explicit exception.
    Under `python -O`, assert statements are silently removed,
    making this validation disappear.
    """

    def fake_quantities2arrays(*args):
        return args, u.s

    monkeypatch.setattr(fh, "_quantities2arrays", fake_quantities2arrays)

    with pytest.raises(RuntimeError, match="does not match stop.unit"):
        fh.wrap_arange_args(
            start=None,
            stop=5 * u.m,
            step=1,
            expected_out_unit=fh.UNIT_FROM_LIKE_ARG,
        )
