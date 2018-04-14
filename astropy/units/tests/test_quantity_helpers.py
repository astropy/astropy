"""
Test ``allclose`` and ``isclose``.
``allclose`` was ``quantity_allclose`` in ``astropy.tests.helper``.
"""
import numpy as np
import pytest

from ... import units as u


@pytest.mark.parametrize(
    ('a', 'b'),
    [([1, 2], [1, 2]),
     ([1, 2] * u.m, [100, 200] * u.cm),
     (1 * u.s, 1000 * u.ms)])
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
