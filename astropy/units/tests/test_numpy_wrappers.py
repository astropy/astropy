from ...tests.helper import pytest
from ... import units as u

from ..numpy_wrappers import assert_allclose, linspace


def test_assert_allclose():

    assert_allclose([1,2], [1,2])

    assert_allclose([1,2] * u.m, [100,200] * u.cm)

    with pytest.raises(AssertionError):
        assert_allclose([1,2] * u.m, [90,200] * u.cm)

    with pytest.raises(TypeError) as exc:
        assert_allclose([1,2] * u.m, [100,200])
    assert exc.value.args[0] == "If `actual` is a Quantity, `desired` should also be a Quantity"

    with pytest.raises(TypeError) as exc:
        assert_allclose([1,2], [100,200] * u.cm)
    assert exc.value.args[0] == "If `desired` is a Quantity, `actual` should also be a Quantity"


def test_linspace():

    assert_allclose(linspace(1, 10, 10), [1,2,3,4,5,6,7,8,9,10])

    assert_allclose(linspace(1 * u.m, 10 * u.m, 10), [1,2,3,4,5,6,7,8,9,10] * u.m)

    assert_allclose(linspace(100 * u.cm, 10 * u.m, 10), [1,2,3,4,5,6,7,8,9,10] * u.m)

    with pytest.raises(TypeError) as exc:
        linspace(1 * u.m, 10, 10)
    assert exc.value.args[0] == "If `start` is a Quantity, `stop` should also be a Quantity"

    with pytest.raises(TypeError) as exc:
        linspace(1, 10 * u.m, 10)
    assert exc.value.args[0] == "If `stop` is a Quantity, `start` should also be a Quantity"
