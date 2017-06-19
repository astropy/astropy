from ... import units as u

from ..helper import assert_quantity_allclose, pytest


def test_assert_quantity_allclose():

    assert_quantity_allclose([1, 2], [1, 2])

    assert_quantity_allclose([1, 2] * u.m, [100, 200] * u.cm)

    assert_quantity_allclose([1, 2] * u.m, [101, 201] * u.cm, atol=2 * u.cm)

    with pytest.raises(AssertionError):
        assert_quantity_allclose([1, 2] * u.m, [90, 200] * u.cm)

    with pytest.raises(AssertionError):
        assert_quantity_allclose([1, 2] * u.m, [101, 201] * u.cm, atol=0.5 * u.cm)

    with pytest.raises(u.UnitsError) as exc:
        assert_quantity_allclose([1, 2] * u.m, [100, 200])
    assert exc.value.args[0] == "Units for 'desired' () and 'actual' (m) are not convertible"

    with pytest.raises(u.UnitsError) as exc:
        assert_quantity_allclose([1, 2], [100, 200] * u.cm)
    assert exc.value.args[0] == "Units for 'desired' (cm) and 'actual' () are not convertible"

    with pytest.raises(u.UnitsError) as exc:
        assert_quantity_allclose([1, 2] * u.m, [100, 200] * u.cm, atol=0.3)
    assert exc.value.args[0] == "Units for 'atol' () and 'actual' (m) are not convertible"

    with pytest.raises(u.UnitsError) as exc:
        assert_quantity_allclose([1, 2], [1, 2], atol=0.3 * u.m)
    assert exc.value.args[0] == "Units for 'atol' (m) and 'actual' () are not convertible"

    with pytest.raises(u.UnitsError) as exc:
        assert_quantity_allclose([1, 2], [1, 2], rtol=0.3 * u.m)
    assert exc.value.args[0] == "`rtol` should be dimensionless"
