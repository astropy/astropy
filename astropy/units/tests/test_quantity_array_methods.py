# The purpose of these tests are to ensure that calling quantities using
# array methods returns quantities with the right units, or raises exceptions.

import numpy as np

from ... import units as u
from ...tests.helper import pytest
from ...tests.compat import assert_allclose


class TestQuantityStatsFuncs(object):
    """
    Test statistical functions
    """

    def test_mean(self):
        q1 = np.array([1., 2., 4., 5., 6.]) * u.m
        assert np.mean(q1) == 3.6 * u.m
        assert np.mean(q1).unit == u.m
        assert np.mean(q1).value == 3.6

    def test_std(self):
        q1 = np.array([1., 2.]) * u.m
        assert np.std(q1) == 0.5 * u.m
        assert np.std(q1).unit == u.m
        assert np.std(q1).value == 0.5

    def test_var(self):
        q1 = np.array([1., 2.]) * u.m
        assert np.var(q1) == 0.25 * u.m ** 2
        assert np.var(q1).unit == u.m ** 2
        assert np.var(q1).value == 0.25

    def test_median(self):
        q1 = np.array([1., 2., 4., 5., 6.]) * u.m
        assert np.median(q1) == 4. * u.m
        assert np.median(q1).unit == u.m
        assert np.median(q1).value == 4.

    def test_min(self):
        q1 = np.array([1., 2., 4., 5., 6.]) * u.m
        assert np.min(q1) == 1. * u.m
        assert np.min(q1).unit == u.m
        assert np.min(q1).value == 1.

    def test_max(self):
        q1 = np.array([1., 2., 4., 5., 6.]) * u.m
        assert np.max(q1) == 6. * u.m
        assert np.max(q1).unit == u.m
        assert np.max(q1).value == 6.

    def test_clip(self):
        q1 = np.array([1., 2., 4., 5., 6.]) * u.km / u.m
        c1 = q1.clip(1500, 5.5 * u.Mm / u.km)
        assert all(c1 == np.array([1.5, 2., 4., 5., 5.5]) * u.km / u.m)

    def test_conj(self):
        q1 = np.array([1., 2., 4., 5., 6.]) * u.km / u.m
        assert all(q1.conj() == q1)

    def test_ptp(self):
        q1 = np.array([1., 2., 4., 5., 6.]) * u.m
        assert np.ptp(q1) == 5. * u.m
        assert np.ptp(q1).unit == u.m
        assert np.ptp(q1).value == 5.

    def test_round(self):
        q1 = np.array([1.2, 2.2, 3.2]) * u.kg
        assert np.all(np.round(q1) == np.array([1, 2, 3]) * u.kg)

    def test_sum(self):

        q1 = np.array([1., 2., 6.]) * u.m
        assert np.all(q1.sum() == 9. * u.m)
        assert np.all(np.sum(q1) == 9. * u.m)

        q2 = np.array([[4., 5., 9.], [1., 1., 1.]]) * u.s
        assert np.all(q2.sum(0) == np.array([5., 6., 10.]) * u.s)
        assert np.all(np.sum(q2, 0) == np.array([5., 6., 10.]) * u.s)

    def test_cumsum(self):

        q1 = np.array([1, 2, 6]) * u.m
        assert np.all(q1.cumsum() == np.array([1, 3, 9]) * u.m)
        assert np.all(np.cumsum(q1) == np.array([1, 3, 9]) * u.m)

        q2 = np.array([4, 5, 9]) * u.s
        assert np.all(q2.cumsum() == np.array([4, 9, 18]) * u.s)
        assert np.all(np.cumsum(q2) == np.array([4, 9, 18]) * u.s)

    def test_prod(self):

        q1 = np.array([1, 2, 6]) * u.m
        with pytest.raises(ValueError) as exc:
            q1.prod()
        assert 'cannot use prod' in exc.value.args[0]

        with pytest.raises(ValueError) as exc:
            np.prod(q1)
        assert 'cannot use prod' in exc.value.args[0]

        q2 = np.array([3., 4., 5.]) * u.Unit(1)
        assert q2.prod() == 60. * u.Unit(1)
        assert np.prod(q2) == 60. * u.Unit(1)

    def test_cumprod(self):

        q1 = np.array([1, 2, 6]) * u.m
        with pytest.raises(ValueError) as exc:
            q1.cumprod()
        assert 'cannot use cumprod' in exc.value.args[0]

        with pytest.raises(ValueError) as exc:
            np.cumprod(q1)
        assert 'cannot use cumprod' in exc.value.args[0]

        q2 = np.array([3, 4, 5]) * u.Unit(1)
        assert np.all(q2.cumprod() == np.array([3, 12, 60]) * u.Unit(1))
        assert np.all(np.cumprod(q2) == np.array([3, 12, 60]) * u.Unit(1))


class TestArrayConversion(object):
    """
    Test array conversion methods
    """

    def test_item(self):
        q1 = np.array([1, 2, 3]) * u.m / u.km
        assert q1.item(1) == 2 * q1.unit
        q1.itemset(1, 1)
        assert q1.item(1) == 1000 * u.m / u.km
        q1.itemset(1, 100 * u.cm / u.km)
        assert q1.item(1) == 1 * u.m / u.km
        with pytest.raises(TypeError):
            q1.itemset(1, 1.5 * u.m / u.km)
        with pytest.raises(ValueError):
            q1.itemset()

        q1[1] = 1
        assert q1[1] == 1000 * u.m / u.km
        q1[1] = 100 * u.cm / u.km
        assert q1[1] == 1 * u.m / u.km
        with pytest.raises(TypeError):
            q1[1] = 1.5 * u.m / u.km

        q1 = np.array([1, 2, 3]) * u.m / u.km
        assert all(q1.take((0, 2)) == np.array([1, 3]) * u.m / u.km)
        q1.put((1, 2), (3, 4))
        assert np.all(q1.take((1, 2)) == np.array([3000, 4000]) * q1.unit)
        q1.put(0, 500 * u.cm / u.km)
        assert q1.item(0) == 5 * u.m / u.km

    def test_slice(self):
        q2 = np.array([[1, 2, 3], [4, 5, 6]]) * u.km / u.m
        q1 = q2.copy()
        q2[0, 0] = 10000
        assert q2.unit == q1.unit
        assert q2[0, 0].value == 10
        q2[0] = 9 * u.Mm / u.km
        assert all(q2.flatten()[:3].value == np.array([9, 9, 9]))
        q2[0, :-1] = 8000
        assert all(q2.flatten()[:3].value == np.array([8, 8, 9]))
        with pytest.raises(u.UnitsException):
            q2[1, 1] = 10 * u.s
        with pytest.raises(TypeError):
            q2[0, 1] = 1.5 * u.km / u.m

    def test_fill(self):
        q1 = np.array([1, 2, 3]) * u.m / u.km
        q1.fill(2)
        assert np.all(q1 == 2000 * u.m / u.km)

    def test_repeat_compress_diagonal(self):
        q1 = np.array([1, 2, 3]) * u.m / u.km
        q2 = q1.repeat(2)
        assert q2.unit == q1.unit
        assert all(q2.value == q1.value.repeat(2))
        q2.sort()
        assert q2.unit == q1.unit
        q2 = q1.compress(np.array([True, True, False, False]))
        assert q2.unit == q1.unit
        assert all(q2.value == q1.value.compress(np.array([True, True,
                                                           False, False])))
        q1 = np.array([[1, 2], [3, 4]]) * u.m / u.km
        q2 = q1.diagonal()
        assert q2.unit == q1.unit
        assert all(q2.value == q1.value.diagonal())

    def test_byte_type_view_field_changes(self):
        q1 = np.array([1, 2, 3]) * u.m / u.km
        q2 = q1.byteswap()
        assert q2.unit == q1.unit
        assert all(q2.value == q1.value.byteswap())
        q2 = q1.astype(np.float64)
        assert all(q2 == q1)
        assert q2.dtype == np.float64
        q2 = q1.view(np.ndarray)
        assert not hasattr(q2, 'unit')
        q2a = q1.getfield(np.int32, offset=0)
        q2b = q1.byteswap().getfield(np.int32, offset=4)
        assert q2a.unit == q1.unit
        assert all(q2b.byteswap() == q2a)

    def test_sort(self):
        q1 = np.array([1., 5., 2., 4.]) * u.km / u.m
        i = q1.argsort()
        assert not hasattr(i, 'unit')
        q1.sort()
        i = q1.searchsorted([1500, 2500])
        assert not hasattr(i, 'unit')
        assert all(i == q1.to(
            u.dimensionless_unscaled).value.searchsorted([1500, 2500]))

    def test_not_implemented(self):
        q1 = np.array([1, 2, 3]) * u.m / u.km

        with pytest.raises(NotImplementedError):
            q1.choose([0, 0, 1])

        with pytest.raises(NotImplementedError):
            q1.list()
        with pytest.raises(NotImplementedError):
            q1.tostring()
        with pytest.raises(NotImplementedError):
            q1.tofile(0)
        with pytest.raises(NotImplementedError):
            q1.dump('a.a')
        with pytest.raises(NotImplementedError):
            q1.dumps()
