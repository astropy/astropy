# The purpose of these tests are to ensure that calling quantities using
# array methods returns quantities with the right units, or raises exceptions.

import numpy as np
from numpy.testing.utils import assert_allclose

from ... import units as u
from ...tests.helper import pytest

NUMPY_LT_1P7 = [int(x) for x in np.__version__.split('.')[:2]] < [1, 7]


class TestQuantityArrayCopy(object):
    """
    Test whether arrays are properly copied/used in place
    """

    def test_copy_on_creation(self):
        v = np.arange(1000.)
        q_nocopy = u.Quantity(v, "km/s", copy=False)
        q_copy = u.Quantity(v, "km/s", copy=True)
        v[0] = -1.
        assert q_nocopy[0].value == v[0]
        assert q_copy[0].value != v[0]

    def test_to_copies(self):
        q = u.Quantity(np.arange(1.,100.), "km/s")
        q2 = q.to(u.m/u.s)
        assert np.all(q.value != q2.value)
        q3 = q.to(u.km/u.s)
        assert np.all(q.value == q3.value)
        q[0] = -1.*u.km/u.s
        assert q[0].value != q3[0].value

    def test_si_copies(self):
        q = u.Quantity(np.arange(100.), "m/s")
        q2 = q.si
        assert np.all(q.value == q2.value)
        q[0] = -1.*u.m/u.s
        assert q[0].value != q2[0].value

    def test_getitem_is_view(self):
        """Check that [keys] work, and that, like ndarray, it returns
        a view, so that changing one changes the other.

        Also test that one can add axes (closes #1422)
        """
        q = u.Quantity(np.arange(100.), "m/s")
        q_sel = q[10:20]
        q_sel[0] = -1.*u.m/u.s
        assert q_sel[0] == q[10]
        # also check that getitem can do new axes
        q2 = q[:, np.newaxis]
        q2[10,0] = -9*u.m/u.s
        assert np.all(q2.flatten() == q)


class TestQuantityStatsFuncs(object):
    """
    Test statistical functions
    """

    def test_mean(self):
        q1 = np.array([1., 2., 4., 5., 6.]) * u.m
        assert np.mean(q1) == 3.6 * u.m

    def test_mean_inplace(self):
        q1 = np.array([1., 2., 4., 5., 6.]) * u.m
        qi = 1.5 * u.s
        np.mean(q1, out=qi)
        assert qi == 3.6 * u.m

    def test_std(self):
        q1 = np.array([1., 2.]) * u.m
        assert np.std(q1) == 0.5 * u.m

    def test_std_inplace(self):

        # can't use decorator since test causes a segfault in Numpy < 1.7, and
        # py.test will run the test anyway to see if it works
        pytest.xfail()

        q1 = np.array([1., 2.]) * u.m
        qi = 1.5 * u.s
        np.std(q1, out=qi)
        assert qi == 0.5 * u.m

    def test_var(self):
        q1 = np.array([1., 2.]) * u.m
        assert np.var(q1) == 0.25 * u.m ** 2

    def test_var_inplace(self):

        # can't use decorator since test causes a segfault in Numpy < 1.7, and
        # py.test will run the test anyway to see if it works
        if NUMPY_LT_1P7:
            pytest.xfail()

        q1 = np.array([1., 2.]) * u.m
        qi = 1.5 * u.s
        np.var(q1, out=qi)
        assert qi == 0.25 * u.m ** 2

    def test_median(self):
        q1 = np.array([1., 2., 4., 5., 6.]) * u.m
        assert np.median(q1) == 4. * u.m

    def test_median_inplace(self):
        q1 = np.array([1., 2., 4., 5., 6.]) * u.m
        qi = 1.5 * u.s
        np.median(q1, out=qi)
        assert qi == 4 * u.m

    def test_min(self):
        q1 = np.array([1., 2., 4., 5., 6.]) * u.m
        assert np.min(q1) == 1. * u.m

    def test_min_inplace(self):
        q1 = np.array([1., 2., 4., 5., 6.]) * u.m
        qi = 1.5 * u.s
        np.min(q1, out=qi)
        assert qi == 1. * u.m

    def test_argmin(self):
        q1 = np.array([6., 2., 4., 5., 6.]) * u.m
        assert np.argmin(q1) == 1

    def test_max(self):
        q1 = np.array([1., 2., 4., 5., 6.]) * u.m
        assert np.max(q1) == 6. * u.m

    def test_max_inplace(self):
        q1 = np.array([1., 2., 4., 5., 6.]) * u.m
        qi = 1.5 * u.s
        np.max(q1, out=qi)
        assert qi == 6. * u.m

    def test_argmax(self):
        q1 = np.array([5., 2., 4., 5., 6.]) * u.m
        assert np.argmax(q1) == 4

    def test_clip(self):
        q1 = np.array([1., 2., 4., 5., 6.]) * u.km / u.m
        c1 = q1.clip(1500, 5.5 * u.Mm / u.km)
        assert np.all(c1 == np.array([1.5, 2., 4., 5., 5.5]) * u.km / u.m)

    def test_clip_inplace(self):
        q1 = np.array([1., 2., 4., 5., 6.]) * u.km / u.m
        c1 = q1.clip(1500, 5.5 * u.Mm / u.km, out=q1)
        assert np.all(q1 == np.array([1.5, 2., 4., 5., 5.5]) * u.km / u.m)
        c1[0] = 10 * u.Mm/u.mm
        assert np.all(c1.value == q1.value)

    def test_conj(self):
        q1 = np.array([1., 2., 4., 5., 6.]) * u.km / u.m
        assert np.all(q1.conj() == q1)

    def test_ptp(self):
        q1 = np.array([1., 2., 4., 5., 6.]) * u.m
        assert np.ptp(q1) == 5. * u.m

    @pytest.mark.xfail
    def test_ptp_inplace(self):
        q1 = np.array([1., 2., 4., 5., 6.]) * u.m
        qi = 1.5 * u.s
        np.ptp(q1, out=qi)
        assert qi == 5. * u.m

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

    def test_sum_inplace(self):
        q1 = np.array([1., 2., 6.]) * u.m
        qi = 1.5 * u.s
        np.sum(q1, out=qi)
        assert qi == 9. * u.m

    def test_cumsum(self):

        q1 = np.array([1, 2, 6]) * u.m
        assert np.all(q1.cumsum() == np.array([1, 3, 9]) * u.m)
        assert np.all(np.cumsum(q1) == np.array([1, 3, 9]) * u.m)

        q2 = np.array([4, 5, 9]) * u.s
        assert np.all(q2.cumsum() == np.array([4, 9, 18]) * u.s)
        assert np.all(np.cumsum(q2) == np.array([4, 9, 18]) * u.s)

    def test_cumsum_inplace(self):
        q1 = np.array([1, 2, 6]) * u.m
        qi = np.ones(3) * u.s
        np.cumsum(q1, out=qi)
        assert np.all(qi == np.array([1, 3, 9]) * u.m)
        q2 = q1
        q1.cumsum(out=q1)
        assert np.all(q2 == qi)

    def test_nansum(self):

        q1 = np.array([1., 2., np.nan]) * u.m
        assert np.all(q1.nansum() == 3. * u.m)
        assert np.all(np.nansum(q1) == 3. * u.m)

        q2 = np.array([[np.nan, 5., 9.], [1., np.nan, 1.]]) * u.s
        assert np.all(q2.nansum(0) == np.array([1., 5., 10.]) * u.s)
        assert np.all(np.nansum(q2, 0) == np.array([1., 5., 10.]) * u.s)

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

    def test_diff(self):

        q1 = np.array([1., 2., 4., 10.]) * u.m
        assert np.all(q1.diff() == np.array([1., 2., 6.]) * u.m)
        assert np.all(np.diff(q1) == np.array([1., 2., 6.]) * u.m)

    def test_ediff1d(self):

        q1 = np.array([1., 2., 4., 10.]) * u.m
        assert np.all(q1.ediff1d() == np.array([1., 2., 6.]) * u.m)
        assert np.all(np.ediff1d(q1) == np.array([1., 2., 6.]) * u.m)

    @pytest.mark.xfail
    def test_dot_func(self):

        q1 = np.array([1., 2., 4., 10.]) * u.m
        q2 = np.array([3., 4., 5., 6.]) * u.s
        q3 = np.dot(q1, q2)
        assert q3.value == np.dot(q1.value, q2.value)
        assert q3.unit == u.m * u.s

    def test_dot_meth(self):

        q1 = np.array([1., 2., 4., 10.]) * u.m
        q2 = np.array([3., 4., 5., 6.]) * u.s
        q3 = q1.dot(q2)
        assert q3.value == np.dot(q1.value, q2.value)
        assert q3.unit == u.m * u.s

    @pytest.mark.xfail
    def test_trace_func(self):

        q = np.array([[1.,2.],[3.,4.]]) * u.m
        assert np.trace(q) == 5. * u.m

    def test_trace_meth(self):

        q1 = np.array([[1.,2.],[3.,4.]]) * u.m
        assert q1.trace() == 5. * u.m

        cont = u.Quantity(4., u.s)

        q2 = np.array([[3.,4.],[5.,6.]]) * u.m
        q2.trace(out=cont)
        assert cont == 9. * u.m

    def test_clip_func(self):

        q = np.arange(10) * u.m
        assert np.all(np.clip(q, 3 * u.m, 6 * u.m) == np.array([3., 3.,3.,3.,4.,5.,6.,6.,6.,6.]) * u.m)

    def test_clip_meth(self):

        expected = np.array([3.,3.,3.,3.,4.,5.,6.,6.,6.,6.]) * u.m

        q1 = np.arange(10) * u.m
        q3 = q1.clip(3 * u.m, 6 * u.m)
        assert np.all(q1.clip(3 * u.m, 6 * u.m) == expected)

        cont = np.zeros(10) * u.s

        q1.clip(3 * u.m, 6 * u.m, out=cont)

        assert np.all(cont == expected)


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
        with pytest.raises(u.UnitsError):
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
        q1 = np.array([1, 2, 3], dtype=np.int64) * u.m / u.km
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
