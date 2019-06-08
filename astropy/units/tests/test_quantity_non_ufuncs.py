import inspect

import numpy as np
from numpy.testing import assert_array_equal

import pytest

from astropy import units as u
from astropy.units.quantity_helper.function_helpers import (
    ARRAY_FUNCTION_ENABLED)


NO_ARRAY_FUNCTION = not ARRAY_FUNCTION_ENABLED


# To get the functions that could be covered, we look for those that
# are wrapped.  Of course, this does not give a full list pre-1.17.
all_wrapped_functions = {name: f for name, f in np.__dict__.items()
                         if callable(f) and hasattr(f, '__wrapped__') and
                         f is not np.printoptions}


class CoverageMeta(type):
    """Meta class that tracks which functions are covered by tests.

    Assumes that a test is called 'test_<function_name>'.
    """
    covered = set()

    def __new__(mcls, name, bases, members):
        for k, v in members.items():
            if inspect.isfunction(v) and k.startswith('test'):
                f = k.replace('test_', '')
                if f in all_wrapped_functions:
                    mcls.covered.add(all_wrapped_functions[f])

        return super().__new__(mcls, name, bases, members)


class BasicTestSetup(metaclass=CoverageMeta):
    """Test setup for functions that should not change the unit.

    Also provides a default Quantity with shape (3, 3) and units of m.
    """
    def setup(self):
        self.q = np.arange(9.).reshape(3, 3) / 4. * u.m


class InvariantUnitTestSetup(BasicTestSetup):
    def check(self, func, *args, **kwargs):
        o = func(self.q, *args, **kwargs)
        expected = func(self.q.value, *args, **kwargs) * self.q.unit
        assert o.shape == expected.shape
        assert np.all(o == expected)


class NoUnitTestSetup(BasicTestSetup):
    def check(self, func, *args, **kwargs):
        out = func(self.q, *args, **kwargs)
        expected = func(self.q.value, *args, *kwargs)
        assert type(out) is type(expected)
        if isinstance(expected, tuple):
            assert all(np.all(o == x) for o, x in zip(out, expected))
        else:
            assert np.all(out == expected)


class TestShapeInformation(BasicTestSetup):
    def test_alen(self):
        assert np.alen(self.q) == 3

    def test_shape(self):
        assert np.shape(self.q) == (3, 3)

    def test_size(self):
        assert np.size(self.q) == 9

    def test_ndim(self):
        assert np.ndim(self.q) == 2


class TestShapeManipulation(InvariantUnitTestSetup):
    # Note: do not parametrize the below, since test names are used
    # to check coverage.
    def test_reshape(self):
        self.check(np.reshape, (9, 1))

    def test_ravel(self):
        self.check(np.ravel)

    def test_moveaxis(self):
        self.check(np.moveaxis, 0, 1)

    def test_rollaxis(self):
        self.check(np.rollaxis, 0, 2)

    def test_swapaxes(self):
        self.check(np.swapaxes, 0, 1)

    def test_transpose(self):
        self.check(np.transpose)

    def test_atleast_1d(self):
        q = 1. * u.m
        o, so = np.atleast_1d(q, self.q)
        assert o.shape == (1,)
        assert o == q
        expected = np.atleast_1d(self.q.value) * u.m
        assert np.all(so == expected)

    def test_atleast_2d(self):
        q = 1. * u.m
        o, so = np.atleast_2d(q, self.q)
        assert o.shape == (1, 1)
        assert o == q
        expected = np.atleast_2d(self.q.value) * u.m
        assert np.all(so == expected)

    def test_atleast_3d(self):
        q = 1. * u.m
        o, so = np.atleast_3d(q, self.q)
        assert o.shape == (1, 1, 1)
        assert o == q
        expected = np.atleast_3d(self.q.value) * u.m
        assert np.all(so == expected)

    def test_expand_dims(self):
        self.check(np.expand_dims, 1)

    def test_squeeze(self):
        o = np.squeeze(self.q[:, np.newaxis, :])
        assert o.shape == (3, 3)
        assert np.all(o == self.q)

    def test_flip(self):
        self.check(np.flip)

    def test_fliplr(self):
        self.check(np.fliplr)

    def test_flipud(self):
        self.check(np.flipud)

    def test_rot90(self):
        self.check(np.rot90)

    def test_broadcast_to(self):
        # TODO: should we change the default for subok?
        self.check(np.broadcast_to, (3, 3, 3), subok=True)

    def test_broadcast_arrays(self):
        # TODO: should we change the default for subok?
        q2 = np.ones((3, 3, 3)) / u.s
        o1, o2 = np.broadcast_arrays(self.q, q2, subok=True)
        assert isinstance(o1, u.Quantity)
        assert isinstance(o2, u.Quantity)
        assert o1.shape == o2.shape == (3, 3, 3)
        assert np.all(o1 == self.q)
        assert np.all(o2 == q2)


class TestArgFunctions(NoUnitTestSetup):
    def test_argmin(self):
        self.check(np.argmin)

    def test_argmax(self):
        self.check(np.argmax)

    def test_argsort(self):
        self.check(np.argsort)

    def test_lexsort(self):
        self.check(np.lexsort)

    def test_searchsorted(self):
        q = self.q.ravel()
        q2 = np.array([150., 350.]) * u.cm
        out = np.searchsorted(q, q2)
        expected = np.searchsorted(q.value, q2.to_value(q.unit))
        assert np.all(out == expected)

    def test_nonzero(self):
        self.check(np.nonzero)

    def test_argwhere(self):
        self.check(np.argwhere)

    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_argpartition(self):
        self.check(np.argpartition, 2)

    def test_flatnonzero(self):
        self.check(np.flatnonzero)


class TestAlongAxis(BasicTestSetup):
    def test_take_along_axis(self):
        indices = np.expand_dims(np.argmax(self.q, axis=0), axis=0)
        out = np.take_along_axis(self.q, indices, axis=0)
        expected = np.take_along_axis(self.q.value, indices,
                                      axis=0) * self.q.unit
        assert np.all(out == expected)

    def test_put_along_axis(self):
        q = self.q.copy()
        indices = np.expand_dims(np.argmax(self.q, axis=0), axis=0)
        np.put_along_axis(q, indices, axis=0, values=-100 * u.cm)
        expected = q.value.copy()
        np.put_along_axis(expected, indices, axis=0, values=-1)
        expected = expected * q.unit
        assert np.all(q == expected)


class TestIndicesFrom(NoUnitTestSetup):
    def test_diag_indices_from(self):
        self.check(np.diag_indices_from)

    def test_triu_indices_from(self):
        self.check(np.triu_indices_from)

    def test_tril_indices_from(self):
        self.check(np.tril_indices_from)


class TestRealImag(InvariantUnitTestSetup):
    def setup(self):
        self.q = (np.arange(9.).reshape(3, 3) + 1j) * u.m

    def test_real(self):
        self.check(np.real)

    def test_imag(self):
        self.check(np.imag)


class TestCopyAndCreation(InvariantUnitTestSetup):
    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_copy(self):
        self.check(np.copy)
        # Also as kwarg
        copy = np.copy(a=self.q)
        assert_array_equal(copy, self.q)

    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_asfarray(self):
        self.check(np.asfarray)
        farray = np.asfarray(a=self.q)
        assert_array_equal(farray, self.q)

    def test_empty_like(self):
        o = np.empty_like(self.q)
        assert o.shape == (3, 3)
        assert isinstance(o, u.Quantity)
        assert o.unit == self.q.unit
        o2 = np.empty_like(prototype=self.q)
        assert o2.shape == (3, 3)
        assert isinstance(o2, u.Quantity)
        assert o2.unit == self.q.unit
        o3 = np.empty_like(self.q, subok=False)
        assert type(o3) is np.ndarray

    def test_zeros_like(self):
        self.check(np.zeros_like)
        o2 = np.zeros_like(a=self.q)
        assert_array_equal(o2, self.q * 0.)

    def test_ones_like(self):
        self.check(np.ones_like)

    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_full_like(self):
        o = np.full_like(self.q, 0.5 * u.km)
        expected = np.empty_like(self.q.value) * u.m
        expected[...] = 0.5 * u.km
        assert np.all(o == expected)
        with pytest.raises(u.UnitsError):
            np.full_like(self.q, 0.5 * u.s)


class TestAccessingParts(InvariantUnitTestSetup):
    def test_diag(self):
        self.check(np.diag)

    def test_diagonal(self):
        self.check(np.diagonal)

    def test_diagflat(self):
        self.check(np.diagflat)

    def test_compress(self):
        o = np.compress([True, False, True], self.q, axis=0)
        expected = np.compress([True, False, True], self.q.value,
                               axis=0) * self.q.unit
        assert np.all(o == expected)

    def test_extract(self):
        o = np.extract([True, False, True], self.q)
        expected = np.extract([True, False, True],
                              self.q.value) * self.q.unit
        assert np.all(o == expected)

    def test_delete(self):
        self.check(np.delete, slice(1, 2), 0)
        self.check(np.delete, [0, 2], 1)

    def test_trim_zeros(self):
        q = self.q.ravel()
        out = np.trim_zeros(q)
        expected = np.trim_zeros(q.value) * u.m
        assert np.all(out == expected)

    def test_roll(self):
        self.check(np.roll, 1)
        self.check(np.roll, 1, axis=0)

    def test_take(self):
        self.check(np.take, [0, 1], axis=1)
        self.check(np.take, 1)


class TestSettingParts(metaclass=CoverageMeta):
    def test_put(self):
        q = np.arange(3.) * u.m
        np.put(q, [0, 2], [50, 150] * u.cm)
        assert q.unit == u.m
        expected = [50, 100, 150] * u.cm
        assert np.all(q == expected)

    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_putmask(self):
        q = np.arange(3.) * u.m
        mask = [True, False, True]
        values = [50, 0, 150] * u.cm
        np.putmask(q, mask, values)
        assert q.unit == u.m
        expected = [50, 100, 150] * u.cm
        assert np.all(q == expected)
        with pytest.raises(u.UnitsError):
            np.putmask(q, mask, values.value)
        with pytest.raises(u.UnitsError):
            np.putmask(q.value, mask, values)
        a = np.arange(3.)
        values = [50, 0, 150] * u.percent
        np.putmask(a, mask, values)
        expected = np.array([0.5, 1., 1.5])
        assert np.all(a == expected)

    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_place(self):
        q = np.arange(3.) * u.m
        np.place(q, [True, False, True], [50, 150] * u.cm)
        assert q.unit == u.m
        expected = [50, 100, 150] * u.cm
        assert np.all(q == expected)

    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_copyto(self):
        q = np.arange(3.) * u.m
        np.copyto(q, [50, 0, 150] * u.cm, where=[True, False, True])
        assert q.unit == u.m
        expected = [50, 100, 150] * u.cm
        assert np.all(q == expected)

    def test_fill_diagonal(self):
        q = np.arange(9.).reshape(3, 3) * u.m
        expected = q.value.copy()
        np.fill_diagonal(expected, 0.25)
        expected = expected * u.m
        np.fill_diagonal(q, 25. * u.cm)
        assert q.unit == u.m
        assert np.all(q == expected)


class TestRepeat(InvariantUnitTestSetup):
    def test_tile(self):
        self.check(np.tile, 2)

    def test_repeat(self):
        self.check(np.repeat, 2)

    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_resize(self):
        self.check(np.resize, (4, 4))


class TestConcatenate(metaclass=CoverageMeta):
    def setup(self):
        self.q1 = np.arange(6.).reshape(2, 3) * u.m
        self.q2 = self.q1.to(u.cm)

    def check(self, func, *args, **kwargs):
        q_list = kwargs.pop('q_list', [self.q1, self.q2])
        o = func(q_list, *args, **kwargs)
        unit = q_list[0].unit
        v_list = [q.to_value(unit) for q in q_list]
        expected = func(v_list, *args, **kwargs) * unit
        assert o.shape == expected.shape
        assert np.all(o == expected)

    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_concatenate(self):
        self.check(np.concatenate)

    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_stack(self):
        self.check(np.stack)

    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_column_stack(self):
        self.check(np.column_stack)

    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_hstack(self):
        self.check(np.hstack)

    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_vstack(self):
        self.check(np.vstack)

    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_dstack(self):
        self.check(np.dstack)

    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_block(self):
        self.check(np.block)

    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_append(self):
        out = np.append(self.q1, self.q2, axis=0)
        assert out.unit == self.q1.unit
        expected = np.append(self.q1.value, self.q2.to_value(self.q1.unit),
                             axis=0) * self.q1.unit
        assert np.all(out == expected)

    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_insert(self):
        # Unit of inserted values is ignored.
        q = np.arange(12.).reshape(6, 2) * u.m
        out = np.insert(q, (3, 5), [50., 25.] * u.cm)
        assert isinstance(out, u.Quantity)
        assert out.unit == q.unit
        expected = np.insert(q.value, (3, 5), [0.5, 0.25]) * u.m
        assert np.all(out == expected)

        with pytest.raises(TypeError):
            np.insert(q, 3 * u.cm, 50. * u.cm)

    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_pad(self):
        q = np.arange(1., 6.) * u.m
        out = np.pad(q, (2, 3), 'constant', constant_values=(0., 150.*u.cm))
        assert out.unit == q.unit
        expected = np.pad(q.value, (2, 3), 'constant',
                          constant_values=(0., 1.5)) * q.unit
        assert np.all(out == expected)
        out = np.pad(q, (2, 3), 'linear_ramp', end_values=(25.*u.cm, 0.))
        assert out.unit == q.unit
        expected = np.pad(q.value, (2, 3), 'linear_ramp',
                          end_values=(0.25, 0.)) * q.unit
        assert np.all(out == expected)


class TestSplit(metaclass=CoverageMeta):
    def setup(self):
        self.q = np.arange(54.).reshape(3, 3, 6) * u.m

    def check(self, func, *args, **kwargs):
        out = func(self.q, *args, **kwargs)
        expected = func(self.q.value, *args, **kwargs)
        expected = [x * self.q.unit for x in expected]
        assert len(out) == len(expected)
        assert all(o.shape == x.shape for o, x in zip(out, expected))
        assert all(np.all(o == x) for o, x in zip(out, expected))

    def test_split(self):
        self.check(np.split, [1])

    def test_array_split(self):
        self.check(np.array_split, 2)

    def test_hsplit(self):
        self.check(np.hsplit, [1, 4])

    def test_vsplit(self):
        self.check(np.vsplit, [1])

    def test_dsplit(self):
        self.check(np.dsplit, [1])


class TestUfuncReductions(InvariantUnitTestSetup):
    def test_amax(self):
        self.check(np.amax)

    def test_amin(self):
        self.check(np.amin)

    def test_sum(self):
        self.check(np.sum)

    def test_cumsum(self):
        self.check(np.cumsum)

    def test_any(self):
        with pytest.raises(NotImplementedError):
            np.any(self.q)

    def test_all(self):
        with pytest.raises(NotImplementedError):
            np.all(self.q)

    def test_sometrue(self):
        with pytest.raises(NotImplementedError):
            np.sometrue(self.q)

    def test_alltrue(self):
        with pytest.raises(NotImplementedError):
            np.alltrue(self.q)

    def test_prod(self):
        with pytest.raises(u.UnitsError):
            np.prod(self.q)

    def test_product(self):
        with pytest.raises(u.UnitsError):
            np.product(self.q)

    def test_cumprod(self):
        with pytest.raises(u.UnitsError):
            np.cumprod(self.q)

    def test_cumproduct(self):
        with pytest.raises(u.UnitsError):
            np.cumproduct(self.q)


class TestUfuncLike(InvariantUnitTestSetup):
    def test_ptp(self):
        self.check(np.ptp)
        self.check(np.ptp, axis=0)

    def test_round_(self):
        self.check(np.round_)

    def test_around(self):
        self.check(np.around)

    def test_fix(self):
        self.check(np.fix)

    def test_angle(self):
        q = np.array([1+0j, 0+1j, 1+1j, 0+0j]) * u.m
        out = np.angle(q)
        expected = np.angle(q.value) * u.radian
        assert np.all(out == expected)

    def test_i0(self):
        q = np.array([0., 10., 20.]) * u.percent
        out = np.i0(q)
        expected = np.i0(q.to_value(u.one)) * u.one
        assert isinstance(out, u.Quantity)
        assert np.all(out == expected)
        with pytest.raises(u.UnitsError):
            np.i0(self.q)

    def test_clip(self):
        qmin = 200 * u.cm
        qmax = [270, 280, 290] * u.cm
        out = np.clip(self.q, qmin, qmax)
        expected = np.clip(self.q.value, qmin.to_value(self.q.unit),
                           qmax.to_value(self.q.unit)) * self.q.unit
        assert np.all(out == expected)

    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_sinc(self):
        q = [0., 3690., -270., 690.] * u.deg
        out = np.sinc(q)
        expected = np.sinc(q.to_value(u.radian)) * u.one
        assert isinstance(out, u.Quantity)
        assert np.all(out == expected)
        with pytest.raises(u.UnitsError):
            np.sinc(1.*u.one)

    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_where(self):
        out = np.where([True, False, True], self.q, 1. * u.km)
        expected = np.where([True, False, True], self.q.value,
                            1000.) * self.q.unit
        assert np.all(out == expected)

    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_choose(self):
        # from np.choose docstring
        a = np.array([0, 1]).reshape((2, 1, 1))
        q1 = np.array([1, 2, 3]).reshape((1, 3, 1)) * u.cm
        q2 = np.array([-1, -2, -3, -4, -5]).reshape((1, 1, 5)) * u.m
        out = np.choose(a, (q1, q2))
        # result is 2x3x5, res[0,:,:]=c1, res[1,:,:]=c2
        expected = np.choose(a, (q1.value, q2.to_value(q1.unit))) * u.cm
        assert np.all(out == expected)

    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_select(self):
        q = self.q
        out = np.select([q < 0.55 * u.m, q > 1. * u.m],
                        [q, q.to(u.cm)], default=-1. * u.km)
        expected = np.select([q.value < 0.55, q.value > 1],
                             [q.value, q.value], default=-1000) * u.m
        assert np.all(out == expected)

    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_real_if_close(self):
        q = np.array([1+0j, 0+1j, 1+1j, 0+0j]) * u.m
        out = np.real_if_close(q)
        expected = np.real_if_close(q.value) * u.m
        assert np.all(out == expected)

    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_tril(self):
        self.check(np.tril)

    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_triu(self):
        self.check(np.triu)

    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_unwrap(self):
        q = [0., 3690., -270., 690.] * u.deg
        out = np.unwrap(q)
        expected = np.rad2deg(np.unwrap(q.to_value(u.rad))) * u.deg
        assert np.allclose(out, expected, atol=1*u.urad, rtol=0)
        with pytest.raises(u.UnitsError):
            np.unwrap([1., 2.]*u.m)
        with pytest.raises(u.UnitsError):
            np.unwrap(q, discont=1.*u.m)

    def test_nan_to_num(self):
        q = np.array([-np.inf, +np.inf, np.nan, 3., 4.]) * u.m
        out = np.nan_to_num(q)
        expected = np.nan_to_num(q.value) * q.unit
        assert np.all(out == expected)

    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_nan_to_num_complex(self):
        q = np.array([-np.inf, +np.inf, np.nan, 3., 4.]) * u.m
        out = np.nan_to_num(q, nan=1.*u.km, posinf=2.*u.km, neginf=-2*u.km)
        expected = [-2000., 2000., 1000., 3., 4.] * u.m
        assert np.all(out == expected)


class TestUfuncLikeTests(metaclass=CoverageMeta):
    def setup(self):
        self.q = np.array([-np.inf, +np.inf, np.nan, 3., 4.]) * u.m

    def check(self, func):
        out = func(self.q)
        expected = func(self.q.value)
        assert type(out) is np.ndarray
        assert out.dtype.kind == 'b'
        assert np.all(out == expected)

    def test_isposinf(self):
        self.check(np.isposinf)

    def test_isneginf(self):
        self.check(np.isneginf)

    def test_isreal(self):
        self.check(np.isreal)
        assert not np.isreal([1. + 1j]*u.m)

    def test_iscomplex(self):
        self.check(np.iscomplex)
        assert np.iscomplex([1. + 1j]*u.m)

    def test_isclose(self):
        q1 = np.arange(3.) * u.m
        q2 = np.array([0., 102., 199.]) * u.cm
        atol = 1.5 * u.cm
        rtol = 1. * u.percent
        out = np.isclose(q1, q2, atol=atol)
        expected = np.isclose(q1.value, q2.to_value(q1.unit),
                              atol=atol.to_value(q1.unit))
        assert type(out) is np.ndarray
        assert out.dtype.kind == 'b'
        assert np.all(out == expected)
        out = np.isclose(q1, q2, atol=0, rtol=rtol)
        expected = np.isclose(q1.value, q2.to_value(q1.unit),
                              atol=0, rtol=0.01)
        assert type(out) is np.ndarray
        assert out.dtype.kind == 'b'
        assert np.all(out == expected)

    @pytest.mark.xfail
    def test_isclose_failure(self):
        q_cm = self.q.to(u.cm)
        # atol does not have units; TODO: should this work by default?
        out = np.isclose(self.q, q_cm)
        expected = np.isclose(self.q.value, q_cm.to_value(u.m))
        assert np.all(out == expected)


class TestReductionLikeFunctions(InvariantUnitTestSetup):
    def test_average(self):
        q1 = np.arange(9.).reshape(3, 3) * u.m
        q2 = np.eye(3) / u.s
        o = np.average(q1, weights=q2)
        expected = np.average(q1.value, weights=q2.value) * u.m
        assert np.all(o == expected)

    def test_mean(self):
        self.check(np.mean)

    def test_std(self):
        self.check(np.std)

    def test_var(self):
        o = np.var(self.q)
        expected = np.var(self.q.value) * self.q.unit ** 2
        assert np.all(o == expected)

    def test_median(self):
        self.check(np.median)

    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_quantile(self):
        self.check(np.quantile, 0.5)
        o = np.quantile(self.q, 50 * u.percent)
        expected = np.quantile(self.q.value, 0.5) * u.m
        assert np.all(o == expected)
        # For ndarray input, we return a Quantity.
        o2 = np.quantile(self.q.value, 50 * u.percent)
        assert o2.unit == u.dimensionless_unscaled
        assert np.all(o2 == expected.value)

    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_percentile(self):
        self.check(np.percentile, 0.5)
        o = np.percentile(self.q, 0.5 * u.one)
        expected = np.percentile(self.q.value, 50) * u.m
        assert np.all(o == expected)

    def test_trace(self):
        self.check(np.trace)

    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_count_nonzero(self):
        q1 = np.arange(9.).reshape(3, 3) * u.m
        o = np.count_nonzero(q1)
        assert type(o) is not u.Quantity
        assert o == 8
        o = np.count_nonzero(q1, axis=1)
        # Returns integer Quantity with units of m
        assert type(o) is np.ndarray
        assert np.all(o == np.array([2, 3, 3]))

    def test_allclose(self):
        q1 = np.arange(3.) * u.m
        q2 = np.array([0., 101., 199.]) * u.cm
        atol = 2 * u.cm
        rtol = 1. * u.percent
        assert np.allclose(q1, q2, atol=atol)
        assert np.allclose(q1, q2, atol=0., rtol=rtol)

    def test_allclose_failures(self):
        q1 = np.arange(3.) * u.m
        q2 = np.array([0., 101., 199.]) * u.cm
        with pytest.raises(u.UnitsError):
            # Default atol breaks code; TODO: should this work?
            assert np.allclose(q1, q2)
        with pytest.raises(u.UnitsError):
            np.allclose(q1, q2, atol=2, rtol=0)
        with pytest.raises(u.UnitsError):
            np.allclose(q1, q2, atol=0, rtol=1. * u.s)

    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_array_equal(self):
        q1 = np.arange(3.) * u.m
        q2 = q1.to(u.cm)
        assert np.array_equal(q1, q2)
        q3 = q1.value * u.cm
        assert not np.array_equal(q1, q3)

    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_array_equiv(self):
        q1 = np.array([[0., 1., 2.]]*3) * u.m
        q2 = q1[0].to(u.cm)
        assert np.array_equiv(q1, q2)
        q3 = q1[0].value * u.cm
        assert not np.array_equiv(q1, q3)


class TestNanFunctions(InvariantUnitTestSetup):
    def setup(self):
        super().setup()
        self.q[1, 1] = np.nan

    def test_nanmax(self):
        self.check(np.nanmax)

    def test_nanmin(self):
        self.check(np.nanmin)

    def test_nanargmin(self):
        out = np.nanargmin(self.q)
        expected = np.nanargmin(self.q.value)
        assert out == expected

    def test_nanargmax(self):
        out = np.nanargmax(self.q)
        expected = np.nanargmax(self.q.value)
        assert out == expected

    def test_nanmean(self):
        self.check(np.nanmean)

    def test_nanmedian(self):
        self.check(np.nanmedian)

    def test_nansum(self):
        self.check(np.nansum)

    def test_nancumsum(self):
        self.check(np.nancumsum)

    def test_nanstd(self):
        self.check(np.nanstd)

    def test_nanvar(self):
        out = np.nanvar(self.q)
        expected = np.nanvar(self.q.value) * self.q.unit ** 2
        assert np.all(out == expected)

    def test_nanprod(self):
        with pytest.raises(u.UnitsError):
            np.nanprod(self.q)

    def test_nancumprod(self):
        with pytest.raises(u.UnitsError):
            np.nancumprod(self.q)

    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_nanquantile(self):
        self.check(np.nanquantile, 0.5)
        o = np.nanquantile(self.q, 50 * u.percent)
        expected = np.nanquantile(self.q.value, 0.5) * u.m
        assert np.all(o == expected)

    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_nanpercentile(self):
        self.check(np.nanpercentile, 0.5)
        o = np.nanpercentile(self.q, 0.5 * u.one)
        expected = np.nanpercentile(self.q.value, 50) * u.m
        assert np.all(o == expected)


class TestVariousProductFunctions(metaclass=CoverageMeta):
    """
    Test functions that are similar to gufuncs
    """
    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_cross(self):
        q1 = np.arange(6.).reshape(2, 3) * u.m
        q2 = np.array([4., 5., 6.]) / u.s
        o = np.cross(q1, q2)
        expected = np.cross(q1.value, q2.value) * u.m / u.s
        assert np.all(o == expected)

    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_outer(self):
        q1 = np.array([1, 2, 3]) * u.m
        q2 = np.array([1, 2]) / u.s
        o = np.outer(q1, q2)
        assert np.all(o == np.array([[1, 2], [2, 4], [3, 6]]) * u.m / u.s)

    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_inner(self):
        q1 = np.array([1, 2, 3]) * u.m
        q2 = np.array([4, 5, 6]) / u.s
        o = np.inner(q1, q2)
        assert o == 32 * u.m / u.s

    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_dot(self):
        q1 = np.array([1., 2., 3.]) * u.m
        q2 = np.array([4., 5., 6.]) / u.s
        o = np.dot(q1, q2)
        assert o == 32. * u.m / u.s

    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_vdot(self):
        q1 = np.array([1j, 2j, 3j]) * u.m
        q2 = np.array([4j, 5j, 6j]) / u.s
        o = np.vdot(q1, q2)
        assert o == (32. + 0j) * u.m / u.s

    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_tensordot(self):
        # From the docstring example
        a = np.arange(60.).reshape(3, 4, 5) * u.m
        b = np.arange(24.).reshape(4, 3, 2) / u.s
        c = np.tensordot(a, b, axes=([1, 0], [0, 1]))
        expected = np.tensordot(a.value, b.value,
                                axes=([1, 0], [0, 1])) * u.m / u.s
        assert np.all(c == expected)

    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_kron(self):
        q1 = np.eye(2) * u.m
        q2 = np.ones(2) / u.s
        o = np.kron(q1, q2)
        expected = np.kron(q1.value, q2.value) * u.m / u.s
        assert np.all(o == expected)

    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_einsum(self):
        q1 = np.arange(9.).reshape(3, 3) * u.m
        o = np.einsum('...i', q1)
        assert np.all(o == q1)
        o = np.einsum('ii', q1)
        expected = np.einsum('ii', q1.value) * u.m
        assert np.all(o == expected)
        q2 = np.eye(3) / u.s
        o = np.einsum('ij,jk', q1, q2)
        assert np.all(o == q1 / u.s)

    def test_einsum_path(self):
        q1 = np.arange(9.).reshape(3, 3) * u.m
        o = np.einsum_path('...i', q1)
        assert o[0] == ['einsum_path', (0,)]
        o = np.einsum_path('ii', q1)
        assert o[0] == ['einsum_path', (0,)]
        q2 = np.eye(3) / u.s
        o = np.einsum_path('ij,jk', q1, q2)
        assert o[0] == ['einsum_path', (0, 1)]


class TestIntDiffFunctions(metaclass=CoverageMeta):
    def test_trapz(self):
        y = np.arange(9.) * u.m / u.s
        out = np.trapz(y)
        expected = np.trapz(y.value) * y.unit
        assert np.all(out == expected)

        dx = 10. * u.s
        out = np.trapz(y, dx=dx)
        expected = np.trapz(y.value, dx=dx.value) * y.unit * dx.unit
        assert np.all(out == expected)

        x = np.arange(9.) * u.s
        out = np.trapz(y, x)
        expected = np.trapz(y.value, x.value) * y.unit * x.unit
        assert np.all(out == expected)

    def test_diff(self):
        # Simple diff works out of the box.
        x = np.arange(10.) * u.m
        out = np.diff(x)
        expected = np.diff(x.value) * u.m
        assert np.all(out == expected)

    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_diff_prepend_append(self):
        x = np.arange(10.) * u.m
        out = np.diff(x, prepend=-12.5*u.cm, append=1*u.km)
        expected = np.diff(x.value, prepend=-0.125, append=1000.) * x.unit
        x = np.arange(10.) * u.m
        out = np.diff(x, prepend=-12.5*u.cm, append=1*u.km, n=2)
        expected = np.diff(x.value, prepend=-0.125, append=1000.,
                           n=2) * x.unit
        assert np.all(out == expected)

    def test_gradient(self):
        # Simple gradient works out of the box.
        x = np.arange(10.) * u.m
        out = np.gradient(x)
        expected = np.gradient(x.value) * u.m
        assert np.all(out == expected)

    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_gradient_spacing(self):
        # Simple gradient works out of the box.
        x = np.arange(10.) * u.m
        spacing = 10. * u.s
        out = np.gradient(x, spacing)
        expected = np.gradient(x.value, spacing.value) * (x.unit /
                                                          spacing.unit)
        assert np.all(out == expected)
        f = np.array([[1, 2, 6], [3, 4, 5]]) * u.m
        dx = 2. * u.s
        y = [1., 1.5, 3.5] * u.GHz
        dfdx, dfdy = np.gradient(f, dx, y)
        exp_dfdx, exp_dfdy = np.gradient(f.value, dx.value, y.value)
        exp_dfdx = exp_dfdx * f.unit / dx.unit
        exp_dfdy = exp_dfdy * f.unit / y.unit
        assert np.all(dfdx == exp_dfdx)
        assert np.all(dfdy == exp_dfdy)

        dfdx2 = np.gradient(f, dx, axis=0)
        assert np.all(dfdx2 == exp_dfdx)
        dfdy2 = np.gradient(f, y, axis=1)
        assert np.all(dfdy2 == exp_dfdy)


class TestSpaceFunctions(metaclass=CoverageMeta):
    def test_linspace(self):
        # Note: linspace gets unit of end point, not superlogical.
        out = np.linspace(1000.*u.m, 10.*u.km, 5)
        expected = np.linspace(1, 10, 5) * u.km
        assert np.all(out == expected)

        q1 = np.arange(6.).reshape(2, 3) * u.m
        q2 = 10000. * u.cm
        out = np.linspace(q1, q2, 5)
        expected = np.linspace(q1.to_value(q2.unit), q2.value, 5) * q2.unit
        assert np.all(out == expected)

    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_logspace(self):
        unit = u.m / u.s**2
        out = np.logspace(10.*u.dex(unit), 20*u.dex(unit), 10)
        expected = np.logspace(10., 20., 10) * unit
        assert np.all(out == expected)
        out = np.logspace(10.*u.STmag, 20*u.STmag, 10)
        expected = np.logspace(10., 20., 10, base=10.**(-0.4)) * u.ST
        assert u.allclose(out, expected)

    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_geomspace(self):
        out = np.geomspace(1000.*u.m, 10.*u.km, 5)
        expected = np.geomspace(1, 10, 5) * u.km
        assert np.all(out == expected)

        q1 = np.arange(1., 7.).reshape(2, 3) * u.m
        q2 = 10000. * u.cm
        out = np.geomspace(q1, q2, 5)
        expected = np.geomspace(q1.to_value(q2.unit), q2.value, 5) * q2.unit
        assert np.all(out == expected)


class TestInterpolationFunctions(metaclass=CoverageMeta):
    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_interp(self):
        x = np.array([1250., 2750.]) * u.m
        xp = np.arange(5.) * u.km
        yp = np.arange(5.) * u.day
        out = np.interp(x, xp, yp)
        expected = np.interp(x.to_value(xp.unit), xp.value, yp.value) * yp.unit
        assert np.all(out == expected)

    @pytest.mark.xfail
    def test_piecewise(self):
        # TODO: this needs a proper own implementation to take care of the
        # unit of the output of possible functions.
        x = np.linspace(-2.5, 2.5, 6) * u.m
        out = np.piecewise(x, [x < 0, x >= 0], [-1*u.s, 1*u.day])
        expected = np.piecewise(x.value, [x.value < 0, x.value >= 0],
                                [-1, 24*3600]) * u.s
        assert np.all(out == expected)


class TestHistogramFunctions(metaclass=CoverageMeta):
    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_bincount(self):
        i = np.array([1, 1, 2, 3, 2, 4])
        weights = np.arange(len(i)) * u.Jy
        out = np.bincount(i, weights)
        expected = np.bincount(i, weights.value) * weights.unit
        assert np.all(out == expected)

    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_digitize(self):
        x = np.array([1500., 2500., 4500.]) * u.m
        bins = np.arange(10.) * u.km
        out = np.digitize(x, bins)
        expected = np.digitize(x.to_value(bins.unit), bins.value)
        assert np.all(out == expected)

    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_histogram(self):
        x = np.array([1.1, 1.2, 1.3, 2.1, 5.1]) * u.m
        out_h, out_b = np.histogram(x)
        assert type(out_h) is np.ndarray
        assert isinstance(out_b, u.Quantity)
        expected_h, expected_b = np.histogram(x.value)
        expected_b = expected_b * x.unit
        assert np.all(out_h == expected_h)
        assert np.all(out_b == expected_b)
        # With bins
        out2_h, out2_b = np.histogram(x, [125, 200] * u.cm)
        assert type(out2_h) is np.ndarray
        assert isinstance(out2_b, u.Quantity)
        expected2_h, expected2_b = np.histogram(x.value, [1.25, 2.])
        expected2_b = expected2_b * x.unit
        assert np.all(out2_h == expected2_h)
        assert np.all(out2_b == expected2_b)
        # With density
        out2d_h, out2d_b = np.histogram(x, [125, 200] * u.cm, density=True)
        assert isinstance(out2d_h, u.Quantity)
        assert isinstance(out2d_b, u.Quantity)
        expected2d_h, expected2d_b = np.histogram(x.value, [1.25, 2.],
                                                  density=True)
        expected2d_h = expected2d_h / x.unit
        expected2d_b = expected2d_b * x.unit
        assert np.all(out2d_h == expected2d_h)
        assert np.all(out2d_b == expected2d_b)
        # With bins and weights
        weights = np.arange(len(x)) / u.s
        out3_h, out3_b = np.histogram(x, [125, 200] * u.cm, weights=weights)
        assert isinstance(out3_h, u.Quantity)
        assert isinstance(out3_b, u.Quantity)
        expected3_h, expected3_b = np.histogram(x.value, [1.25, 2.],
                                                weights=weights.value)
        expected3_h = expected3_h * weights.unit
        expected3_b = expected3_b * x.unit
        assert np.all(out3_h == expected3_h)
        assert np.all(out3_b == expected3_b)
        out3d_h, out3d_b = np.histogram(x, [125, 200] * u.cm, weights=weights,
                                        density=True)
        assert isinstance(out3_h, u.Quantity)
        assert isinstance(out3_b, u.Quantity)
        expected3d_h, expected3d_b = np.histogram(x.value, [1.25, 2.],
                                                  weights=weights.value,
                                                  density=True)
        expected3d_h = expected3d_h * weights.unit / x.unit
        expected3d_b = expected3d_b * x.unit
        assert np.all(out3_h == expected3_h)
        assert np.all(out3_b == expected3_b)

        with pytest.raises(u.UnitsError):
            np.histogram(x, [125, 200] * u.s)

        with pytest.raises(u.UnitsError):
            np.histogram(x, [125, 200])

        with pytest.raises(u.UnitsError):
            np.histogram(x.value, [125, 200] * u.s)

    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_histogram_bin_edges(self):
        x = np.array([1.1, 1.2, 1.3, 2.1, 5.1]) * u.m
        out_b = np.histogram_bin_edges(x)
        expected_b = np.histogram_bin_edges(x.value) * x.unit
        assert np.all(out_b == expected_b)
        # With bins
        out2_b = np.histogram_bin_edges(x, [125, 200] * u.cm)
        expected2_b = np.histogram_bin_edges(x.value, [1.25, 2.]) * x.unit
        assert np.all(out2_b == expected2_b)
        with pytest.raises(u.UnitsError):
            np.histogram_bin_edges(x, [125, 200] * u.s)

        with pytest.raises(u.UnitsError):
            np.histogram_bin_edges(x, [125, 200])

        with pytest.raises(u.UnitsError):
            np.histogram_bin_edges(x.value, [125, 200] * u.s)

    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_histogram2d(self):
        # TODO: add test cases for passing in bins
        x = np.array([1.1, 1.2, 1.3, 2.1, 5.1]) * u.m
        y = np.array([1.2, 2.2, 2.4, 3.0, 4.0]) * u.cm
        out_h, out_bx, out_by = np.histogram2d(x, y)
        expected_h, expected_bx, expected_by = np.histogram2d(x.value, y.value)
        expected_bx = expected_bx * x.unit
        expected_by = expected_by * y.unit
        assert np.all(out_h == expected_h)
        assert np.all(out_bx == expected_bx)
        assert np.all(out_by == expected_by)
        outd_h, outd_bx, outd_by = np.histogram2d(x, y, density=True)
        expectedd_h, expectedd_bx, expectedd_by = np.histogram2d(
            x.value, y.value, density=True)
        expectedd_h = expectedd_h / x.unit / y.unit
        expectedd_bx = expectedd_bx * x.unit
        expectedd_by = expectedd_by * y.unit
        assert np.all(outd_h == expectedd_h)
        assert np.all(outd_bx == expectedd_bx)
        assert np.all(outd_by == expectedd_by)

    @pytest.mark.xfail
    def test_histogramdd(self):
        # Postponing given various forms of passing in data, which
        # are not all dispatched correctly anyway:
        # https://github.com/numpy/numpy/issues/13728
        xyz = np.random.normal(size=(10, 3)) * u.m
        out_h, out_b = np.histogramdd(xyz)
        expected_h, expected_b = np.histogramdd(xyz.value)
        expected_b = expected_b * xyz.unit
        assert np.all(out_h == expected_h)
        assert np.all(out_b == expected_b)

    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_correlate(self):
        x1 = [1, 2, 3] * u.m
        x2 = [0, 1, 0.5] * u.m
        out = np.correlate(x1, x2)
        expected = np.correlate(x1.value, x2.value) * u.m ** 2
        assert np.all(out == expected)

    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_convolve(self):
        x1 = [1, 2, 3] * u.m
        x2 = [0, 1, 0.5] * u.m
        out = np.convolve(x1, x2)
        expected = np.convolve(x1.value, x2.value) * u.m ** 2
        assert np.all(out == expected)

    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_cov(self):
        # Do not see how we can use cov with Quantity
        x = np.array([[0, 2], [1, 1], [2, 0]]).T * u.m
        with pytest.raises(TypeError):
            np.cov(x)

    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_corrcoef(self):
        # Do not see how we can use cov with Quantity
        x = np.array([[0, 2], [1, 1], [2, 0]]).T * u.m
        with pytest.raises(TypeError):
            np.corrcoef(x)


class TestSortFunctions(InvariantUnitTestSetup):
    def test_sort(self):
        self.check(np.sort)

    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_sort_complex(self):
        self.check(np.sort_complex)

    def test_msort(self):
        self.check(np.msort)

    def test_partition(self):
        self.check(np.partition, 2)


class TestStringFunctions(metaclass=CoverageMeta):
    # For all these functions, we could change it to work on Quantity,
    # but it would mean deviating from the docstring.  Not clear whether
    # that is worth it.
    def setup(self):
        self.q = np.arange(3.) * u.Jy

    @pytest.mark.xfail
    def test_array2string(self):
        out = np.array2string(self.q)
        expected = str(self.q)
        assert out == expected

    @pytest.mark.xfail
    def test_array_repr(self):
        out = np.array_repr(self.q)
        expected = (np.array_repr(self.q.value)[:-1] +
                    ', {!r})'.format(str(self.q.unit)))
        assert out == expected

    @pytest.mark.xfail
    def test_array_str(self):
        out = np.array_str(self.q)
        expected = str(self.q)
        assert out == expected


class TestBitAndIndexFunctions(metaclass=CoverageMeta):
    # Index/bit functions generally fail for floats, so the usual
    # float quantity are safe, but the integer ones are not.
    def setup(self):
        self.q = np.arange(3) * u.m
        self.uint_q = u.Quantity(np.arange(3), 'm', dtype='u1')

    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_packbits(self):
        with pytest.raises(TypeError):
            np.packbits(self.q)
        with pytest.raises(TypeError):
            np.packbits(self.uint_q)

    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_unpackbits(self):
        with pytest.raises(TypeError):
            np.unpackbits(self.q)
        with pytest.raises(TypeError):
            np.unpackbits(self.uint_q)

    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_unravel_index(self):
        with pytest.raises(TypeError):
            np.unravel_index(self.q, 3)
        with pytest.raises(TypeError):
            np.unravel_index(self.uint_q, 3)

    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_ravel_multi_index(self):
        with pytest.raises(TypeError):
            np.ravel_multi_index((self.q,), 3)
        with pytest.raises(TypeError):
            np.ravel_multi_index((self.uint_q,), 3)

    @pytest.mark.xfail(NO_ARRAY_FUNCTION,
                       reason="Needs __array_function__ support")
    def test_ix_(self):
        with pytest.raises(TypeError):
            np.ix_(self.q)
        with pytest.raises(TypeError):
            np.ix_(self.uint_q)


class TestDtypeFunctions(NoUnitTestSetup):
    def test_common_type(self):
        self.check(np.common_type)

    def test_result_type(self):
        self.check(np.result_type)

    def test_can_cast(self):
        self.check(np.can_cast, self.q.dtype)
        self.check(np.can_cast, 'f4')

    def test_min_scalar_type(self):
        out = np.min_scalar_type(self.q[0])
        expected = np.min_scalar_type(self.q.value[0])
        assert out == expected

    def test_iscomplexobj(self):
        self.check(np.iscomplexobj)

    def test_isrealobj(self):
        self.check(np.isrealobj)


class TestMeshGrid(metaclass=CoverageMeta):
    def test_meshgrid(self):
        q1 = np.arange(3.) * u.m
        q2 = np.arange(5.) * u.s
        o1, o2 = np.meshgrid(q1, q2)
        e1, e2 = np.meshgrid(q1.value, q2.value)
        assert np.all(o1 == e1 * q1.unit)
        assert np.all(o2 == e2 * q2.unit)


class TestMemoryFunctions(NoUnitTestSetup):
    def test_shares_memory(self):
        self.check(np.shares_memory, self.q.value)

    def test_may_share_memory(self):
        self.check(np.may_share_memory, self.q.value)


class TestDatetimeFunctions(BasicTestSetup):
    def test_busday_count(self):
        with pytest.raises(TypeError):
            np.busday_count(self.q, self.q)

    def test_busday_offset(self):
        with pytest.raises(TypeError):
            np.busday_offset(self.q, self.q)

    def test_datetime_as_string(self):
        with pytest.raises(TypeError):
            np.datetime_as_string(self.q)

    def test_is_busday(self):
        with pytest.raises(TypeError):
            np.is_busday(self.q)


function_functions = {
    np.apply_along_axis, np.apply_over_axes,
    }
CoverageMeta.covered |= function_functions

financial_functions = {f for f in all_wrapped_functions.values()
                       if f in np.lib.financial.__dict__.values()}
CoverageMeta.covered |= financial_functions

deprecated_functions = {
    np.asscalar, np.rank
    }
CoverageMeta.covered |= deprecated_functions

io_functions = {np.save, np.savez, np.savetxt, np.savez_compressed}
CoverageMeta.covered |= io_functions

poly_functions = {
    np.poly, np.polyadd, np.polyder, np.polydiv, np.polyfit, np.polyint,
    np.polymul, np.polysub, np.polyval, np.roots, np.vander
    }
CoverageMeta.covered |= poly_functions

setops_functions = {f for f in all_wrapped_functions.values()
                    if f in np.lib.arraysetops.__dict__.values()}
CoverageMeta.covered |= setops_functions


@pytest.mark.xfail(NO_ARRAY_FUNCTION,
                   reason="no __array_function__ wrapping in numpy<1.17")
def test_completeness():
    assert set(all_wrapped_functions.values()) == CoverageMeta.covered
