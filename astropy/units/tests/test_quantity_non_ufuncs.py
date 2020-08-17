# Licensed under a 3-clause BSD style license - see LICENSE.rst
import itertools
import inspect

import numpy as np
from numpy.testing import assert_array_equal

import pytest

from astropy import units as u
from astropy.units.quantity_helper.function_helpers import (
    ARRAY_FUNCTION_ENABLED, SUBCLASS_SAFE_FUNCTIONS, UNSUPPORTED_FUNCTIONS,
    FUNCTION_HELPERS, DISPATCHED_FUNCTIONS, IGNORED_FUNCTIONS)
from astropy.utils.compat import NUMPY_LT_1_18, NUMPY_LT_1_20


needs_array_function = pytest.mark.xfail(
    not ARRAY_FUNCTION_ENABLED,
    reason="Needs __array_function__ support")


# To get the functions that could be covered, we look for those that
# are wrapped.  Of course, this does not give a full list pre-1.17.
def get_wrapped_functions(*modules):
    wrapped_functions = {}
    for mod in modules:
        for name, f in mod.__dict__.items():
            if f is np.printoptions:
                continue
            if callable(f) and hasattr(f, '__wrapped__'):
                wrapped_functions[name] = f
    return wrapped_functions


all_wrapped_functions = get_wrapped_functions(np, np.fft, np.linalg)
all_wrapped = set(all_wrapped_functions.values())


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
    # alen is deprecated in Numpy 1.8
    if NUMPY_LT_1_18:
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
        # Decided *not* to change default for subok for Quantity, since
        # that would be contrary to the docstring and might break code.
        self.check(np.broadcast_to, (3, 3, 3), subok=True)
        out = np.broadcast_to(self.q, (3, 3, 3))
        assert type(out) is np.ndarray  # NOT Quantity

    def test_broadcast_arrays(self):
        # Decided *not* to change default for subok for Quantity, since
        # that would be contrary to the docstring and might break code.
        q2 = np.ones((3, 3, 3)) / u.s
        o1, o2 = np.broadcast_arrays(self.q, q2, subok=True)
        assert isinstance(o1, u.Quantity)
        assert isinstance(o2, u.Quantity)
        assert o1.shape == o2.shape == (3, 3, 3)
        assert np.all(o1 == self.q)
        assert np.all(o2 == q2)
        a1, a2 = np.broadcast_arrays(self.q, q2)
        assert type(a1) is np.ndarray
        assert type(a2) is np.ndarray


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

    @needs_array_function
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

    @pytest.mark.parametrize('axis', (0, 1))
    def test_apply_along_axis(self, axis):
        out = np.apply_along_axis(np.square, axis, self.q)
        expected = np.apply_along_axis(np.square, axis,
                                       self.q.value) * self.q.unit ** 2
        assert_array_equal(out, expected)

    @needs_array_function
    @pytest.mark.parametrize('axes', ((1,), (0,), (0, 1)))
    def test_apply_over_axes(self, axes):
        def function(x, axis):
            return np.sum(np.square(x), axis)

        out = np.apply_over_axes(function, self.q, axes)
        expected = np.apply_over_axes(function, self.q.value, axes)
        expected = expected * self.q.unit ** (2 * len(axes))
        assert_array_equal(out, expected)


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
    @needs_array_function
    def test_copy(self):
        self.check(np.copy)
        # Also as kwarg
        copy = np.copy(a=self.q)
        assert_array_equal(copy, self.q)

    @needs_array_function
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

    @needs_array_function
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

    @needs_array_function
    def test_diag_1d_input(self):
        # Also check 1-D case; drops unit w/o __array_function__.
        q = self.q.ravel()
        o = np.diag(q)
        expected = np.diag(q.value) << q.unit
        assert o.unit == self.q.unit
        assert o.shape == expected.shape
        assert_array_equal(o, expected)

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

    @needs_array_function
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

    @needs_array_function
    def test_place(self):
        q = np.arange(3.) * u.m
        np.place(q, [True, False, True], [50, 150] * u.cm)
        assert q.unit == u.m
        expected = [50, 100, 150] * u.cm
        assert np.all(q == expected)

        a = np.arange(3.)
        np.place(a, [True, False, True], [50, 150] * u.percent)
        assert type(a) is np.ndarray
        expected = np.array([0.5, 1., 1.5])
        assert np.all(a == expected)

    @needs_array_function
    def test_copyto(self):
        q = np.arange(3.) * u.m
        np.copyto(q, [50, 0, 150] * u.cm, where=[True, False, True])
        assert q.unit == u.m
        expected = [50, 100, 150] * u.cm
        assert np.all(q == expected)

        a = np.arange(3.)
        np.copyto(a, [50, 0, 150] * u.percent, where=[True, False, True])
        assert type(a) is np.ndarray
        expected = np.array([0.5, 1., 1.5])
        assert np.all(a == expected)

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

    @needs_array_function
    def test_resize(self):
        self.check(np.resize, (4, 4))


class TestConcatenate(metaclass=CoverageMeta):
    def setup(self):
        self.q1 = np.arange(6.).reshape(2, 3) * u.m
        self.q2 = self.q1.to(u.cm)

    def check(self, func, *args, **kwargs):
        q_list = kwargs.pop('q_list', [self.q1, self.q2])
        q_ref = kwargs.pop('q_ref', q_list[0])
        o = func(q_list, *args, **kwargs)
        v_list = [q_ref._to_own_unit(q) for q in q_list]
        expected = func(v_list, *args, **kwargs) * q_ref.unit
        assert o.shape == expected.shape
        assert np.all(o == expected)

    @needs_array_function
    def test_concatenate(self):
        self.check(np.concatenate)
        self.check(np.concatenate, axis=1)

        self.check(np.concatenate, q_list=[np.zeros(self.q1.shape), self.q1, self.q2],
                   q_ref=self.q1)

        out = np.empty((4, 3)) * u.dimensionless_unscaled
        result = np.concatenate([self.q1, self.q2], out=out)
        assert out is result
        assert out.unit == self.q1.unit
        expected = np.concatenate(
            [self.q1.value, self.q2.to_value(self.q1.unit)]) * self.q1.unit
        assert np.all(result == expected)

        with pytest.raises(TypeError):
            np.concatenate([self.q1, object()])

    @needs_array_function
    def test_stack(self):
        self.check(np.stack)

    @needs_array_function
    def test_column_stack(self):
        self.check(np.column_stack)

    @needs_array_function
    def test_hstack(self):
        self.check(np.hstack)

    @needs_array_function
    def test_vstack(self):
        self.check(np.vstack)

    @needs_array_function
    def test_dstack(self):
        self.check(np.dstack)

    @needs_array_function
    def test_block(self):
        self.check(np.block)

        result = np.block([[0., 1.*u.m], [1.*u.cm, 2.*u.km]])
        assert np.all(result == np.block([[0, 1.], [.01, 2000.]]) << u.m)

    @needs_array_function
    def test_append(self):
        out = np.append(self.q1, self.q2, axis=0)
        assert out.unit == self.q1.unit
        expected = np.append(self.q1.value, self.q2.to_value(self.q1.unit),
                             axis=0) * self.q1.unit
        assert np.all(out == expected)

        a = np.arange(3.)
        result = np.append(a, 50. * u.percent)
        assert isinstance(result, u.Quantity)
        assert result.unit == u.dimensionless_unscaled
        expected = np.append(a, 0.5) * u.dimensionless_unscaled
        assert np.all(result == expected)

    @needs_array_function
    def test_insert(self):
        # Unit of inserted values is not ignored.
        q = np.arange(12.).reshape(6, 2) * u.m
        out = np.insert(q, (3, 5), [50., 25.] * u.cm)
        assert isinstance(out, u.Quantity)
        assert out.unit == q.unit
        expected = np.insert(q.value, (3, 5), [0.5, 0.25]) << q.unit
        assert np.all(out == expected)
        # 0 can have any unit.
        out2 = np.insert(q, (3, 5), 0)
        expected2 = np.insert(q.value, (3, 5), 0) << q.unit
        assert np.all(out2 == expected2)

        a = np.arange(3.)
        result = np.insert(a, (2,), 50. * u.percent)
        assert isinstance(result, u.Quantity)
        assert result.unit == u.dimensionless_unscaled
        expected = np.insert(a, (2,), 0.5) * u.dimensionless_unscaled
        assert np.all(result == expected)

        with pytest.raises(TypeError):
            np.insert(q, 3 * u.cm, 50. * u.cm)
        with pytest.raises(u.UnitsError):
            np.insert(q, (3, 5), 0. * u.s)

    @needs_array_function
    def test_pad(self):
        q = np.arange(1., 6.) * u.m
        out = np.pad(q, (2, 3), 'constant', constant_values=(0., 150.*u.cm))
        assert out.unit == q.unit
        expected = np.pad(q.value, (2, 3), 'constant',
                          constant_values=(0., 1.5)) * q.unit
        assert np.all(out == expected)
        out2 = np.pad(q, (2, 3), 'constant', constant_values=150.*u.cm)
        assert out2.unit == q.unit
        expected2 = np.pad(q.value, (2, 3), 'constant',
                           constant_values=1.5) * q.unit
        assert np.all(out2 == expected2)
        out3 = np.pad(q, (2, 3), 'linear_ramp', end_values=(25.*u.cm, 0.))
        assert out3.unit == q.unit
        expected3 = np.pad(q.value, (2, 3), 'linear_ramp',
                           end_values=(0.25, 0.)) * q.unit
        assert np.all(out3 == expected3)


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
        with pytest.raises(TypeError):
            np.any(self.q)

    def test_all(self):
        with pytest.raises(TypeError):
            np.all(self.q)

    def test_sometrue(self):
        with pytest.raises(TypeError):
            np.sometrue(self.q)

    def test_alltrue(self):
        with pytest.raises(TypeError):
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

    @needs_array_function
    def test_sinc(self):
        q = [0., 3690., -270., 690.] * u.deg
        out = np.sinc(q)
        expected = np.sinc(q.to_value(u.radian)) * u.one
        assert isinstance(out, u.Quantity)
        assert np.all(out == expected)
        with pytest.raises(u.UnitsError):
            np.sinc(1.*u.one)

    @needs_array_function
    def test_where(self):
        out = np.where([True, False, True], self.q, 1. * u.km)
        expected = np.where([True, False, True], self.q.value,
                            1000.) * self.q.unit
        assert np.all(out == expected)

    @needs_array_function
    def test_choose(self):
        # from np.choose docstring
        a = np.array([0, 1]).reshape((2, 1, 1))
        q1 = np.array([1, 2, 3]).reshape((1, 3, 1)) * u.cm
        q2 = np.array([-1, -2, -3, -4, -5]).reshape((1, 1, 5)) * u.m
        out = np.choose(a, (q1, q2))
        # result is 2x3x5, res[0,:,:]=c1, res[1,:,:]=c2
        expected = np.choose(a, (q1.value, q2.to_value(q1.unit))) * u.cm
        assert np.all(out == expected)

    @needs_array_function
    def test_select(self):
        q = self.q
        out = np.select([q < 0.55 * u.m, q > 1. * u.m],
                        [q, q.to(u.cm)], default=-1. * u.km)
        expected = np.select([q.value < 0.55, q.value > 1],
                             [q.value, q.value], default=-1000) * u.m
        assert np.all(out == expected)

    @needs_array_function
    def test_real_if_close(self):
        q = np.array([1+0j, 0+1j, 1+1j, 0+0j]) * u.m
        out = np.real_if_close(q)
        expected = np.real_if_close(q.value) * u.m
        assert np.all(out == expected)

    @needs_array_function
    def test_tril(self):
        self.check(np.tril)

    @needs_array_function
    def test_triu(self):
        self.check(np.triu)

    @needs_array_function
    def test_unwrap(self):
        q = [0., 3690., -270., 690.] * u.deg
        out = np.unwrap(q)
        expected = (np.unwrap(q.to_value(u.rad)) * u.rad).to(q.unit)
        assert out.unit == expected.unit
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

    @needs_array_function
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

    @needs_array_function
    def test_allclose_atol_default_unit(self):
        q_cm = self.q.to(u.cm)
        out = np.isclose(self.q, q_cm)
        expected = np.isclose(self.q.value, q_cm.to_value(u.m))
        assert np.all(out == expected)
        q1 = np.arange(3.) * u.m
        q2 = np.array([0., 101., 198.]) * u.cm
        out = np.isclose(q1, q2, atol=0.011, rtol=0)
        expected = np.isclose(q1.value, q2.to_value(q1.unit),
                              atol=0.011, rtol=0)
        assert np.all(out == expected)
        out2 = np.isclose(q2, q1, atol=0.011, rtol=0)
        expected2 = np.isclose(q2.value, q1.to_value(q2.unit),
                               atol=0.011, rtol=0)
        assert np.all(out2 == expected2)


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

    @needs_array_function
    def test_quantile(self):
        self.check(np.quantile, 0.5)
        o = np.quantile(self.q, 50 * u.percent)
        expected = np.quantile(self.q.value, 0.5) * u.m
        assert np.all(o == expected)
        # For ndarray input, we return a Quantity.
        o2 = np.quantile(self.q.value, 50 * u.percent)
        assert o2.unit == u.dimensionless_unscaled
        assert np.all(o2 == expected.value)
        o3 = 0 * o2
        result = np.quantile(self.q, 50 * u.percent, out=o3)
        assert result is o3
        assert np.all(o3 == expected)
        o4 = 0 * o2
        result = np.quantile(self.q, 50 * u.percent, None, o4)
        assert result is o4
        assert np.all(o4 == expected)

    @needs_array_function
    def test_percentile(self):
        self.check(np.percentile, 0.5)
        o = np.percentile(self.q, 0.5 * u.one)
        expected = np.percentile(self.q.value, 50) * u.m
        assert np.all(o == expected)

    def test_trace(self):
        self.check(np.trace)

    @needs_array_function
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

    @needs_array_function
    def test_allclose_atol_default_unit(self):
        q1 = np.arange(3.) * u.m
        q2 = np.array([0., 101., 199.]) * u.cm
        assert np.allclose(q1, q2, atol=0.011, rtol=0)
        assert not np.allclose(q2, q1, atol=0.011, rtol=0)

    def test_allclose_failures(self):
        q1 = np.arange(3.) * u.m
        q2 = np.array([0., 101., 199.]) * u.cm
        with pytest.raises(u.UnitsError):
            np.allclose(q1, q2, atol=2*u.s, rtol=0)
        with pytest.raises(u.UnitsError):
            np.allclose(q1, q2, atol=0, rtol=1.*u.s)

    @needs_array_function
    def test_array_equal(self):
        q1 = np.arange(3.) * u.m
        q2 = q1.to(u.cm)
        assert np.array_equal(q1, q2)
        q3 = q1.value * u.cm
        assert not np.array_equal(q1, q3)

    @needs_array_function
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

    @needs_array_function
    def test_nanquantile(self):
        self.check(np.nanquantile, 0.5)
        o = np.nanquantile(self.q, 50 * u.percent)
        expected = np.nanquantile(self.q.value, 0.5) * u.m
        assert np.all(o == expected)

    @needs_array_function
    def test_nanpercentile(self):
        self.check(np.nanpercentile, 0.5)
        o = np.nanpercentile(self.q, 0.5 * u.one)
        expected = np.nanpercentile(self.q.value, 50) * u.m
        assert np.all(o == expected)


class TestVariousProductFunctions(metaclass=CoverageMeta):
    """
    Test functions that are similar to gufuncs
    """
    @needs_array_function
    def test_cross(self):
        q1 = np.arange(6.).reshape(2, 3) * u.m
        q2 = np.array([4., 5., 6.]) / u.s
        o = np.cross(q1, q2)
        expected = np.cross(q1.value, q2.value) * u.m / u.s
        assert np.all(o == expected)

    @needs_array_function
    def test_outer(self):
        q1 = np.array([1, 2, 3]) * u.m
        q2 = np.array([1, 2]) / u.s
        o = np.outer(q1, q2)
        assert np.all(o == np.array([[1, 2], [2, 4], [3, 6]]) * u.m / u.s)

        o2 = 0 * o
        result = np.outer(q1, q2, out=o2)
        assert result is o2
        assert np.all(o2 == o)

        with pytest.raises(TypeError):
            np.outer(q1, q2, out=object())

    @needs_array_function
    def test_inner(self):
        q1 = np.array([1, 2, 3]) * u.m
        q2 = np.array([4, 5, 6]) / u.s
        o = np.inner(q1, q2)
        assert o == 32 * u.m / u.s

    @needs_array_function
    def test_dot(self):
        q1 = np.array([1., 2., 3.]) * u.m
        q2 = np.array([4., 5., 6.]) / u.s
        o = np.dot(q1, q2)
        assert o == 32. * u.m / u.s

    @needs_array_function
    def test_vdot(self):
        q1 = np.array([1j, 2j, 3j]) * u.m
        q2 = np.array([4j, 5j, 6j]) / u.s
        o = np.vdot(q1, q2)
        assert o == (32. + 0j) * u.m / u.s

    @needs_array_function
    def test_tensordot(self):
        # From the docstring example
        a = np.arange(60.).reshape(3, 4, 5) * u.m
        b = np.arange(24.).reshape(4, 3, 2) / u.s
        c = np.tensordot(a, b, axes=([1, 0], [0, 1]))
        expected = np.tensordot(a.value, b.value,
                                axes=([1, 0], [0, 1])) * u.m / u.s
        assert np.all(c == expected)

    @needs_array_function
    def test_kron(self):
        q1 = np.eye(2) * u.m
        q2 = np.ones(2) / u.s
        o = np.kron(q1, q2)
        expected = np.kron(q1.value, q2.value) * u.m / u.s
        assert np.all(o == expected)

    @needs_array_function
    def test_einsum(self):
        q1 = np.arange(9.).reshape(3, 3) * u.m
        o = np.einsum('...i', q1)
        assert np.all(o == q1)
        o = np.einsum('ii', q1)
        expected = np.einsum('ii', q1.value) * u.m
        assert np.all(o == expected)
        q2 = np.eye(3) / u.s
        o2 = np.einsum('ij,jk', q1, q2)
        assert np.all(o2 == q1 / u.s)
        o3 = 0 * o2
        result = np.einsum('ij,jk', q1, q2, out=o3)
        assert result is o3
        assert np.all(o3 == o2)

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

    @needs_array_function
    def test_diff_prepend_append(self):
        x = np.arange(10.) * u.m
        out = np.diff(x, prepend=-12.5*u.cm, append=1*u.km)
        expected = np.diff(x.value, prepend=-0.125, append=1000.) * x.unit
        assert np.all(out == expected)
        x = np.arange(10.) * u.m
        out = np.diff(x, prepend=-12.5*u.cm, append=1*u.km, n=2)
        expected = np.diff(x.value, prepend=-0.125, append=1000.,
                           n=2) * x.unit
        assert np.all(out == expected)

        with pytest.raises(TypeError):
            np.diff(x, prepend=object())

    def test_gradient(self):
        # Simple gradient works out of the box.
        x = np.arange(10.) * u.m
        out = np.gradient(x)
        expected = np.gradient(x.value) * u.m
        assert np.all(out == expected)

    @needs_array_function
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
        dfdy2 = np.gradient(f, y, axis=(1,))
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

    @needs_array_function
    def test_logspace(self):
        unit = u.m / u.s**2
        out = np.logspace(10.*u.dex(unit), 20*u.dex(unit), 10)
        expected = np.logspace(10., 20., 10) * unit
        assert np.all(out == expected)
        out = np.logspace(10.*u.STmag, 20*u.STmag, 10)
        expected = np.logspace(10., 20., 10, base=10.**(-0.4)) * u.ST
        assert u.allclose(out, expected)

    @needs_array_function
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
    @needs_array_function
    def test_interp(self):
        x = np.array([1250., 2750.]) * u.m
        xp = np.arange(5.) * u.km
        yp = np.arange(5.) * u.day
        out = np.interp(x, xp, yp)
        expected = np.interp(x.to_value(xp.unit), xp.value, yp.value) * yp.unit
        assert np.all(out == expected)

        out = np.interp(x, xp, yp.value)
        assert type(out) is np.ndarray
        assert np.all(out == expected.value)

    @needs_array_function
    def test_piecewise(self):
        x = np.linspace(-2.5, 2.5, 6) * u.m
        out = np.piecewise(x, [x < 0, x >= 0], [-1*u.s, 1*u.day])
        expected = np.piecewise(x.value, [x.value < 0, x.value >= 0],
                                [-1, 24*3600]) * u.s
        assert out.unit == expected.unit
        assert np.all(out == expected)

        out2 = np.piecewise(x, [x < 1 * u.m, x >= 0],
                            [-1*u.s, 1*u.day, lambda x: 1*u.hour])
        expected2 = np.piecewise(x.value, [x.value < 1, x.value >= 0],
                                 [-1, 24*3600, 3600]) * u.s
        assert out2.unit == expected2.unit
        assert np.all(out2 == expected2)

        out3 = np.piecewise(x, [x < 1 * u.m, x >= 0],
                            [0, 1*u.percent, lambda x: 1*u.one])
        expected3 = np.piecewise(x.value, [x.value < 1, x.value >= 0],
                                 [0, 0.01, 1]) * u.one
        assert out3.unit == expected3.unit
        assert np.all(out3 == expected3)

        with pytest.raises(TypeError):  # no Quantity in condlist.
            np.piecewise(x, [x], [0.])

        with pytest.raises(TypeError):  # no Quantity in condlist.
            np.piecewise(x.value, [x], [0.])


class TestBincountDigitize(metaclass=CoverageMeta):
    @needs_array_function
    def test_bincount(self):
        i = np.array([1, 1, 2, 3, 2, 4])
        weights = np.arange(len(i)) * u.Jy
        out = np.bincount(i, weights)
        expected = np.bincount(i, weights.value) * weights.unit
        assert_array_equal(out, expected)

        with pytest.raises(TypeError):
            np.bincount(weights)

    @needs_array_function
    def test_digitize(self):
        x = np.array([1500., 2500., 4500.]) * u.m
        bins = np.arange(10.) * u.km
        out = np.digitize(x, bins)
        expected = np.digitize(x.to_value(bins.unit), bins.value)
        assert_array_equal(out, expected)


class TestHistogramFunctions(metaclass=CoverageMeta):

    def setup(self):
        self.x = np.array([1.1, 1.2, 1.3, 2.1, 5.1]) * u.m
        self.y = np.array([1.2, 2.2, 2.4, 3.0, 4.0]) * u.cm
        self.weights = np.arange(len(self.x)) / u.s

    def check(self, function, *args, value_args=None, value_kwargs=None,
              expected_units=None, **kwargs):
        """Check quanties are treated correctly in the histogram function.
        Test is done by applying ``function(*args, **kwargs)``, where
        the argument can be quantities, and comparing the result to
        ``function(*value_args, **value_kwargs)``, with the outputs
        converted to quantities using the ``expected_units`` (where `None`
        indicates the output is expected to be a regular array).

        For ``**value_kwargs``, any regular ``kwargs`` are treated as
        defaults, i.e., non-quantity arguments do not have to be repeated.
        """
        if value_kwargs is None:
            value_kwargs = kwargs
        else:
            for k, v in kwargs.items():
                value_kwargs.setdefault(k, v)
        # Get the result, using the Quantity override.
        out = function(*args, **kwargs)
        # Get the comparison, with non-Quantity arguments.
        expected = function(*value_args, **value_kwargs)
        # All histogram functions return a tuple of the actual histogram
        # and the bin edges.  First, check the actual histogram.
        out_h = out[0]
        expected_h = expected[0]
        if expected_units[0] is not None:
            expected_h = expected_h * expected_units[0]
        assert_array_equal(out_h, expected_h)
        # Check bin edges.  Here, histogramdd returns an interable of the
        # bin edges as the second return argument, while histogram and
        # histogram2d return the bin edges directly.
        if function is np.histogramdd:
            bin_slice = 1
        else:
            bin_slice = slice(1, None)

        for o_bin, e_bin, e_unit in zip(out[bin_slice],
                                        expected[bin_slice],
                                        expected_units[bin_slice]):
            if e_unit is not None:
                e_bin = e_bin * e_unit
            assert_array_equal(o_bin, e_bin)

    @needs_array_function
    def test_histogram(self):
        x = self.x
        weights = self.weights
        # Plain histogram.
        self.check(np.histogram, x,
                   value_args=(x.value,),
                   expected_units=(None, x.unit))
        # With bins.
        self.check(np.histogram, x, [125, 200] * u.cm,
                   value_args=(x.value, [1.25, 2.]),
                   expected_units=(None, x.unit))
        # With density.
        self.check(np.histogram, x, [125, 200] * u.cm, density=True,
                   value_args=(x.value, [1.25, 2.]),
                   expected_units=(1/x.unit, x.unit))
        # With weights.
        self.check(np.histogram, x, [125, 200] * u.cm, weights=weights,
                   value_args=(x.value, [1.25, 2.]),
                   value_kwargs=dict(weights=weights.value),
                   expected_units=(weights.unit, x.unit))
        # With weights and density.
        self.check(np.histogram, x, [125, 200] * u.cm,
                   weights=weights, density=True,
                   value_args=(x.value, [1.25, 2.]),
                   value_kwargs=dict(weights=weights.value),
                   expected_units=(weights.unit/x.unit, x.unit))

        with pytest.raises(u.UnitsError):
            np.histogram(x, [125, 200] * u.s)

        with pytest.raises(u.UnitsError):
            np.histogram(x, [125, 200])

        with pytest.raises(u.UnitsError):
            np.histogram(x.value, [125, 200] * u.s)

    @needs_array_function
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

    @needs_array_function
    def test_histogram2d(self):
        x, y = self.x, self.y
        weights = self.weights
        # Basic tests with X, Y.
        self.check(np.histogram2d, x, y,
                   value_args=(x.value, y.value),
                   expected_units=(None, x.unit, y.unit))
        # Check units with density.
        self.check(np.histogram2d, x, y, density=True,
                   value_args=(x.value, y.value),
                   expected_units=(1/(x.unit*y.unit), x.unit, y.unit))
        # Check units with weights.
        self.check(np.histogram2d, x, y, weights=weights,
                   value_args=(x.value, y.value),
                   value_kwargs=dict(weights=weights.value),
                   expected_units=(weights.unit, x.unit, y.unit))
        # Check quantity bin sizes.
        inb_y = [0, 0.025, 1.] * u.m
        self.check(np.histogram2d, x, y, [5, inb_y],
                   value_args=(x.value, y.value,
                               [5, np.array([0, 2.5, 100.])]),
                   expected_units=(None, x.unit, y.unit))
        # Check we dispatch on bin sizes (and check kwarg as well).
        inb2_y = [0, 250, 10000.] * u.percent
        self.check(np.histogram2d, x.value, y.value, bins=[5, inb2_y],
                   value_args=(x.value, y.value),
                   value_kwargs=dict(bins=[5, np.array([0, 2.5, 100.])]),
                   expected_units=(None, u.one, u.one))

        # Single-item bins should be integer, not Quantity.
        with pytest.raises(TypeError):
            np.histogram2d(x, y, 125 * u.s)

        with pytest.raises(TypeError):
            np.histogram2d(x.value, y.value, 125 * u.s)

        # Bin units need to match units of x, y.
        with pytest.raises(u.UnitsError):
            np.histogram2d(x, y, [125, 200] * u.s)

        with pytest.raises(u.UnitsError):
            np.histogram2d(x, y, ([125, 200], [125, 200]))

        with pytest.raises(u.UnitsError):
            np.histogram2d(x.value, y.value, [125, 200] * u.s)

    @needs_array_function
    def test_histogramdd(self):
        # First replicates of the histogram2d tests, but using the
        # histogramdd override.  Normally takes the sample as a tuple
        # with a given number of dimensions, and returns the histogram
        # as well as a tuple of bin edges.
        sample = self.x, self.y
        sample_units = self.x.unit, self.y.unit
        sample_values = (self.x.value, self.y.value)
        weights = self.weights
        # Basic tests with X, Y
        self.check(np.histogramdd, sample,
                   value_args=(sample_values,),
                   expected_units=(None, sample_units))
        # Check units with density.
        self.check(np.histogramdd, sample, density=True,
                   value_args=(sample_values,),
                   expected_units=(1/(self.x.unit*self.y.unit),
                                   sample_units))
        # Check units with weights.
        self.check(np.histogramdd, sample, weights=weights,
                   value_args=(sample_values,),
                   value_kwargs=dict(weights=weights.value),
                   expected_units=(weights.unit, sample_units))
        # Check quantity bin sizes.
        inb_y = [0, 0.025, 1.] * u.m
        self.check(np.histogramdd, sample, [5, inb_y],
                   value_args=(sample_values, [5, np.array([0, 2.5, 100.])]),
                   expected_units=(None, sample_units))
        # Check we dispatch on bin sizes (and check kwarg as well).
        inb2_y = [0, 250, 10000.] * u.percent
        self.check(np.histogramdd, sample_values, bins=[5, inb2_y],
                   value_args=(sample_values,),
                   value_kwargs=dict(bins=[5, np.array([0, 2.5, 100.])]),
                   expected_units=(None, (u.one, u.one)))
        # For quantities, it is probably not that likely one would pass
        # in the sample as an array, but check that it works anyway.
        # This also gives a 3-D check.
        xyz = np.random.normal(size=(10, 3)) * u.m
        self.check(np.histogramdd, xyz,
                   value_args=(xyz.value,),
                   expected_units=(None, (xyz.unit,)*3))
        # Passing it in as a tuple should work just as well; note the
        # *last* axis contains the sample dimension.
        self.check(np.histogramdd, (xyz[:, 0], xyz[:, 1], xyz[:, 2]),
                   value_args=(xyz.value,),
                   expected_units=(None, (xyz.unit,)*3))

        # Single-item bins should be integer, not Quantity.
        with pytest.raises(TypeError):
            np.histogramdd(sample, 125 * u.s)

        # Sequence of single items should be integer.
        with pytest.raises(TypeError):
            np.histogramdd(sample, [125, 200] * u.s)

        with pytest.raises(TypeError):
            np.histogramdd(sample_values, [125, 200] * u.s)

        # Units of bins should match.
        with pytest.raises(u.UnitsError):
            np.histogramdd(sample, ([125, 200], [125, 200]))

        with pytest.raises(u.UnitsError):
            np.histogramdd(sample_values, ([125, 200] * u.s, [125, 200]))

    @needs_array_function
    def test_correlate(self):
        x1 = [1, 2, 3] * u.m
        x2 = [0, 1, 0.5] * u.m
        out = np.correlate(x1, x2)
        expected = np.correlate(x1.value, x2.value) * u.m ** 2
        assert np.all(out == expected)

    @needs_array_function
    def test_convolve(self):
        x1 = [1, 2, 3] * u.m
        x2 = [0, 1, 0.5] * u.m
        out = np.convolve(x1, x2)
        expected = np.convolve(x1.value, x2.value) * u.m ** 2
        assert np.all(out == expected)

    @needs_array_function
    def test_cov(self):
        # Do not see how we can use cov with Quantity
        x = np.array([[0, 2], [1, 1], [2, 0]]).T * u.m
        with pytest.raises(TypeError):
            np.cov(x)

    @needs_array_function
    def test_corrcoef(self):
        # Do not see how we can use cov with Quantity
        x = np.array([[0, 2], [1, 1], [2, 0]]).T * u.m
        with pytest.raises(TypeError):
            np.corrcoef(x)


class TestSortFunctions(InvariantUnitTestSetup):
    def test_sort(self):
        self.check(np.sort)

    @needs_array_function
    def test_sort_complex(self):
        self.check(np.sort_complex)

    def test_msort(self):
        self.check(np.msort)

    def test_partition(self):
        self.check(np.partition, 2)


class TestStringFunctions(metaclass=CoverageMeta):
    # For these, making behaviour work means deviating only slightly from
    # the docstring, and by default they fail miserably.  So, might as well.
    def setup(self):
        self.q = np.arange(3.) * u.Jy

    @needs_array_function
    def test_array2string(self):
        # The default formatters cannot handle units, so if we do not pass
        # a relevant formatter, we are better off just treating it as an
        # array (which happens for all subtypes).
        out0 = np.array2string(self.q)
        expected0 = str(self.q.value)
        assert out0 == expected0
        # Arguments are interpreted as usual.
        out1 = np.array2string(self.q, separator=', ')
        expected1 = '[0., 1., 2.]'
        assert out1 == expected1
        # If we do pass in a formatter, though, it should be used.
        out2 = np.array2string(self.q, separator=', ', formatter={'all': str})
        expected2 = '[0.0 Jy, 1.0 Jy, 2.0 Jy]'
        assert out2 == expected2
        # Also as positional argument (no, nobody will do this!)
        out3 = np.array2string(self.q, None, None, None, ', ', '',
                               np._NoValue, {'float': str})
        assert out3 == expected2
        # But not if the formatter is not relevant for us.
        out4 = np.array2string(self.q, separator=', ', formatter={'int': str})
        assert out4 == expected1

    @needs_array_function
    def test_array_repr(self):
        out = np.array_repr(self.q)
        assert out == "Quantity([0., 1., 2.], unit='Jy')"
        q2 = self.q.astype('f4')
        out2 = np.array_repr(q2)
        assert out2 == "Quantity([0., 1., 2.], unit='Jy', dtype=float32)"

    @needs_array_function
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

    @needs_array_function
    def test_packbits(self):
        with pytest.raises(TypeError):
            np.packbits(self.q)
        with pytest.raises(TypeError):
            np.packbits(self.uint_q)

    @needs_array_function
    def test_unpackbits(self):
        with pytest.raises(TypeError):
            np.unpackbits(self.q)
        with pytest.raises(TypeError):
            np.unpackbits(self.uint_q)

    @needs_array_function
    def test_unravel_index(self):
        with pytest.raises(TypeError):
            np.unravel_index(self.q, 3)
        with pytest.raises(TypeError):
            np.unravel_index(self.uint_q, 3)

    @needs_array_function
    def test_ravel_multi_index(self):
        with pytest.raises(TypeError):
            np.ravel_multi_index((self.q,), 3)
        with pytest.raises(TypeError):
            np.ravel_multi_index((self.uint_q,), 3)

    @needs_array_function
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


class TestSetOpsFcuntions(metaclass=CoverageMeta):
    def setup(self):
        self.q = np.array([[0., 1., -1.],
                           [3., 5., 3.],
                           [0., 1., -1]]) * u.m
        self.q2 = np.array([0., 100., 150., 200.]) * u.cm

    def check(self, function, qs, *args, **kwargs):
        unit = kwargs.pop('unit', self.q.unit)
        out = function(*qs, *args, **kwargs)
        qv = tuple(q.to_value(self.q.unit) for q in qs)
        expected = function(*qv, *args, **kwargs)
        if isinstance(expected, tuple):
            if unit:
                expected = (expected[0] * unit,) + expected[1:]
            for o, e in zip(out, expected):
                assert_array_equal(o, e)
        else:
            if unit:
                expected = expected * unit
            assert_array_equal(out, expected)

    def check1(self, function, *args, **kwargs):
        self.check(function, (self.q,), *args, **kwargs)

    def check2(self, function, *args, **kwargs):
        self.check(function, (self.q, self.q2), *args, **kwargs)

    @pytest.mark.parametrize('kwargs', (
        dict(return_index=True, return_inverse=True),
        dict(return_counts=True),
        dict(return_index=True, return_inverse=True, return_counts=True)))
    def test_unique(self, kwargs):
        self.check1(np.unique, **kwargs)

    @needs_array_function
    @pytest.mark.parametrize('kwargs', (
        dict(axis=0),
        dict(axis=1),
        dict(return_counts=True, return_inverse=False, axis=1)))
    def test_unique_more_complex(self, kwargs):
        self.check1(np.unique, **kwargs)

    @needs_array_function
    @pytest.mark.parametrize('kwargs', (
        dict(),
        dict(return_indices=True)))
    def test_intersect1d(self, kwargs):
        self.check2(np.intersect1d, **kwargs)

    @needs_array_function
    def test_setxor1d(self):
        self.check2(np.setxor1d)

    @needs_array_function
    def test_union1d(self):
        self.check2(np.union1d)
        result = np.union1d(np.array([0., np.nan]), np.arange(3) << u.m)
        assert result.unit is u.m
        assert_array_equal(result.value, np.array([0., 1., 2., np.nan]))

    @needs_array_function
    def test_setdiff1d(self):
        self.check2(np.setdiff1d)

    @needs_array_function
    def test_in1d(self):
        self.check2(np.in1d, unit=None)
        # Check zero is treated as having any unit.
        assert np.in1d(np.zeros(1), self.q2)
        with pytest.raises(u.UnitsError):
            np.in1d(np.ones(1), self.q2)

    @needs_array_function
    def test_isin(self):
        self.check2(np.isin, unit=None)

    def test_ediff1d(self):
        # ediff1d works always as it calls the Quantity method.
        self.check1(np.ediff1d)
        x = np.arange(10.) * u.m
        out = np.ediff1d(x, to_begin=-12.5*u.cm, to_end=1*u.km)
        expected = np.ediff1d(x.value, to_begin=-0.125, to_end=1000.) * x.unit
        assert_array_equal(out, expected)


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


# These functions always worked; ensure they do not regress.
# Note that they are *not* wrapped so no need to check coverage.
@pytest.mark.parametrize('function', [np.fft.fftfreq, np.fft.rfftfreq])
def test_fft_frequencies(function):
    out = function(128, d=0.1*u.s)
    expected = function(128, d=0.1) / u.s
    assert_array_equal(out, expected)


@needs_array_function
class TestFFT(InvariantUnitTestSetup):
    # These are all trivial, just preserve the unit.
    def setup(self):
        # Use real input; gets turned into complex as needed.
        self.q = np.arange(128.).reshape(8, -1) * u.s

    def test_fft(self):
        self.check(np.fft.fft)

    def test_ifft(self):
        self.check(np.fft.ifft)

    def test_rfft(self):
        self.check(np.fft.rfft)

    def test_irfft(self):
        self.check(np.fft.irfft)

    def test_fft2(self):
        self.check(np.fft.fft2)

    def test_ifft2(self):
        self.check(np.fft.ifft2)

    def test_rfft2(self):
        self.check(np.fft.rfft2)

    def test_irfft2(self):
        self.check(np.fft.irfft2)

    def test_fftn(self):
        self.check(np.fft.fftn)

    def test_ifftn(self):
        self.check(np.fft.ifftn)

    def test_rfftn(self):
        self.check(np.fft.rfftn)

    def test_irfftn(self):
        self.check(np.fft.irfftn)

    def test_hfft(self):
        self.check(np.fft.hfft)

    def test_ihfft(self):
        self.check(np.fft.ihfft)

    def test_fftshift(self):
        self.check(np.fft.fftshift)

    def test_ifftshift(self):
        self.check(np.fft.ifftshift)


class TestLinAlg(metaclass=CoverageMeta):
    def setup(self):
        # Use a matrix safe for inversion, etc.
        self.q = np.array([[1., -1., 2.],
                           [0., 3., -1.],
                           [-1., -1., 1.]]) << u.m

    def test_cond(self):
        c = np.linalg.cond(self.q)
        expected = np.linalg.cond(self.q.value)
        assert c == expected

    def test_matrix_rank(self):
        r = np.linalg.matrix_rank(self.q)
        x = np.linalg.matrix_rank(self.q.value)
        assert r == x

    @needs_array_function
    def test_matrix_rank_with_tol(self):
        # Use a matrix that is not so good, so tol=1 and tol=0.01 differ.
        q = np.arange(9.).reshape(3, 3) / 4 * u.m
        tol = 1. * u.cm
        r2 = np.linalg.matrix_rank(q, tol)
        x2 = np.linalg.matrix_rank(q.value, tol.to_value(q.unit))
        assert r2 == x2

    def test_matrix_power(self):
        q1 = np.linalg.matrix_power(self.q, 1)
        assert_array_equal(q1, self.q)
        q2 = np.linalg.matrix_power(self.q, 2)
        assert_array_equal(q2, self.q @ self.q)
        q2 = np.linalg.matrix_power(self.q, 4)
        assert_array_equal(q2, self.q @ self.q @ self.q @ self.q)

    @needs_array_function
    def test_matrix_inv_power(self):
        qinv = np.linalg.inv(self.q.value) / self.q.unit
        qm1 = np.linalg.matrix_power(self.q, -1)
        assert_array_equal(qm1, qinv)
        qm3 = np.linalg.matrix_power(self.q, -3)
        assert_array_equal(qm3, qinv @ qinv @ qinv)

    @needs_array_function
    def test_multi_dot(self):
        q2 = np.linalg.multi_dot([self.q, self.q])
        q2x = self.q @ self.q
        assert_array_equal(q2, q2x)
        q3 = np.linalg.multi_dot([self.q, self.q, self.q])
        q3x = self.q @ self.q @ self.q
        assert_array_equal(q3, q3x)

    @needs_array_function
    def test_svd(self):
        m = np.arange(10.) * np.arange(5.)[:, np.newaxis] * u.m
        svd_u, svd_s, svd_vt = np.linalg.svd(m, full_matrices=False)
        svd_ux, svd_sx, svd_vtx = np.linalg.svd(m.value, full_matrices=False)
        svd_sx <<= m.unit
        assert_array_equal(svd_u, svd_ux)
        assert_array_equal(svd_vt, svd_vtx)
        assert_array_equal(svd_s, svd_sx)
        assert u.allclose(svd_u @ np.diag(svd_s) @ svd_vt, m)

        s2 = np.linalg.svd(m, compute_uv=False)
        svd_s2x = np.linalg.svd(m.value, compute_uv=False) << m.unit
        assert_array_equal(s2, svd_s2x)

    @needs_array_function
    def test_inv(self):
        inv = np.linalg.inv(self.q)
        expected = np.linalg.inv(self.q.value) / self.q.unit
        assert_array_equal(inv, expected)

    @needs_array_function
    def test_pinv(self):
        pinv = np.linalg.pinv(self.q)
        expected = np.linalg.pinv(self.q.value) / self.q.unit
        assert_array_equal(pinv, expected)
        rcond = 0.01 * u.cm
        pinv2 = np.linalg.pinv(self.q, rcond)
        expected2 = np.linalg.pinv(self.q.value,
                                   rcond.to_value(self.q.unit)) / self.q.unit
        assert_array_equal(pinv2, expected2)

    @needs_array_function
    def test_tensorinv(self):
        inv = np.linalg.tensorinv(self.q, ind=1)
        expected = np.linalg.tensorinv(self.q.value, ind=1) / self.q.unit
        assert_array_equal(inv, expected)

    @needs_array_function
    def test_det(self):
        det = np.linalg.det(self.q)
        expected = np.linalg.det(self.q.value)
        expected <<= self.q.unit ** self.q.shape[-1]
        assert_array_equal(det, expected)
        with pytest.raises(np.linalg.LinAlgError):
            np.linalg.det(self.q[0])  # Not 2-D
        with pytest.raises(np.linalg.LinAlgError):
            np.linalg.det(self.q[:-1])  # Not square.

    @needs_array_function
    def test_slogdet(self):
        # TODO: Could be supported if we had a natural logarithm unit.
        with pytest.raises(TypeError):
            logdet = np.linalg.slogdet(self.q)
            assert hasattr(logdet, 'unit')

    @needs_array_function
    def test_solve(self):
        b = np.array([1., 2., 4.]) * u.m / u.s
        x = np.linalg.solve(self.q, b)
        xx = np.linalg.solve(self.q.value, b.value)
        xx <<= b.unit / self.q.unit
        assert_array_equal(x, xx)
        assert u.allclose(self.q @ x, b)

    @needs_array_function
    def test_tensorsolve(self):
        b = np.array([1., 2., 4.]) * u.m / u.s
        x = np.linalg.tensorsolve(self.q, b)
        xx = np.linalg.tensorsolve(self.q.value, b.value)
        xx <<= b.unit / self.q.unit
        assert_array_equal(x, xx)
        assert u.allclose(self.q @ x, b)

    @needs_array_function
    def test_lstsq(self):
        b = np.array([1., 2., 4.]) * u.m / u.s
        x, residuals, rank, s = np.linalg.lstsq(self.q, b, rcond=None)
        xx, residualsx, rankx, sx = np.linalg.lstsq(self.q.value, b.value,
                                                    rcond=None)
        xx <<= b.unit / self.q.unit
        residualsx <<= b.unit ** 2
        sx <<= self.q.unit

        assert_array_equal(x, xx)
        assert_array_equal(residuals, residualsx)
        assert_array_equal(s, sx)
        assert rank == rankx
        assert u.allclose(self.q @ x, b)

        # Also do one where we can check the answer...
        m = np.eye(3)
        b = np.arange(3) * u.m
        x, residuals, rank, s = np.linalg.lstsq(m, b, rcond=1.*u.percent)
        assert_array_equal(x, b)
        assert np.all(residuals == 0 * u.m**2)
        assert rank == 3
        assert_array_equal(s, np.array([1., 1., 1.]) << u.one)

        with pytest.raises(u.UnitsError):
            np.linalg.lstsq(m, b, rcond=1.*u.s)

    @needs_array_function
    def test_norm(self):
        n = np.linalg.norm(self.q)
        expected = np.linalg.norm(self.q.value) << self.q.unit
        assert_array_equal(n, expected)
        # Special case: 1-D, ord=0.
        n1 = np.linalg.norm(self.q[0], ord=0)
        expected1 = np.linalg.norm(self.q[0].value, ord=0) << u.one
        assert_array_equal(n1, expected1)

    @needs_array_function
    def test_cholesky(self):
        # Numbers from np.linalg.cholesky docstring.
        q = np.array([[1, -2j], [2j, 5]]) * u.m
        cd = np.linalg.cholesky(q)
        cdx = np.linalg.cholesky(q.value) << q.unit ** 0.5
        assert_array_equal(cd, cdx)
        assert u.allclose(cd @ cd.T.conj(), q)

    @needs_array_function
    def test_qr(self):
        # This is not exhaustive...
        a = np.array([[1, -2j], [2j, 5]]) * u.m
        q, r = np.linalg.qr(a)
        qx, rx = np.linalg.qr(a.value)
        qx <<= u.one
        rx <<= a.unit
        assert_array_equal(q, qx)
        assert_array_equal(r, rx)
        assert u.allclose(q @ r, a)

    @needs_array_function
    def test_eig(self):
        w, v = np.linalg.eig(self.q)
        wx, vx = np.linalg.eig(self.q.value)
        wx <<= self.q.unit
        vx <<= u.one
        assert_array_equal(w, wx)
        assert_array_equal(v, vx)

        # Comprehensible example
        q = np.diag((1, 2, 3) * u.m)
        w, v = np.linalg.eig(q)
        assert_array_equal(w, np.arange(1, 4) * u.m)
        assert_array_equal(v, np.eye(3))

    @needs_array_function
    def test_eigvals(self):
        w = np.linalg.eigvals(self.q)
        wx = np.linalg.eigvals(self.q.value) << self.q.unit
        assert_array_equal(w, wx)
        # Comprehensible example
        q = np.diag((1, 2, 3) * u.m)
        w = np.linalg.eigvals(q)
        assert_array_equal(w, np.arange(1, 4) * u.m)

    @needs_array_function
    def test_eigh(self):
        w, v = np.linalg.eigh(self.q)
        wx, vx = np.linalg.eigh(self.q.value)
        wx <<= self.q.unit
        vx <<= u.one
        assert_array_equal(w, wx)
        assert_array_equal(v, vx)

    @needs_array_function
    def test_eigvalsh(self):
        w = np.linalg.eigvalsh(self.q)
        wx = np.linalg.eigvalsh(self.q.value) << self.q.unit
        assert_array_equal(w, wx)


untested_functions = set()
if NUMPY_LT_1_20:
    financial_functions = {f for f in all_wrapped_functions.values()
                           if f in np.lib.financial.__dict__.values()}
    untested_functions |= financial_functions

deprecated_functions = {
    np.asscalar
    }
if NUMPY_LT_1_18:
    deprecated_functions |= {np.rank}
else:
    deprecated_functions |= {np.alen}

untested_functions |= deprecated_functions
io_functions = {np.save, np.savez, np.savetxt, np.savez_compressed}
untested_functions |= io_functions

poly_functions = {
    np.poly, np.polyadd, np.polyder, np.polydiv, np.polyfit, np.polyint,
    np.polymul, np.polysub, np.polyval, np.roots, np.vander
    }
untested_functions |= poly_functions


@needs_array_function
def test_testing_completeness():
    assert not CoverageMeta.covered.intersection(untested_functions)
    assert all_wrapped == (CoverageMeta.covered | untested_functions)


class TestFunctionHelpersCompleteness:
    @pytest.mark.parametrize('one, two', itertools.combinations(
        (SUBCLASS_SAFE_FUNCTIONS,
         UNSUPPORTED_FUNCTIONS,
         set(FUNCTION_HELPERS.keys()),
         set(DISPATCHED_FUNCTIONS.keys())), 2))
    def test_no_duplicates(self, one, two):
        assert not one.intersection(two)

    @needs_array_function
    def test_all_included(self):
        included_in_helpers = (SUBCLASS_SAFE_FUNCTIONS |
                               UNSUPPORTED_FUNCTIONS |
                               set(FUNCTION_HELPERS.keys()) |
                               set(DISPATCHED_FUNCTIONS.keys()))
        assert all_wrapped == included_in_helpers

    # untested_function is created using all_wrapped_functions
    @needs_array_function
    def test_ignored_are_untested(self):
        assert IGNORED_FUNCTIONS == untested_functions
