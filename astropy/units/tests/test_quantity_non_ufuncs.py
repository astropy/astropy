import inspect

import numpy as np

import pytest

from astropy import units as u
from astropy.utils.compat import NUMPY_LT_1_17


def get_covered_functions(cls):
    """Helper function to extract the numpy functions covered.

    Assumes that a test is called 'test_<function_name>'.
    """
    covered = []
    for k, v in cls.__dict__.items():
        if inspect.isfunction(v) and k.startswith('test'):
            f = k.replace('test_', '')
            if f in all_wrapped_functions:
                covered.append(all_wrapped_functions[f])

    return set(covered)


# To get the functions that could be covered, we look for those that
# are wrapped.  Of course, this does not give a full list pre-1.17.
all_wrapped_functions = {name: f for name, f in np.__dict__.items()
                         if callable(f) and hasattr(f, '__wrapped__') and
                         f is not np.printoptions}
check = set()


class BasicTestSetup:
    def setup(self):
        self.q = np.arange(9.).reshape(3, 3) * u.m

    def check(self, func, *args, **kwargs):
        o = func(self.q, *args, **kwargs)
        expected = func(self.q.value, *args, **kwargs) * self.q.unit
        assert o.shape == expected.shape
        assert np.all(o == expected)


class TestShapeInformation(BasicTestSetup):
    def test_alen(self):
        assert np.alen(self.q) == 3

    def test_shape(self):
        assert np.shape(self.q) == (3, 3)

    def test_size(self):
        assert np.size(self.q) == 9

    def test_ndim(self):
        assert np.ndim(self.q) == 2


check |= get_covered_functions(TestShapeInformation)


class TestShapeManipulation(BasicTestSetup):
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

    @pytest.mark.xfail
    def test_broadcast_to(self):
        self.check(np.broadcast_to, (3, 3, 3))

    @pytest.mark.xfail
    def test_broadcast_arrays(self):
        q2 = np.ones(3, 3, 3) / u.s
        o1, o2 = np.broadcast_arrays(self.q, q2)
        assert isinstance(o1, u.Quantity)
        assert isinstance(o2, u.Quantity)
        assert o1.shape == o2.shape == (3, 3, 3)
        assert np.all(o1 == self.q)
        assert np.all(o2 == q2)


check |= get_covered_functions(TestShapeManipulation)

shape_indexing = {
    np.diag_indices_from, np.tril_indices_from, np.triu_indices_from
    }
check |= shape_indexing


class TestRealImag(BasicTestSetup):
    def setup(self):
        self.q = (np.arange(9.).reshape(3, 3) + 1j) * u.m

    def test_real(self):
        self.check(np.real)

    def test_imag(self):
        self.check(np.imag)


check |= get_covered_functions(TestRealImag)


class TestCopyAndCreation(BasicTestSetup):
    @pytest.mark.xfail
    def test_copy(self):
        self.check(np.copy)

    @pytest.mark.xfail
    def test_asfarray(self):
        self.check(np.asfarray)

    def test_empty_like(self):
        o = np.empty_like(self.q)
        assert o.shape == (3, 3)
        assert isinstance(o, u.Quantity)
        assert o.unit == self.q.unit

    def test_zeros_like(self):
        self.check(np.zeros_like)

    def test_ones_like(self):
        self.check(np.ones_like)

    @pytest.mark.xfail
    def test_full_like(self):
        o = np.full_like(self.q, 0.5 * u.km)
        expected = np.empty_like(self.q.value) * u.m
        expected[...] = 0.5 * u.km
        assert np.all(o == expected)


check |= get_covered_functions(TestCopyAndCreation)


class TestAccessingParts(BasicTestSetup):
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

    @pytest.mark.xfail
    def test_place(self):
        q = np.arange(3.) * u.m
        np.place(q, [True, False, True], [50, 150] * u.cm)
        assert q.unit == u.m
        expected = [50, 100, 150] * u.cm
        assert np.all(q == expected)


check |= get_covered_functions(TestAccessingParts)

joining_and_splitting = {
    np.concatenate, np.stack, np.column_stack, np.dstack, np.hstack,
    np.vstack, np.block,
    np.split, np.array_split, np.dsplit, np.hsplit, np.vsplit,
    np.tile, np.repeat, np.pad}
check |= joining_and_splitting

accessing_elements = {
    np.delete, np.insert, np.append, np.resize, np.trim_zeros, np.flatnonzero,
    np.roll, np.put, np.putmask, np.take, np.fill_diagonal,
    }
check |= accessing_elements

ufunc_like = {
    np.angle, np.around, np.round_, np.ptp, np.fix, np.i0,
    np.clip, np.sinc, np.where, np.choose, np.select,
    np.isneginf, np.isposinf, np.isclose, np.nan_to_num,
    np.isreal, np.iscomplex, np.real_if_close,
    np.tril, np.triu, np.unwrap, np.copyto,
}
check |= ufunc_like

ufunc_reductions = {
    np.any, np.all, np.amax, np.amin, np.sum, np.prod, np.product,
    np.alltrue, np.sometrue
    }
check |= ufunc_reductions

ufunc_nanreductions = {
    np.nanmax, np.nanmin, np.nanmean, np.nanmedian, np.nansum, np.nanprod,
    np.nanpercentile, np.nanquantile, np.nanstd, np.nanvar,
    }
check |= ufunc_nanreductions


class TestReductionLikeFunctions(BasicTestSetup):
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

    @pytest.mark.xfail
    def test_quantile(self):
        self.check(np.quantile, 0.5)
        o = np.quantile(self.q, 50 * u.percent)
        expected = np.quantile(self.q.value, 0.5) * u.m
        assert np.all(o == expected)

    @pytest.mark.xfail
    def test_percentile(self):
        self.check(np.percentile, 0.5)
        o = np.percentile(self.q, 0.5 * u.one)
        expected = np.percentile(self.q.value, 50) * u.m
        assert np.all(o == expected)

    def test_trace(self):
        self.check(np.trace)

    @pytest.mark.xfail
    def test_count_nonzero(self):
        q1 = np.arange(9.).reshape(3, 3) * u.m
        o = np.count_nonzero(q1)
        assert o == 8
        o = np.count_nonzero(q1, axis=1)
        # Returns integer Quantity with units of m
        assert o == np.array([2, 3, 3])

    @pytest.mark.xfail
    def test_allclose(self):
        q1 = np.arange(3.) * u.m
        q2 = np.array([0., 101., 199.]) * u.cm
        atol = 2 * u.cm
        rtol = 1. * u.percent
        assert np.allclose(q1, q2, atol=atol)
        # Default atol breaks code; everything else works.
        assert np.allclose(q1, q2, rtol=rtol)

    def test_allclose_failures(self):
        q1 = np.arange(3.) * u.m
        q2 = np.array([0., 101., 199.]) * u.cm
        with pytest.raises(u.UnitsError):
            np.allclose(q1, q2, atol=2, rtol=0)
        with pytest.raises(u.UnitsError):
            np.allclose(q1, q2, atol=0, rtol=1. * u.s)

    @pytest.mark.xfail
    def test_array_equal(self):
        q1 = np.arange(3.) * u.m
        q2 = q1.to(u.cm)
        assert np.array_equal(q1, q2)
        q3 = q1.value * u.cm
        assert not np.array_equal(q1, q3)

    @pytest.mark.xfail
    def test_array_equiv(self):
        q1 = np.array([[0., 1., 2.]]*3) * u.m
        q2 = q1[0].to(u.cm)
        assert np.array_equiv(q1, q2)
        q3 = q1[0].value * u.cm
        assert not np.array_equiv(q1, q3)


check |= get_covered_functions(TestReductionLikeFunctions)


ufunc_accumulations = {
    np.cumsum, np.cumprod, np.cumproduct
    }
check |= ufunc_accumulations

ufunc_accumulation_like = {
    np.nancumsum, np.nancumprod
    }
check |= ufunc_accumulation_like


class TestVariousProductFunctions:
    """
    Test functions that are similar to gufuncs
    """
    @pytest.mark.xfail
    def test_cross(self):
        q1 = np.arange(6.).reshape(2, 3) * u.m
        q2 = np.array([4., 5., 6.]) / u.s
        o = np.cross(q1, q2)
        expected = np.cross(q1.value, q2.value) * u.m / u.s
        assert np.all(o == expected)

    @pytest.mark.xfail
    def test_outer(self):
        q1 = np.array([1, 2, 3]) * u.m
        q2 = np.array([1, 2]) / u.s
        o = np.outer(q1, q2)
        assert np.all(o == np.array([[1, 2], [2, 4], [3, 6]]) * u.m / u.s)

    @pytest.mark.xfail
    def test_inner(self):
        q1 = np.array([1, 2, 3]) * u.m
        q2 = np.array([4, 5, 6]) / u.s
        o = np.inner(q1, q2)
        assert o == 32 * u.m / u.s

    @pytest.mark.xfail
    def test_dot(self):
        q1 = np.array([1., 2., 3.]) * u.m
        q2 = np.array([4., 5., 6.]) / u.s
        o = np.dot(q1, q2)
        assert o == 32. * u.m / u.s

    @pytest.mark.xfail
    def test_vdot(self):
        q1 = np.array([1j, 2j, 3j]) * u.m
        q2 = np.array([4j, 5j, 6j]) / u.s
        o = np.vdot(q1, q2)
        assert o == (32. + 0j) * u.m / u.s

    @pytest.mark.xfail
    def test_tensordot(self):
        # From the docstring example
        a = np.arange(60.).reshape(3, 4, 5) * u.m
        b = np.arange(24.).reshape(4, 3, 2) / u.s
        c = np.tensordot(a, b, axes=([1, 0], [0, 1]))
        expected = np.tensordot(a.value, b.value,
                                axes=([1, 0], [0, 1])) * u.m / u.s
        assert np.all(c == expected)

    @pytest.mark.xfail
    def test_kron(self):
        q1 = np.eye(2) * u.m
        q2 = np.ones(2) / u.s
        o = np.kron(q1, q2)
        expected = np.kron(q1.value, q2.value) * u.m / u.s
        assert np.all(o == expected)

    @pytest.mark.xfail
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


check |= get_covered_functions(TestVariousProductFunctions)

gufunc_like = {
    np.interp,
    np.linspace, np.logspace, np.geomspace,
    np.bincount, np.digitize, np.histogram, np.histogram2d, np.histogramdd,
    np.histogram_bin_edges, np.diff, np.gradient,
    np.cov, np.corrcoef, np.piecewise, np.convolve, np.correlate,
    np.sort, np.sort_complex, np.lexsort, np.msort, np.partition, np.trapz,
    np.searchsorted
    }
check |= gufunc_like

arg_functions = {
    np.argmax, np.argmin, np.argpartition, np.argsort, np.argwhere,
    np.nonzero, np.nanargmin, np.nanargmax,
    np.take_along_axis, np.put_along_axis
    }
check |= arg_functions

string_functions = {np.array_repr, np.array_str, np.array2string}
check |= string_functions

bit_functions = {np.packbits, np.unpackbits}
check |= bit_functions

index_functions = {np.ravel_multi_index, np.unravel_index}
check |= index_functions

dtype_functions = {
    np.common_type, np.result_type, np.can_cast, np.min_scalar_type,
    np.iscomplexobj, np.isrealobj,
    }
check |= dtype_functions

mesh_functions = {np.ix_, np.meshgrid}
check |= mesh_functions

function_functions = {
    np.apply_along_axis, np.apply_over_axes,
    }
check |= function_functions

financial_functions = {f for f in all_wrapped_functions.values()
                       if f in np.lib.financial.__dict__.values()}
check |= financial_functions

datetime_functions = {
    np.busday_count, np.busday_offset, np.datetime_as_string,
    np.is_busday,
    }
check |= datetime_functions.copy()

deprecated_functions = {
    np.asscalar, np.rank
    }
check |= deprecated_functions

memory_functions = {
    np.shares_memory, np.may_share_memory
    }
check |= memory_functions

io_functions = {np.save, np.savez, np.savetxt, np.savez_compressed}
check |= io_functions

poly_functions = {
    np.poly, np.polyadd, np.polyder, np.polydiv, np.polyfit, np.polyint,
    np.polymul, np.polysub, np.polyval, np.roots, np.vander
    }
check |= poly_functions

setops_functions = {f for f in all_wrapped_functions.values()
                    if f in np.lib.arraysetops.__dict__.values()}
check |= setops_functions


@pytest.mark.xfail(NUMPY_LT_1_17,
                   reason="no __array_function__ wrapping in numpy<1.17")
def test_completeness():
    assert set(all_wrapped_functions.values()) == check
