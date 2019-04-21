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

copying_and_creation = {
    np.copy, np.asfarray,
    np.empty_like, np.ones_like, np.zeros_like, np.full_like,
    }
check |= copying_and_creation


class TestShapeInformation:
    def setup(self):
        self.q = np.arange(9.).reshape(3, 3) * u.m

    def test_alen(self):
        assert np.alen(self.q) == 3

    def test_shape(self):
        assert np.shape(self.q) == (3, 3)

    def test_size(self):
        assert np.size(self.q) == 9

    def test_ndim(self):
        assert np.ndim(self.q) == 2


check |= get_covered_functions(TestShapeInformation)


shape_manipulation = {
    np.reshape, np.ravel,
    np.moveaxis, np.rollaxis, np.swapaxes, np.transpose,
    np.atleast_1d, np.atleast_2d, np.atleast_3d,
    np.expand_dims, np.squeeze,
    np.flip, np.fliplr, np.flipud, np.rot90}
check |= shape_manipulation

shape_indexing = {
    np.diag_indices_from, np.tril_indices_from, np.triu_indices_from
    }
check |= shape_indexing


class TestRealImag:
    def setup(self):
        self.q = (np.arange(9.).reshape(3, 3) + 1j) * u.m

    def test_real(self):
        o = np.real(self.q)
        expected = np.real(self.q.value) * self.q.unit
        assert np.all(o == expected)

    def test_imag(self):
        o = np.imag(self.q)
        expected = np.imag(self.q.value) * self.q.unit
        assert np.all(o == expected)


check |= get_covered_functions(TestRealImag)

different_views = {
    np.diag, np.diagonal, np.compress, np.extract, np.diagflat, np.place
    }
check |= different_views

broadcast_functions = {
    np.broadcast_to, np.broadcast_arrays,
    }
check |= broadcast_functions

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


class TestReductionLikeFunctions:
    def test_average(self):
        q1 = np.arange(9.).reshape(3, 3) * u.m
        q2 = np.eye(3) / u.s
        o = np.average(q1, weights=q2)
        expected = np.average(q1.value, weights=q2.value) * u.m
        assert np.all(o == expected)

    def test_mean(self):
        q1 = np.arange(9.).reshape(3, 3) * u.m
        o = np.mean(q1)
        expected = np.mean(q1.value) * u.m
        assert np.all(o == expected)

    def test_std(self):
        q1 = np.arange(9.).reshape(3, 3) * u.m
        o = np.std(q1)
        expected = np.std(q1.value) * u.m
        assert np.all(o == expected)

    def test_var(self):
        q1 = np.arange(9.).reshape(3, 3) * u.m
        o = np.var(q1)
        expected = np.var(q1.value) * u.m ** 2
        assert np.all(o == expected)

    def test_median(self):
        q1 = np.arange(9.).reshape(3, 3) * u.m
        o = np.median(q1)
        expected = np.median(q1.value) * u.m
        assert np.all(o == expected)

    @pytest.mark.xfail
    def test_quantile(self):
        q1 = np.arange(9.).reshape(3, 3) * u.m
        o = np.quantile(q1, 0.5)
        expected = np.quantile(q1.value, 0.5) * u.m
        assert np.all(o == expected)

    @pytest.mark.xfail
    def test_percentile(self):
        q1 = np.arange(9.).reshape(3, 3) * u.m
        o = np.percentile(q1, 50)
        expected = np.percentile(q1.value, 50) * u.m
        assert np.all(o == expected)

    def test_trace(self):
        q1 = np.arange(9.).reshape(3, 3) * u.m
        o = np.trace(q1)
        expected = np.trace(q1.value) * u.m
        assert np.all(o == expected)

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
