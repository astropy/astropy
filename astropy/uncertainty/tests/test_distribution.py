# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import absolute_import, division, print_function

import numpy as np

from ... import units as u
from ..core import Distribution
from .. import distributions as ds
from ...utils import NumpyRNGContext
from ...tests.helper import assert_quantity_allclose, pytest

try:
    from scipy.stats import norm  # pylint: disable=W0611
    SMAD_FACTOR = 1 / norm.ppf(0.75)
except ImportError:
    HAS_SCIPY = False
else:
    HAS_SCIPY = True


def test_numpy_init():
    # Test that we can initialize directly from a Numpy array
    rates = np.array([1, 5, 30, 400])[:, np.newaxis]
    parr = np.random.poisson(rates, (4, 1000))
    Distribution(parr)


def test_numpy_init_T():
    rates = np.array([1, 5, 30, 400])
    parr = np.random.poisson(rates, (1000, 4))
    Distribution(parr.T)


def test_quantity_init():
    # Test that we can initialize directly from a Quantity
    pq = np.random.poisson(np.array([1, 5, 30, 400])[:, np.newaxis],
                           (4, 1000)) * u.ct
    Distribution(pq)


def test_quantity_init_T():
    # Test that we can initialize directly from a Quantity
    pq = np.random.poisson(np.array([1, 5, 30, 400]), (1000, 4)) * u.ct
    Distribution(pq.T)


def test_init_scalar():
    parr = np.random.poisson(np.array([1, 5, 30, 400])[:, np.newaxis],
                             (4, 1000))
    with pytest.raises(TypeError) as exc:
        Distribution(parr.ravel()[0])
    assert exc.value.args[0] == "Attempted to initialize a Distribution with a scalar"


class TestDistributionStatistics():
    def setup_class(self):
        with NumpyRNGContext(12345):
            self.data = np.random.normal(np.array([1, 2, 3, 4])[:, np.newaxis],
                                         np.array([3, 2, 4, 5])[:, np.newaxis],
                                         (4, 10000))

        self.distr = Distribution(self.data * u.kpc)

    def test_shape(self):
        # Distribution shape
        assert self.distr.shape == (4, )
        assert self.distr.distribution.shape == (4, 10000)

    def test_size(self):
        # Total number of values
        assert self.distr.size == 4
        assert self.distr.distribution.size == 40000

    def test_n_samples(self):
        # Number of samples
        assert self.distr.n_samples == 10000

    def test_n_distr(self):
        assert self.distr.shape == (4,)

    def test_pdf_mean(self):
        # Mean of each PDF
        expected = np.mean(self.data, axis=-1) * self.distr.unit
        assert_quantity_allclose(self.distr.pdf_mean, expected)
        assert_quantity_allclose(self.distr.pdf_mean, [1, 2, 3, 4] * self.distr.unit, rtol=0.05)

        # make sure the right type comes out - should be a Quantity because it's
        # now a summary statistic
        assert not isinstance(self.distr.pdf_mean, Distribution)
        assert isinstance(self.distr.pdf_mean, u.Quantity)

    def test_pdf_std(self):
        # Standard deviation of each PDF
        expected = np.std(self.data, axis=-1) * self.distr.unit
        assert_quantity_allclose(self.distr.pdf_std, expected)
        assert_quantity_allclose(self.distr.pdf_std, [3, 2, 4, 5] * self.distr.unit, rtol=0.05)

        # make sure the right type comes out - should be a Quantity because it's
        # now a summary statistic
        assert not isinstance(self.distr.pdf_std, Distribution)
        assert isinstance(self.distr.pdf_std, u.Quantity)

    def test_pdf_var(self):
        # Variance of each PDF
        expected = np.var(self.data, axis=-1) * self.distr.unit**2
        assert_quantity_allclose(self.distr.pdf_var, expected)
        assert_quantity_allclose(self.distr.pdf_var, [9, 4, 16, 25] * self.distr.unit**2, rtol=0.1)

        # make sure the right type comes out - should be a Quantity because it's
        # now a summary statistic
        assert not isinstance(self.distr.pdf_var, Distribution)
        assert isinstance(self.distr.pdf_var, u.Quantity)

    def test_pdf_median(self):
        # Median of each PDF
        expected = np.median(self.data, axis=-1) * self.distr.unit
        assert_quantity_allclose(self.distr.pdf_median, expected)
        assert_quantity_allclose(self.distr.pdf_median, [1, 2, 3, 4] * self.distr.unit, rtol=0.1)

        # make sure the right type comes out - should be a Quantity because it's
        # now a summary statistic
        assert not isinstance(self.distr.pdf_median, Distribution)
        assert isinstance(self.distr.pdf_median, u.Quantity)

    @pytest.mark.skipif(not HAS_SCIPY, reason='no scipy')
    def test_pdf_mad_smad(self):
        # Median absolute deviation of each PDF
        median = np.median(self.data, axis=-1, keepdims=True)
        expected = np.median(np.abs(self.data - median), axis=-1) * self.distr.unit
        assert_quantity_allclose(self.distr.pdf_mad, expected)
        assert_quantity_allclose(self.distr.pdf_smad, self.distr.pdf_mad * SMAD_FACTOR, rtol=1e-5)
        assert_quantity_allclose(self.distr.pdf_smad, [3, 2, 4, 5] * self.distr.unit, rtol=0.05)

        # make sure the right type comes out - should be a Quantity because it's
        # now a summary statistic
        assert not isinstance(self.distr.pdf_mad, Distribution)
        assert isinstance(self.distr.pdf_mad, u.Quantity)
        assert not isinstance(self.distr.pdf_smad, Distribution)
        assert isinstance(self.distr.pdf_smad, u.Quantity)

    def test_percentile(self):
        expected = np.percentile(self.data, [10, 50, 90], axis=-1) * self.distr.unit
        percs = self.distr.pdf_percentiles([10, 50, 90])
        assert_quantity_allclose(percs, expected)
        assert percs.shape == (3, 4)

        # make sure the right type comes out - should be a Quantity because it's
        # now a summary statistic
        assert not isinstance(percs, Distribution)
        assert isinstance(percs, u.Quantity)

    def test_add_quantity(self):
        distrplus = self.distr + [2000, 0, 0, 500] * u.pc
        expected = (np.median(self.data, axis=-1) + np.array([2, 0, 0, 0.5])) * self.distr.unit
        assert_quantity_allclose(distrplus.pdf_median, expected)
        expected = np.var(self.data, axis=-1) * self.distr.unit**2
        assert_quantity_allclose(distrplus.pdf_var, expected)

    def test_add_distribution(self):
        another_data = (np.random.randn(4, 10000)
                        * np.array([1000, .01, 80, 10])[:, np.newaxis]
                        + np.array([2000, 0, 0, 500])[:, np.newaxis])
        # another_data is in pc, but main distr is in kpc
        another_distr = Distribution(another_data * u.pc)
        combined_distr = self.distr + another_distr

        expected = np.median(self.data + another_data/1000,
                             axis=-1) * self.distr.unit
        assert_quantity_allclose(combined_distr.pdf_median, expected)

        expected = np.var(self.data + another_data/1000, axis=-1) * self.distr.unit**2
        assert_quantity_allclose(combined_distr.pdf_var, expected)


def test_helper_normal_samples():
    centerq = [1, 5, 30, 400] * u.kpc

    with NumpyRNGContext(12345):
        n_dist = ds.normal(centerq, std=[0.2, 1.5, 4, 1]*u.kpc, n_samples=100)
        assert n_dist.distribution.shape == (4, 100)
        assert n_dist.shape == (4, )
        assert n_dist.unit == u.kpc
        assert np.all(n_dist.pdf_std > 100*u.pc)

        n_dist2 = ds.normal(centerq, std=[0.2, 1.5, 4, 1]*u.pc, n_samples=20000)
        assert n_dist2.distribution.shape == (4, 20000)
        assert n_dist2.shape == (4, )
        assert n_dist2.unit == u.kpc
        assert np.all(n_dist2.pdf_std < 100*u.pc)


def test_helper_poisson_samples():
    centerqcounts = [1, 5, 30, 400] * u.count

    with NumpyRNGContext(12345):
        p_dist = ds.poisson(centerqcounts, n_samples=100)
        assert p_dist.shape == (4,)
        assert p_dist.distribution.shape == (4, 100)
        assert p_dist.unit == u.count
        p_min = np.min(p_dist)
        assert isinstance(p_min, Distribution)
        assert p_min.shape == ()
        assert np.all(p_min >= 0)
        assert np.all(np.abs(p_dist.pdf_mean - centerqcounts) < centerqcounts)


def test_helper_uniform_samples():
    udist = ds.uniform(lower=[1, 2]*u.kpc, upper=[3, 4]*u.kpc, n_samples=1000)
    assert udist.shape == (2, )
    assert udist.distribution.shape == (2, 1000)
    assert np.all(np.min(udist.distribution, axis=-1) > [1, 2]*u.kpc)
    assert np.all(np.max(udist.distribution, axis=-1) < [3, 4]*u.kpc)

    # try the alternative creator
    udist = ds.uniform(center=[1, 3, 2] * u.pc, width=[5, 4, 3] * u.pc, n_samples=1000)
    assert udist.shape == (3, )
    assert udist.distribution.shape == (3, 1000)
    assert np.all(np.min(udist.distribution, axis=-1) > [-1.5, 1, 0.5]*u.pc)
    assert np.all(np.max(udist.distribution, axis=-1) < [3.5, 5, 3.5]*u.pc)


def test_helper_normal_exact():
    pytest.skip('distribution stretch goal not yet implemented')
    centerq = [1, 5, 30, 400] * u.kpc
    ds.normal(centerq, std=[0.2, 1.5, 4, 1]*u.kpc)
    ds.normal(centerq, var=[0.04, 2.25, 16, 1]*u.kpc**2)
    ds.normal(centerq, ivar=[25, 0.44444444, 0.625, 1]*u.kpc**-2)


def test_helper_poisson_exact():
    pytest.skip('distribution stretch goal not yet implemented')
    centerq = [1, 5, 30, 400] * u.one
    ds.poisson(centerq, n_samples=1000)

    with pytest.raises(u.UnitsError) as exc:
        centerq = [1, 5, 30, 400] * u.kpc
        ds.poisson(centerq, n_samples=1000)
    assert exc.value.args[0] == ("Poisson distribution can only be computed "
                                 "for dimensionless quantities")


def test_reprs():
    darr = np.arange(30).reshape(3, 10)
    distr = Distribution(darr * u.kpc)

    assert 'n_samples=10' in repr(distr)
    assert 'n_samples=10' in str(distr)

    assert r'n_{\rm samp}=10' in distr._repr_latex_()


@pytest.mark.parametrize("func, kws", [
    (ds.normal, {'center': 0, 'std': 2}),
    (ds.uniform, {'lower': 0, 'upper': 2}),
    (ds.poisson, {'center': 2}),
    (ds.normal, {'center': 0*u.count, 'std': 2*u.count}),
    (ds.uniform, {'lower': 0*u.count, 'upper': 2*u.count}),
    (ds.poisson, {'center': 2*u.count})
])
def test_wrong_kw_fails(func, kws):
    with pytest.raises(Exception):
        kw_temp = kws.copy()
        kw_temp['n_sample'] = 100  # note the missing "s"
        assert func(**kw_temp).n_samples == 100
    kw_temp = kws.copy()
    kw_temp['n_samples'] = 100
    assert func(**kw_temp).n_samples == 100


def test_index_assignment_quantity():
    arr = np.random.randn(2, 1000)
    distr = Distribution(arr*u.kpc)
    d1q, d2q = distr
    assert isinstance(d1q, Distribution)
    assert isinstance(d2q, Distribution)

    ndistr = ds.normal(center=[1, 2]*u.kpc, std=[3, 4]*u.kpc, n_samples=1000)
    n1, n2 = ndistr
    assert isinstance(n1, ds.Distribution)
    assert isinstance(n2, ds.Distribution)


def test_index_assignment_array():
    arr = np.random.randn(2, 1000)
    distr = Distribution(arr)
    d1a, d2a = distr
    assert isinstance(d1a, Distribution)
    assert isinstance(d2a, Distribution)

    ndistr = ds.normal(center=[1, 2], std=[3, 4], n_samples=1000)
    n1, n2 = ndistr
    assert isinstance(n1, ds.Distribution)
    assert isinstance(n2, ds.Distribution)


def test_histogram():
    arr = np.random.randn(2, 3, 1000)
    distr = Distribution(arr)

    hist, bins = distr.pdf_histogram(bins=10)
    assert hist.shape == (2, 3, 10)
    assert bins.shape == (2, 3, 11)


def test_array_repr_latex():
    # as of this writing ndarray does not have a _repr_latex_, and this test
    # ensure distributions account for that. However, if in the future ndarray
    # gets a _repr_latex_, we can skip this.


    arr = np.random.randn(4, 1000)

    if hasattr(arr, '_repr_latex_'):
        pytest.skip('in this version of numpy, ndarray has a _repr_latex_')

    distr = Distribution(arr)
    assert distr._repr_latex_() is None
