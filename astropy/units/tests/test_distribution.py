# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import absolute_import, division, print_function

import numpy as np

from ... import units as u
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
    # Test that we can initialize directly from a Numpy array if we provide a unit
    rates = np.array([1, 5, 30, 400])[:, np.newaxis]
    parr = np.random.poisson(rates, (4, 1000))
    u.Distribution(parr, unit=u.kpc)

def test_numpy_init_T():
    rates = np.array([1, 5, 30, 400])
    parr = np.random.poisson(rates, (1000, 4))
    u.Distribution(parr.T, unit=u.kpc)


def test_quantity_init():
    # Test that we can initialize directly from a Quantity
    pq = np.random.poisson(np.array([1, 5, 30, 400])[:, np.newaxis],
                           (4, 1000)) * u.kpc
    u.Distribution(pq)


def test_init_scalar():
    parr = np.random.poisson(np.array([1, 5, 30, 400])[:, np.newaxis],
                             (4, 1000))
    with pytest.raises(TypeError) as exc:
        u.Distribution(parr.ravel()[0])
    assert exc.value.args[0] == "Attempted to initialize a Distribution with a scalar"


class TestDistributionStatistics():
    def setup_class(self):
        with NumpyRNGContext(12345):
            self.data = np.random.normal(np.array([1, 2, 3, 4])[:, np.newaxis],
                                         np.array([3, 2, 4, 5])[:, np.newaxis],
                                         (4, 10000))

        self.distr = u.Distribution(self.data * u.kpc)

    def test_shape(self):
        # u.Distribution shape
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
        # Shape of the PDF (note, this is actually the number of values regardless of samples, needs a better name?)
        assert self.distr.shape == (4,)

    def test_pdf_mean(self):
        # Mean of each PDF
        expected = np.mean(self.data, axis=-1) * self.distr.unit
        assert_quantity_allclose(self.distr.pdf_mean, expected)
        assert_quantity_allclose(self.distr.pdf_mean, [1, 2, 3, 4] * self.distr.unit, rtol=0.05)

        # make sure the right type comes out - should be a Quantity because it's
        # now a summary statistic
        assert not isinstance(self.distr.pdf_mean, u.Distribution)
        assert isinstance(self.distr.pdf_mean, u.Quantity)

    def test_pdf_std(self):
        # Standard deviation of each PDF
        expected = np.std(self.data, axis=-1) * self.distr.unit
        assert_quantity_allclose(self.distr.pdf_std, expected)
        assert_quantity_allclose(self.distr.pdf_std, [3, 2, 4, 5] * self.distr.unit, rtol=0.05)

        # make sure the right type comes out - should be a Quantity because it's
        # now a summary statistic
        assert not isinstance(self.distr.pdf_std, u.Distribution)
        assert isinstance(self.distr.pdf_std, u.Quantity)

    def test_pdf_var(self):
        # Variance of each PDF
        expected = np.var(self.data, axis=-1) * self.distr.unit**2
        assert_quantity_allclose(self.distr.pdf_var, expected)
        assert_quantity_allclose(self.distr.pdf_var, [9, 4, 16, 25] * self.distr.unit**2, rtol=0.1)

        # make sure the right type comes out - should be a Quantity because it's
        # now a summary statistic
        assert not isinstance(self.distr.pdf_var, u.Distribution)
        assert isinstance(self.distr.pdf_var, u.Quantity)

    def test_pdf_median(self):
        # Median of each PDF
        expected = np.median(self.data, axis=-1) * self.distr.unit
        assert_quantity_allclose(self.distr.pdf_median, expected)
        assert_quantity_allclose(self.distr.pdf_median, [1, 2, 3, 4] * self.distr.unit, rtol=0.1)

        # make sure the right type comes out - should be a Quantity because it's
        # now a summary statistic
        assert not isinstance(self.distr.pdf_median, u.Distribution)
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
        assert not isinstance(self.distr.pdf_mad, u.Distribution)
        assert isinstance(self.distr.pdf_mad, u.Quantity)
        assert not isinstance(self.distr.pdf_smad, u.Distribution)
        assert isinstance(self.distr.pdf_smad, u.Quantity)

    def test_percentile(self):
        expected = np.percentile(self.data, [10, 50, 90], axis=-1) * self.distr.unit
        percs = self.distr.percentiles([10, 50, 90])
        assert_quantity_allclose(percs, expected)
        assert percs.shape == (3, 4)

        # make sure the right type comes out - should be a Quantity because it's
        # now a summary statistic
        assert not isinstance(percs, u.Distribution)
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
        another_distr = u.Distribution(another_data, unit=u.pc)
        combined_distr = self.distr + another_distr

        expected = np.median(self.data + another_data/1000,
                             axis=-1) * self.distr.unit
        assert_quantity_allclose(combined_distr.pdf_median, expected)

        expected = np.var(self.data + another_data/1000, axis=-1) * self.distr.unit**2
        assert_quantity_allclose(combined_distr.pdf_var, expected)


def test_helper_normal_samples():
    centerq = [1, 5, 30, 400] * u.kpc

    with NumpyRNGContext(12345):
        n_dist = u.NormalDistribution(centerq, std=[0.2, 1.5, 4, 1]*u.kpc, n_samples=100)
        assert n_dist.distribution.shape == (4, 100)
        assert n_dist.shape == (4, )
        assert n_dist.unit == u.kpc
        assert np.all(n_dist.pdf_std > 100*u.pc)

        n_dist2 = u.NormalDistribution(centerq, std=[0.2, 1.5, 4, 1]*u.pc, n_samples=20000)
        assert n_dist2.distribution.shape == (4, 20000)
        assert n_dist2.shape == (4, )
        assert n_dist2.unit == u.kpc
        assert np.all(n_dist2.pdf_std < 100*u.pc)


def test_helper_poisson_samples():
    centerqadu = [1, 5, 30, 400] * u.adu

    with NumpyRNGContext(12345):
        p_dist = u.PoissonDistribution(centerqadu, n_samples=100)
        assert p_dist.shape == (4,)
        assert p_dist.distribution.shape == (4, 100)
        assert p_dist.unit == u.adu
        p_min = np.min(p_dist)
        assert isinstance(p_min, u.Distribution)
        assert p_min.shape == ()
        assert np.all(p_min >= 0)
        assert np.all(np.abs(p_dist.pdf_mean - centerqadu) < centerqadu)


def test_helper_uniform_samples():
    udist = u.UniformDistribution([1, 2]*u.kpc, [3, 4]*u.kpc)
    assert udist.shape == (2, )
    assert udist.distribution.shape == (2, 1000)
    assert np.all(np.min(udist.distribution, axis=-1) > [1, 2]*u.kpc)
    assert np.all(np.max(udist.distribution, axis=-1) < [3, 4]*u.kpc)

    # try the alternative creator
    udist = u.UniformDistribution.from_center_width([1, 3, 2] * u.pc, [5, 4, 3] * u.pc)
    assert udist.shape == (3, )
    assert udist.distribution.shape == (3, 1000)
    assert np.all(np.min(udist.distribution, axis=-1) > [-1.5, 1, 0.5]*u.pc)
    assert np.all(np.max(udist.distribution, axis=-1) < [3.5, 5, 3.5]*u.pc)


def test_helper_normal_exact():
    pytest.skip('distribution stretch goal not yet implemented')
    centerq = [1, 5, 30, 400] * u.kpc
    u.NormalDistribution(centerq, std=[0.2, 1.5, 4, 1]*u.kpc)
    u.NormalDistribution(centerq, var=[0.04, 2.25, 16, 1]*u.kpc**2)
    u.NormalDistribution(centerq, ivar=[25, 0.44444444, 0.625, 1]*u.kpc**-2)


def test_helper_poisson_exact():
    pytest.skip('distribution stretch goal not yet implemented')
    # TODO: what does Poisson u.Distribution mean for a quantity in e.g. kpc?
    # Should we not restrict to dimensionless? Have updated the test to reflect
    # that here.
    centerq = [1, 5, 30, 400] * u.one
    u.PoissonDistribution(centerq)

    with pytest.raises(u.UnitsError) as exc:
        centerq = [1, 5, 30, 400] * u.kpc
        u.PoissonDistribution(centerq)
    assert exc.value.args[0] == ("Poisson distribution can only be computed "
                                 "for dimensionless quantities")


def test_arithmetic_exact():
    pytest.skip('distribution stretch goal not yet implemented')
    dist = (u.NormalDistribution(3 * u.kpc)
            * u.PoissonDistribution(5 * u.one)
            + u.UniformDistribution(3 * u.pc, 5 * u.pc))


def test_reprs():
    darr = np.arange(30).reshape(3, 10)
    distr = u.Distribution(darr, unit=u.kpc)

    assert 'n_samples=10' in repr(distr)
    assert 'n_samples=10' in str(distr)

    assert r'n_{\rm samp}=10' in distr._repr_latex_()


@pytest.mark.parametrize("klass, kws", [
    (u.NormalDistribution, {'center': 0, 'std': 2}),
    (u.UniformDistribution, {'lower': 0, 'upper': 2}),
    (u.PoissonDistribution, {'poissonval': 2}),
    (u.NormalDistribution, {'center': 0*u.count, 'std': 2*u.count}),
    (u.UniformDistribution, {'lower': 0*u.count, 'upper': 2*u.count}),
    (u.PoissonDistribution, {'poissonval': 2*u.count})
])
def test_wrong_kw_fails(klass, kws):
    with pytest.raises(Exception):
        kw_temp = kws.copy()
        kw_temp['n_sample'] = 100  # note the missing "s"
        assert klass(**kw_temp).n_samples == 100
    kw_temp = kws.copy()
    kw_temp['n_samples'] = 100
    assert klass(**kw_temp).n_samples == 100


def test_index_assignment_quantity():
    arr = np.random.randn(2, 1000)
    distr = u.Distribution(arr*u.kpc)
    d1q, d2q = distr
    assert isinstance(d1q, u.Distribution)
    assert isinstance(d2q, u.Distribution)

    ndistr = u.NormalDistribution(center=[1, 2]*u.kpc, std=[3, 4]*u.kpc)
    n1, n2 = ndistr
    assert isinstance(n1, u.NormalDistribution)
    assert isinstance(n2, u.NormalDistribution)


def test_index_assignment_array():
    arr = np.random.randn(2, 1000)
    distr = u.Distribution(arr)
    d1a, d2a = distr
    assert isinstance(d1a, u.Distribution)
    assert isinstance(d2a, u.Distribution)

    ndistr = u.NormalDistribution(center=[1, 2], std=[3, 4])
    n1, n2 = ndistr
    assert isinstance(n1, u.NormalDistribution)
    assert isinstance(n2, u.NormalDistribution)
