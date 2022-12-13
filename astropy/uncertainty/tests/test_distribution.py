# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
import pytest
from numpy.testing import assert_array_equal

from astropy import units as u
from astropy.coordinates import Angle
from astropy.tests.helper import assert_quantity_allclose
from astropy.uncertainty import distributions as ds
from astropy.uncertainty.core import Distribution
from astropy.utils import NumpyRNGContext
from astropy.utils.compat.optional_deps import HAS_SCIPY

if HAS_SCIPY:
    from scipy.stats import norm  # pylint: disable=W0611

    SMAD_FACTOR = 1 / norm.ppf(0.75)


class TestInit:
    @classmethod
    def setup_class(self):
        self.rates = np.array([1, 5, 30, 400])[:, np.newaxis]
        self.parr = np.random.poisson(self.rates, (4, 1000))
        self.parr_t = np.random.poisson(self.rates.squeeze(), (1000, 4))

    def test_numpy_init(self):
        # Test that we can initialize directly from a Numpy array
        Distribution(self.parr)

    def test_numpy_init_T(self):
        Distribution(self.parr_t.T)

    def test_quantity_init(self):
        # Test that we can initialize directly from a Quantity
        pq = self.parr << u.ct
        pqd = Distribution(pq)
        assert isinstance(pqd, u.Quantity)
        assert isinstance(pqd, Distribution)
        assert isinstance(pqd.value, Distribution)
        assert_array_equal(pqd.value.distribution, self.parr)

    def test_quantity_init_T(self):
        # Test that we can initialize directly from a Quantity
        pq = self.parr_t << u.ct
        Distribution(pq.T)

    def test_quantity_init_with_distribution(self):
        # Test that we can initialize a Quantity from a Distribution.
        pd = Distribution(self.parr)
        qpd = pd << u.ct
        assert isinstance(qpd, u.Quantity)
        assert isinstance(qpd, Distribution)
        assert qpd.unit == u.ct
        assert_array_equal(qpd.value.distribution, pd.distribution.astype(float))


def test_init_scalar():
    parr = np.random.poisson(np.array([1, 5, 30, 400])[:, np.newaxis], (4, 1000))
    with pytest.raises(
        TypeError, match=r"Attempted to initialize a Distribution with a scalar"
    ):
        Distribution(parr.ravel()[0])


class TestDistributionStatistics:
    def setup_class(self):
        with NumpyRNGContext(12345):
            self.data = np.random.normal(
                np.array([1, 2, 3, 4])[:, np.newaxis],
                np.array([3, 2, 4, 5])[:, np.newaxis],
                (4, 10000),
            )

        self.distr = Distribution(self.data * u.kpc)

    def test_shape(self):
        # Distribution shape
        assert self.distr.shape == (4,)
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
        pdf_mean = self.distr.pdf_mean()
        assert_quantity_allclose(pdf_mean, expected)
        assert_quantity_allclose(pdf_mean, [1, 2, 3, 4] * self.distr.unit, rtol=0.05)

        # make sure the right type comes out - should be a Quantity because it's
        # now a summary statistic
        assert not isinstance(pdf_mean, Distribution)
        assert isinstance(pdf_mean, u.Quantity)

        # Check with out argument.
        out = pdf_mean * 0.0
        pdf_mean2 = self.distr.pdf_mean(out=out)
        assert pdf_mean2 is out
        assert np.all(pdf_mean2 == pdf_mean)

    def test_pdf_std(self):
        # Standard deviation of each PDF
        expected = np.std(self.data, axis=-1) * self.distr.unit
        pdf_std = self.distr.pdf_std()
        assert_quantity_allclose(pdf_std, expected)
        assert_quantity_allclose(pdf_std, [3, 2, 4, 5] * self.distr.unit, rtol=0.05)

        # make sure the right type comes out - should be a Quantity because it's
        # now a summary statistic
        assert not isinstance(pdf_std, Distribution)
        assert isinstance(pdf_std, u.Quantity)

        # Check with proper ddof, using out argument.
        out = pdf_std * 0.0
        expected = np.std(self.data, axis=-1, ddof=1) * self.distr.unit
        pdf_std2 = self.distr.pdf_std(ddof=1, out=out)
        assert pdf_std2 is out
        assert np.all(pdf_std2 == expected)

    def test_pdf_var(self):
        # Variance of each PDF
        expected = np.var(self.data, axis=-1) * self.distr.unit**2
        pdf_var = self.distr.pdf_var()
        assert_quantity_allclose(pdf_var, expected)
        assert_quantity_allclose(
            pdf_var, [9, 4, 16, 25] * self.distr.unit**2, rtol=0.1
        )

        # make sure the right type comes out - should be a Quantity because it's
        # now a summary statistic
        assert not isinstance(pdf_var, Distribution)
        assert isinstance(pdf_var, u.Quantity)

        # Check with proper ddof, using out argument.
        out = pdf_var * 0.0
        expected = np.var(self.data, axis=-1, ddof=1) * self.distr.unit**2
        pdf_var2 = self.distr.pdf_var(ddof=1, out=out)
        assert pdf_var2 is out
        assert np.all(pdf_var2 == expected)

    def test_pdf_median(self):
        # Median of each PDF
        expected = np.median(self.data, axis=-1) * self.distr.unit
        pdf_median = self.distr.pdf_median()
        assert_quantity_allclose(pdf_median, expected)
        assert_quantity_allclose(pdf_median, [1, 2, 3, 4] * self.distr.unit, rtol=0.1)

        # make sure the right type comes out - should be a Quantity because it's
        # now a summary statistic
        assert not isinstance(pdf_median, Distribution)
        assert isinstance(pdf_median, u.Quantity)

        # Check with out argument.
        out = pdf_median * 0.0
        pdf_median2 = self.distr.pdf_median(out=out)
        assert pdf_median2 is out
        assert np.all(pdf_median2 == expected)

    @pytest.mark.skipif(not HAS_SCIPY, reason="no scipy")
    def test_pdf_mad_smad(self):
        # Median absolute deviation of each PDF
        median = np.median(self.data, axis=-1, keepdims=True)
        expected = np.median(np.abs(self.data - median), axis=-1) * self.distr.unit
        pdf_mad = self.distr.pdf_mad()
        assert_quantity_allclose(pdf_mad, expected)
        pdf_smad = self.distr.pdf_smad()
        assert_quantity_allclose(pdf_smad, pdf_mad * SMAD_FACTOR, rtol=1e-5)
        assert_quantity_allclose(pdf_smad, [3, 2, 4, 5] * self.distr.unit, rtol=0.05)

        # make sure the right type comes out - should be a Quantity because it's
        # now a summary statistic
        assert not isinstance(pdf_mad, Distribution)
        assert isinstance(pdf_mad, u.Quantity)
        assert not isinstance(pdf_smad, Distribution)
        assert isinstance(pdf_smad, u.Quantity)

        # Check out argument for smad (which checks mad too).
        out = pdf_smad * 0.0
        pdf_smad2 = self.distr.pdf_smad(out=out)
        assert pdf_smad2 is out
        assert np.all(pdf_smad2 == pdf_smad)

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
        expected = (
            np.median(self.data, axis=-1) + np.array([2, 0, 0, 0.5])
        ) * self.distr.unit
        assert_quantity_allclose(distrplus.pdf_median(), expected)
        expected = np.var(self.data, axis=-1) * self.distr.unit**2
        assert_quantity_allclose(distrplus.pdf_var(), expected)

    def test_add_distribution(self):
        another_data = (
            np.random.randn(4, 10000) * np.array([1000, 0.01, 80, 10])[:, np.newaxis]
            + np.array([2000, 0, 0, 500])[:, np.newaxis]
        )
        # another_data is in pc, but main distr is in kpc
        another_distr = Distribution(another_data * u.pc)
        combined_distr = self.distr + another_distr

        expected = np.median(self.data + another_data / 1000, axis=-1) * self.distr.unit
        assert_quantity_allclose(combined_distr.pdf_median(), expected)

        expected = (
            np.var(self.data + another_data / 1000, axis=-1) * self.distr.unit**2
        )
        assert_quantity_allclose(combined_distr.pdf_var(), expected)


def test_helper_normal_samples():
    centerq = [1, 5, 30, 400] * u.kpc

    with NumpyRNGContext(12345):
        n_dist = ds.normal(centerq, std=[0.2, 1.5, 4, 1] * u.kpc, n_samples=100)
        assert n_dist.distribution.shape == (4, 100)
        assert n_dist.shape == (4,)
        assert n_dist.unit == u.kpc
        assert np.all(n_dist.pdf_std() > 100 * u.pc)

        n_dist2 = ds.normal(centerq, std=[0.2, 1.5, 4, 1] * u.pc, n_samples=20000)
        assert n_dist2.distribution.shape == (4, 20000)
        assert n_dist2.shape == (4,)
        assert n_dist2.unit == u.kpc
        assert np.all(n_dist2.pdf_std() < 100 * u.pc)


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
        assert np.all(np.abs(p_dist.pdf_mean() - centerqcounts) < centerqcounts)


def test_helper_uniform_samples():
    udist = ds.uniform(lower=[1, 2] * u.kpc, upper=[3, 4] * u.kpc, n_samples=1000)
    assert udist.shape == (2,)
    assert udist.distribution.shape == (2, 1000)
    assert np.all(np.min(udist.distribution, axis=-1) > [1, 2] * u.kpc)
    assert np.all(np.max(udist.distribution, axis=-1) < [3, 4] * u.kpc)

    # try the alternative creator
    udist = ds.uniform(center=[1, 3, 2] * u.pc, width=[5, 4, 3] * u.pc, n_samples=1000)
    assert udist.shape == (3,)
    assert udist.distribution.shape == (3, 1000)
    assert np.all(np.min(udist.distribution, axis=-1) > [-1.5, 1, 0.5] * u.pc)
    assert np.all(np.max(udist.distribution, axis=-1) < [3.5, 5, 3.5] * u.pc)


def test_helper_normal_exact():
    pytest.skip("distribution stretch goal not yet implemented")
    centerq = [1, 5, 30, 400] * u.kpc
    ds.normal(centerq, std=[0.2, 1.5, 4, 1] * u.kpc)
    ds.normal(centerq, var=[0.04, 2.25, 16, 1] * u.kpc**2)
    ds.normal(centerq, ivar=[25, 0.44444444, 0.625, 1] * u.kpc**-2)


def test_helper_poisson_exact():
    pytest.skip("distribution stretch goal not yet implemented")
    centerq = [1, 5, 30, 400] * u.one
    ds.poisson(centerq, n_samples=1000)

    with pytest.raises(
        u.UnitsError,
        match=r"Poisson distribution can only be computed for dimensionless quantities",
    ):
        centerq = [1, 5, 30, 400] * u.kpc
        ds.poisson(centerq, n_samples=1000)


def test_reprs():
    darr = np.arange(30).reshape(3, 10)
    distr = Distribution(darr * u.kpc)

    assert "n_samples=10" in repr(distr)
    assert "n_samples=10" in str(distr)

    assert r"n_{\rm samp}=10" in distr._repr_latex_()


@pytest.mark.parametrize(
    "func, kws",
    [
        (ds.normal, {"center": 0, "std": 2}),
        (ds.uniform, {"lower": 0, "upper": 2}),
        (ds.poisson, {"center": 2}),
        (ds.normal, {"center": 0 * u.count, "std": 2 * u.count}),
        (ds.uniform, {"lower": 0 * u.count, "upper": 2 * u.count}),
        (ds.poisson, {"center": 2 * u.count}),
    ],
)
def test_wrong_kw_fails(func, kws):
    with pytest.raises(Exception):
        kw_temp = kws.copy()
        kw_temp["n_sample"] = 100  # note the missing "s"
        assert func(**kw_temp).n_samples == 100
    kw_temp = kws.copy()
    kw_temp["n_samples"] = 100
    assert func(**kw_temp).n_samples == 100


def test_index_assignment_quantity():
    arr = np.random.randn(2, 1000)
    distr = Distribution(arr * u.kpc)
    d1q, d2q = distr
    assert isinstance(d1q, Distribution)
    assert isinstance(d2q, Distribution)

    ndistr = ds.normal(center=[1, 2] * u.kpc, std=[3, 4] * u.kpc, n_samples=1000)
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

    if hasattr(arr, "_repr_latex_"):
        pytest.skip("in this version of numpy, ndarray has a _repr_latex_")

    distr = Distribution(arr)
    assert distr._repr_latex_() is None


def test_distr_to():
    distr = ds.normal(10 * u.cm, n_samples=100, std=1 * u.cm)
    todistr = distr.to(u.m)
    assert_quantity_allclose(distr.pdf_mean().to(u.m), todistr.pdf_mean())


def test_distr_noq_to():
    # this is an array distribution not a quantity
    distr = ds.normal(10, n_samples=100, std=1)
    with pytest.raises(AttributeError):
        distr.to(u.m)


def test_distr_to_value():
    distr = ds.normal(10 * u.cm, n_samples=100, std=1 * u.cm)
    tovdistr = distr.to_value(u.m)
    assert np.allclose(distr.pdf_mean().to_value(u.m), tovdistr.pdf_mean())


def test_distr_noq_to_value():
    distr = ds.normal(10, n_samples=100, std=1)
    with pytest.raises(AttributeError):
        distr.to_value(u.m)


def test_distr_angle():
    # Check that Quantity subclasses decay to Quantity appropriately.
    distr = Distribution([2.0, 3.0, 4.0])
    ad = Angle(distr, "deg")
    ad_plus_ad = ad + ad
    assert isinstance(ad_plus_ad, Angle)
    assert isinstance(ad_plus_ad, Distribution)

    ad_times_ad = ad * ad
    assert not isinstance(ad_times_ad, Angle)
    assert isinstance(ad_times_ad, u.Quantity)
    assert isinstance(ad_times_ad, Distribution)

    ad += ad
    assert isinstance(ad, Angle)
    assert isinstance(ad, Distribution)
    assert_array_equal(ad.distribution, ad_plus_ad.distribution)

    with pytest.raises(u.UnitTypeError):
        ad *= ad


def test_distr_angle_view_as_quantity():
    # Check that Quantity subclasses decay to Quantity appropriately.
    distr = Distribution([2.0, 3.0, 4.0])
    ad = Angle(distr, "deg")
    qd = ad.view(u.Quantity)
    assert not isinstance(qd, Angle)
    assert isinstance(qd, u.Quantity)
    assert isinstance(qd, Distribution)
    # View directly as DistributionQuantity class.
    qd2 = ad.view(qd.__class__)
    assert not isinstance(qd2, Angle)
    assert isinstance(qd2, u.Quantity)
    assert isinstance(qd2, Distribution)
    assert_array_equal(qd2.distribution, qd.distribution)
    qd3 = ad.view(qd.dtype, qd.__class__)
    assert not isinstance(qd3, Angle)
    assert isinstance(qd3, u.Quantity)
    assert isinstance(qd3, Distribution)
    assert_array_equal(qd3.distribution, qd.distribution)


def test_distr_cannot_view_new_dtype():
    # A Distribution has a very specific structured dtype with just one
    # element that holds the array of samples.  As it is not clear what
    # to do with a view as a new dtype, we just error on it.
    # TODO: with a lot of thought, this restriction can likely be relaxed.
    distr = Distribution([2.0, 3.0, 4.0])
    with pytest.raises(ValueError, match="with a new dtype"):
        distr.view(np.dtype("f8"))

    # Check subclass just in case.
    ad = Angle(distr, "deg")
    with pytest.raises(ValueError, match="with a new dtype"):
        ad.view(np.dtype("f8"))

    with pytest.raises(ValueError, match="with a new dtype"):
        ad.view(np.dtype("f8"), Distribution)


def test_scalar_quantity_distribution():
    # Regression test for gh-12336
    angles = Distribution([90.0, 30.0, 0.0] * u.deg)
    sin_angles = np.sin(angles)  # This failed in 4.3.
    assert isinstance(sin_angles, Distribution)
    assert isinstance(sin_angles, u.Quantity)
    assert_array_equal(sin_angles, Distribution(np.sin([90.0, 30.0, 0.0] * u.deg)))
