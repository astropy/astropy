# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import absolute_import, division, print_function

import numpy as np

from scipy.stats import norm

from ... import units as u
from ...utils import NumpyRNGContext
from ...tests.helper import assert_quantity_allclose, pytest


SMAD_FACTOR = 1 / norm.ppf(0.75)


def test_numpy_init():
    # Test that we can initialize directly from a Numpy array if we provide a unit
    parr = np.random.poisson([1, 5, 30, 400], (1000, 4))
    u.Distribution(parr, u.kpc)


def test_quantity_init():
    # Test that we can initialize directly from a Quantity
    pq = np.random.poisson([1, 5, 30, 400], (1000, 4)) * u.kpc
    u.Distribution(pq)


def test_numpy_init_nounit():
    parr = np.random.poisson([1, 5, 30, 400], (1000, 4))
    with pytest.raises(u.UnitsError) as exc:
        u.Distribution(parr)
    # FIXME: update error message as needed
    assert exc.value.args[0] == "No unit specified when creating distribution"


class TestDistributionStatistics():

    def setup_class(self):

        with NumpyRNGContext(12345):
            self.data = np.random.normal([1, 2, 3, 4], [3, 2, 4, 5], (10000, 4))

        self.distr = u.Distribution(self.data  * u.kpc)

    def test_shape(self):
        # u.Distribution shape
        assert self.distr.shape == (10000, 4)

    def test_size(self):
        # Total number of values
        assert self.distr.size == 40000

    def test_n_samples(self):
        # Number of samples
        assert self.distr.n_samples == 10000

    def test_n_distr(self):
        # Shape of the PDF (note, this is actually the number of values regardless of samples, needs a better name?)
        assert self.distr.n_distr == (4,)

    def test_pdf_mean(self):
        # Mean of each PDF
        expected = np.mean(self.data, axis=0) * u.kpc
        assert_quantity_allclose(self.distr.pdf_mean, expected)
        assert_quantity_allclose(self.distr.pdf_mean, [1, 2, 3, 4] * u.kpc, rtol=0.01)

    def test_pdf_std(self):
        # Standard deviation of each PDF
        expected = np.std(self.data, axis=0) * u.kpc
        assert_quantity_allclose(self.distr.pdf_std, expected)
        assert_quantity_allclose(self.distr.pdf_std, [3, 2, 4, 5] * u.kpc, rtol=0.01)

    def test_pdf_var(self):
        # Variance of each PDF
        expected = np.var(self.data, axis=0) * u.kpc
        assert_quantity_allclose(self.distr.pdf_var, expected)
        assert_quantity_allclose(self.distr.pdf_var, [9, 4, 16, 25] * u.kpc, rtol=0.01)

    def test_pdf_median(self):
        # Median of each PDF
        expected = np.median(self.data, axis=0) * u.kpc
        assert_quantity_allclose(self.distr.pdf_median, expected)
        assert_quantity_allclose(self.distr.pdf_median, [1, 2, 3, 4] * u.kpc, rtol=0.01)

    def test_pdf_mad(self):
        # Median absolute deviation of each PDF
        median = np.median(self.data, axis=0)
        expected = np.median(np.abs(self.data - median), axis=0) * u.kpc
        assert_quantity_allclose(self.distr.pdf_mad, expected)
        assert_quantity_allclose(self.distr.pdf_smad, self.distr.pdf_mad * SMAD_FACTOR)
        assert_quantity_allclose(self.distr.pdf_smad, [1, 2, 3, 4] * u.kpc, rtol=0.01)

    def test_percentile(self):
        # TODO: numpy defines percentiles from 0 to 100, we should do the same
        # for consistency.
        expected = np.percentile(self.data, [10, 50, 90], axis=0) * u.kpc
        assert_quantity_allclose(self.distr.percentiles([10, 50, 90]), expected)

    def test_add_quantity(self):
        distrplus = self.distr + [2000, 0, 0, 500] * u.pc
        expected = (np.median(self.data, axis=0) + np.array([2, 0, 0, 0.5])) * u.kpc
        assert_quantity_allclose(distrplus.pdf_median, expected)
        expected = np.var(self.data, axis=0) * u.kpc
        assert_quantity_allclose(distrplus.pdf_var, expected)

    def test_add_distribution(self):
        another_data = (np.random.randn(10000, 4)
                                      * [1000,.01 , 3000, 10]
                                      + [2000, 0, 0, 500])
        another_distr = u.Distribution(another_data, u.pc)
        combined_distr = self.distr + another_distr
        expected = np.median(self.data + another_data, axis=0)
        assert_quantity_allclose(combined_distr.pdf_median, expected)
        expected = np.var(self.data + another_data, axis=0)
        assert_quantity_allclose(combined_distr.pdf_var, expected)

def test_helper_normal(self):
    centerq = [1, 5, 30, 400] * u.kpc
    u.NormalDistribution(centerq, std=[0.2, 1.5, 4, 1]*u.kpc)
    u.NormalDistribution(centerq, var=[0.04, 2.25, 16, 1]*u.kpc**2)
    u.NormalDistribution(centerq, ivar=[25, 0.44444444, 0.625, 1]*u.kpc**-2)

def test_helper_poisson(self):

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

def test_helper_uniform(self):

    # TODO: I think we should define it using low/high instead of center/width,
    # for consistency with Numpy, and because if we want to later do a
    # LogUniformDistribution, it will be easier to express as a range.
    u.UniformDistribution([1, 3, 2] * u.pc, [5, 4, 3] * u.pc)


def test_helper_samples(self):

    centerq = [1, 5, 30, 400] * u.one

    n_dist = u.NormalDistribution(centerq, std=[0.2, 1.5, 4, 1]*u.kpc, n_samples=100)
    assert n_dist.shape == (100, 4)

    n_dist = u.NormalDistribution(centerq, std=[0.2, 1.5, 4, 1]*u.kpc, n_samples=20000)
    assert n_dist.shape == (20000, 4)


def test_arithmetic(self):
    dist = (u.NormalDistribution(3 * u.kpc)
            * u.PoissonDistribution(5 * u.one)
            + u.UniformDistribution(3 * u.pc, 5 * u.pc))


# TODO: add a test to check default number of samples and how to change it
