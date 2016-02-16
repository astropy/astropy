import numpy as np
from ... import units as u
from ..distribution import Distribution
from ...utils import NumpyRNGContext
from ...tests.helper import assert_quantity_allclose

def test_numpy_init():
    # Test that we can initialize directly from a Numpy array if we provide a unit
    parr = np.random.poisson([1, 5, 30, 400], (1000, 4))
    distr = Distribution(parr, u.kpc)
    
    
def test_quantity_init():
    # Test that we can initialize directly from a Quantity
    pq = np.random.poisson([1, 5, 30, 400], (1000, 4)) * u.kpc
    distr = Distribution(pq)


class TestDistributionStatistics():
    
    def setup_class(self):
    
        with NumpyRNGContext(12345):
            self.data = np.random.normal([1, 2, 3, 4], [3, 2, 4, 5], (10000, 4))
    
        self.distr = Distribution(self.data  * u.kpc)

    def test_shape(self):
        # Distribution shape
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
        assert_quantity_allclose(self.distr.pdf_mean, np.mean(self.data, axis=0) * u.kpc)
        assert_quantity_allclose(self.distr.pdf_mean, [1, 2, 3, 4] * u.kpc, rtol=0.01)

    def test_pdf_std(self):
        # Standard deviation of each PDF
        assert_quantity_allclose(self.distr.pdf_std, np.std(self.data, axis=0) * u.kpc)
        assert_quantity_allclose(self.distr.pdf_std, [3, 2, 4, 5] * u.kpc, rtol=0.01)

    def test_pdf_var(self):
        # Variance of each PDF
        assert_quantity_allclose(self.distr.pdf_var, np.var(self.data, axis=0) * u.kpc)
        assert_quantity_allclose(self.distr.pdf_var, [9, 4, 16, 25] * u.kpc, rtol=0.01)

    def test_pdf_median(self):
        # Median of each PDF
        assert_quantity_allclose(self.distr.pdf_median, np.median(self.data, axis=0) * u.kpc)
        assert_quantity_allclose(self.distr.pdf_median, [1, 2, 3, 4] * u.kpc, rtol=0.01)

    def test_pdf_mad(self):
        # Median absolute deviation of each PDF
        mad = np.median(np.abs(self.data - np.median(self.data, axis=0)), axis=0)
        assert_quantity_allclose(self.distr.pdf_median, mad * u.kpc)
        # FIXME: need to use a more accurate scalefactor
        assert_quantity_allclose(self.distr.pdf_median, [1, 2, 3, 4] * u.kpc / 1.4826, rtol=0.01)

    def test_pdf_smad(self):
        # Median absolute deviation of each PDF, rescaled to match std for normal
        # FIXME: need to use a more accurate scalefactor
        mad = np.median(np.abs(self.data - np.median(self.data, axis=0)), axis=0) * 1.4826
        assert_quantity_allclose(self.distr.pdf_median, mad * u.kpc)
        assert_quantity_allclose(self.distr.pdf_median, [1, 2, 3, 4] * u.kpc, rtol=0.01)





    
        


