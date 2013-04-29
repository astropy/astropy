"""
An example of how to define a user-defined new fitter
Note that such a chi^2 fitter is available as
astropy.models.fitting.SLSQPFitter, the point here is
to illustrate how you could write your own if you e.g.
want to use a different optimizer or different fit statistic.
"""
import numpy as np
from scipy.optimize import fmin_slsqp
from astropy.models.fitting import Fitter

from astropy.models.fitting import constraintsdef
constraintsdef['Chi2Fitter'] = []

class Chi2Fitter(Fitter):
    """Performs a chi^2 fit using scipy.optimize.fmin_slsqp"""
    def __init__(self, model):
        Fitter.__init__(self, model)
                        
    def errorfunc(self, fps, *args):
        """Here you defined the fit statistic (chi ^ 2 in this example)"""
        meas = args[0]
        self.fitpars = fps
        res = self.model(*args[1:]) - meas
        return np.sum(res ** 2)
    
    def __call__(self, x, y , maxiter=100, epsilon=10 ** (-12)):
        """Here you run the fit, i.e. call the optimizer"""
        bounds = [self.model.constraints.bounds[key] for key in self.model.parnames]
        self.fitpars = fmin_slsqp(self.errorfunc, x0=self.model.parameters[:],
                                  args=(y, x), bounds=bounds)

if __name__ == '__main__':
    # Example fit: 1D Gauss
    from astropy.models.builtin_models import Gauss1DModel
    from astropy.models.fitting import SLSQPFitter
    
    # Set up model and data
    gauss_model = Gauss1DModel(10., xsigma=2.1, xcen=4.2)
    np.random.seed(0)
    x = np.arange(1, 10, .1)
    y = gauss_model(x) + np.random.randn(len(x))

    # Run fit with user-defined fitter
    gauss_chi2_fit = Chi2Fitter(gauss_model)
    gauss_chi2_fit(x, y)  # This call of the fitter actually performs the fit
    print(gauss_model)

    # Run fit with build-in chi2 fitter ... gives identical results
    gauss_chi2_fit = SLSQPFitter(gauss_model)
    gauss_chi2_fit(x, y)  # This call of the fitter actually performs the fit
    print(gauss_model)
