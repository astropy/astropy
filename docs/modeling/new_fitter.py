__doctest_requires__ = {('Chi2Fitter.*'): ['scipy']}
from astropy.modeling.fitting import Fitter

class Chi2Fitter(Fitter):
    """Chi-square (a.k.a. least squares) fitter.
    
    Wraps the `scipy.optimize.minimize` function.
    """
    supported_constraints = ['bounds', 'fixed']

    def errorfunc(self, fitparams, *args):
        """Compute chi-square fit statistic."""
        model, x, y = args
        self._fitter_to_model_params(model, fitparams)
        return ((y - model(x)) ** 2).sum()

    def __call__(self, model, x, y , maxiter=100):
        """Minimize objective function."""
        from scipy.optimize import minimize
        self._validate_constraints(model)
        model_copy = model.copy()
        x0, _ = model_copy._model_to_fit_params()
        fitparams = minimize(self.errorfunc, x0=x0,
                            args=(model_copy, x, y),
                            iter=maxiter)
        self._fitter_to_model_params(model_copy, fitparams)
        return model_copy
