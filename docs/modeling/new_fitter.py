from astropy.modeling.fitting import Fitter

class Chi2Fitter(Fitter):
    """Chi-square (a.k.a. least squares) fitter.
    
    Uses the SLSQP minimiser from scipy.
    """
    supported_constraints = ['bounds', 'eqcons', 'ineqcons', 'fixed', 'tied']

    def errorfunc(self, fitparams, *args):
        model, x, y = args
        self._fitter_to_model_params(model, fitparams)
        return ((y - model(x)) ** 2).sum()

    def __call__(self, model, x, y , maxiter=100):
        from scipy.optimize import fmin_slsqp
        if not model.fittable:
            raise ValueError("Model must be a subclass of ParametricModel")
        self._validate_constraints(model)
        model_copy = model.copy()
        x0, _ = model_copy._model_to_fit_params()
        fitparams = fmin_slsqp(self.errorfunc, x0=x0,
                               args=(model_copy, x, y),
                               iter=maxiter)
        self._fitter_to_model_params(model_copy, fitparams)
        return model_copy
