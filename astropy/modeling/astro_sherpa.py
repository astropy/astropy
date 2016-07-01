import numpy as np
from collections import OrderedDict
from sherpa.fit import Fit
from sherpa.data import Data1D, Data1DInt, Data2D, Data2DInt, DataSimulFit
from sherpa.models import UserModel, Parameter, SimulFitModel
from .fitting import Fitter
from ..utils.exceptions import AstropyUserWarning
from sherpa.stats import Chi2, Chi2ConstVar, Chi2DataVar, Chi2Gehrels, Chi2ModVar, Chi2XspecVar, LeastSq
from sherpa.optmethods import GridSearch, LevMar, MonCar, NelderMead
from sherpa.estmethods import Confidence, Covariance, Projection
# from astropy.modeling

__all__ = ('SherpaFitter',)


class SherpaFitter(Fitter):

    """
    Sherpa Fitter for astropy models. Yay :)

    Parameters
        ----------
        optimizer : string
            the name of a sherpa optimizer.
        statistic : string
            the name of a sherpa statistic.
        estmethod : string
            the name of a sherpa estmethod.
    """

    def __init__(self, optimizer="levmar", statistic="leastsq", estmethod="covariance"):
        try:
            optimizer = optimizer.optmethod
        except AttributeError:
            optimizer = OptMethod(optimizer).optmethod

        try:
            statistic = statistic.stat
        except AttributeError:
            statistic = Stat(statistic).stat

        super(SherpaFitter, self).__init__(optimizer=optimizer, statistic=statistic)

        try:
            self._est_method = estmethod.estmethod
        except AttributeError:
            self._est_method = EstMethod(estmethod).estmethod

        self.fit_info = {}
        self._fitter = None  # a handle for sherpa fit function
        self._fitmodel = None  # a handle for sherpa fit model
        self._data = None  # a handle for sherpa dataset

    def __call__(self, models, x, y, z=None, xerr=None, yerr=None, zerr=None, **kwargs):
        """
        Fit the astropy model with a the sherpa fit routines.

        Parameters
        ----------
        model : `~astropy.modeling.FittableModel`
            model to fit to x, y, z
        x : array
            input coordinates
        y : array
            input coordinates
        z : array (optional)
            input coordinates
        xerr : array (optional)
            an array of errors in x
        yerr : array (optional)
            an array of errors in y
        zerr : array (optional)
            an array of errors in z
        **kwargs:
            keyword arguments will be passed on to sherpa fit routine

        Returns
        -------
        model_copy : `~astropy.modeling.FittableModel` or a list of models.
            a copy of the input model with parameters set by the fitter
        """

        tie_list = []
        try:
            self._data = Dataset(models[0].n_inputs, x, y, z, xerr, yerr, zerr)
        except TypeError:
            self._data = Dataset(models.n_inputs, x, y, z, xerr, yerr, zerr)

        if self._data.ndata > 1:
            if len(models) == 1:
                self._fitmodel = ConvertedModel([models.copy() for _ in xrange(self._data.ndata)], tie_list)
                # Copy the model so each data set has the same model!
            elif len(models) == self._data.ndata:
                self._fitmodel = ConvertedModel(models, tie_list)
            else:
                raise Exception("Don't know how to handle multiple models unless there is one foreach dataset")
        else:
            if len(models) > 1:
                self._data.make_simfit(len(models))
                self._fitmodel = ConvertedModel(models, tie_list)
            else:
                self._fitmodel = ConvertedModel(models)

        self._fitter = Fit(self._data.data, self._fitmodel.sherpa_model, self._stat_method, self._opt_method, self._est_method, **kwargs)
        self.fit_info = self._fitter.fit()

        return self._fitmodel.get_astropy_model()

    def est_errors(self, sigma=None, maxiters=None, numcores=1, methoddict=None, parlist=None):
        """
        Use sherpa error estimators based on the last fit.

        Parameters:
            sigma: float
                this will be set as the confidance interval for which the errors are found too.
            maxiters: int
                the maximum number of iterations the error estimator will run before giving up.
            methoddict: dict
                !not quite sure couldn't figure this one out yet!
            parlist: list
                a list of parameters to find the confidance interval of if none are provided all free
                parameters will be estimated.
        """
        if self._fitter is None:
            ValueError("Must complete a valid fit before errors can be calculated")
        if sigma is not None:
            self._fitter.estmethod.config['sigma'] = sigma
        if maxiters is not None:
            self._fitter.estmethod.config['maxiters'] = maxiters
        if 'numcores' in self._fitter.estmethod.config:
            if not numcores == self._fitter.estmethod.config['numcores']:
                self._fitter.estmethod.config['numcores'] = numcores

        return self._fitter.est_errors(methoddict=methoddict, parlist=parlist)


class SherpaWrapper(object):
    pass

sherpa_stats = {'chi2': Chi2, 'chi2constvar': Chi2ConstVar, 'chi2datavar': Chi2DataVar, 'chi2gehrels': Chi2Gehrels, 'chi2modvar': Chi2ModVar, 'chi2xspecvar': Chi2XspecVar, 'leastsq': LeastSq}


class Stat(SherpaWrapper):
    """
    A wrapper for the fit statistics of sherpa

    Parameter:
        statname: String
            the name of a sherpa statistics.
    """

    stat = None

    def __init__(self, statname=None):
        if statname is not None:
            self.set_stat(statname)

    def set_stat(self, statname):
        """
        Set the sherpa fit statistic

        statname: String
            the name of a sherpa fit statistic.
        """
        if statname.lower() in sherpa_stats:
            self.stat = sherpa_stats[statname.lower()]

    def list_stats(self,):
        """
        Returns a list of the sherpa fit statistics supported.
        """
        return sherpa_stats.keys()


sherpa_optmethods = {'simplex': GridSearch, 'levmar': LevMar, 'moncar': MonCar, 'neldermead': NelderMead}


class OptMethod(SherpaWrapper):

    """
    A wrapper for the optimization methods of sherpa

    Parameter:
        optmethoname: String
            the name of a sherpa optimization method.
    """

    optmethod = None

    def __init__(self, optmethodname=None):
        if optmethodname is not None:
            self.set_optmethod(optmethodname)

    def set_optmethod(self, optmethodname):
        """
        Set the sherpa optimization method

        estmethoname: String
            the name of a sherpa optimization method.
        """

        if optmethodname.lower() in sherpa_optmethods:
            self.optmethod = sherpa_optmethods[optmethodname.lower()]

    def list_optmethods(self,):
        """
        Returns a list of the sherpa optimization methods supported.
        """
        return sherpa_optmethods.keys()

sherpa_estmethods = {'confidence': Confidence, 'covariance': Covariance, 'projection': Projection}


class EstMethod(SherpaWrapper):
    """
    A wrapper for the error estimation methods of sherpa

    Parameter:
        estmethoname: String
            the name of a sherpa error estimation method.
    """

    estmethod = None

    def __init__(self, estmethodname=None):
        if estmethodname is not None:
            self.set_estmethod(estmethodname)

    def set_estmethod(self, estmethodname):
        """
        Set the sherpa error estimation method
        Parameter:
            estmethoname: String
                the name of a sherpa error estimation method.
        """

        if estmethodname.lower() in sherpa_estmethods:
            self.estmethod = sherpa_estmethods[estmethodname.lower()]

    def list_estmethods(self,):
        """
        Returns a list of the sherpa error estimation methods supported.
        """
        return sherpa_estmethods.keys()


class Dataset(SherpaWrapper):
    """
    Parameters
        ----------
        n_dim: int
            Used to veirfy required number of dimentions.
        x : array (or list of arrays)
            input coordinates
        y : array (or list of arrays)
            input coordinates
        z : array (or list of arrays) (optional)
            input coordinates
        xerr : array (or list of arrays) (optional)
            an array of errors in x
        yerr : array (or list of arrays) (optional)
            an array of errors in y
        zerr : array (or list of arrays) (optional)
            an array of errors in z

    returns:
        _data: a sherpa dataset
    """

    def __init__(self, n_dim, x, y, z=None, xerr=None, yerr=None, zerr=None, bkg=None):

        x = np.array(x)
        y = np.array(y)

        if x.ndim == 2 or (x.dtype == np.object or y.dtype == np.object):
            data = []
            if z is None:
                z = len(x) * [None]

            if xerr is None:
                xerr = len(x) * [None]

            if yerr is None:
                yerr = len(y) * [None]

            if zerr is None:
                zerr = len(z) * [None]

            if bkg is None:
                bkg = len(x) * [None]

            for nn, (xx, yy, zz, xxe, yye, zze, bkg) in enumerate(zip(x, y, z, xerr, yerr, zerr, bkg)):
                data.append(self._make_dataset(n_dim, x=xx, y=yy, z=zz, xerr=xxe, yerr=yye, zerr=zze, bkg=bkg, n=nn))
            self.data = DataSimulFit("wrapped_data", data)
            self.ndata = nn + 1
        else:
            self.data = self._make_dataset(n_dim, x=x, y=y, z=z, xerr=xerr, yerr=yerr, zerr=zerr, bkg=bkg)
            self.ndata = 1

    @staticmethod
    def _make_dataset(n_dim, x, y, z=None, xerr=None, yerr=None, zerr=None, bkg=None, n=0):
        """
        Parameters
            ----------
            n_dim: int
                Used to veirfy required number of dimentions.
            x : array
                input coordinates
            y : array
                input coordinates
            z : array (optional)
                input coordinatesbkg
            xerr : array (optional)
                an array of errors in x
            yerr : array (optional)
                an array of errors in y
            zerr : array (optional)
                an array of errors in z
            n  : int
                used in error reporting

        returns:
            _data: a sherpa dataset
        """

        if (z is None and n_dim > 1) or (z is not None and n_dim == 1):
            raise ValueError("Model and data dimentions don't match in dataset %i" % n)

        if z is None:
            assert x.shape == y.shape, "shape of x and y don't match in dataset %i" % n
        else:
            z = np.asarray(z)
            assert x.shape == y.shape == z.shape, "shapes x,y and z don't match in dataset %i" % n

        if yerr is not None:
            yerr = np.array(yerr)
            assert y.shape == yerr.shape, "y's and yerr's shapes do not match in dataset %i" % n

        if xerr is not None:
            xerr = np.array(xerr)
            assert x.shape == xerr.shape, "x's and xerr's shapes do not match in dataset %i" % n

        if z is not None and zerr is not None:
            zerr = np.array(zerr)
            assert z.shape == zerr.shape, "z's and zerr's shapes do not match in dataset %i" % n

        if z is None:
            if xerr is None:
                if yerr is None:
                    data = Data1D("wrapped_data", x=x, y=y)
                else:
                    data = Data1D("wrapped_data", x=x, y=y, staterror=yerr)
            else:
                if yerr is None:
                    data = Data1DInt("wrapped_data", xlo=x - xerr, xhi=x + xerr, y=y)
                else:
                    data = Data1DInt("wrapped_data", xlo=x - xerr, xhi=x + xerr, y=y, staterror=yerr)
        else:
            if xerr is None and yerr is None:
                if zerr is None:
                    data = Data2D("wrapped_data", x0=x, x1=y, y=z)
                else:
                    data = Data2D("wrapped_data", x0=x, x1=y, y=z, staterror=zerr)
            elif xerr is not None and yerr is not None:
                if zerr is None:
                    data = Data2DInt("wrapped_data", x0lo=x - xerr, x0hi=x + xerr, x1lo=y - yerr, x1hi=y + yerr, y=z)
                else:
                    data = Data2DInt("wrapped_data", x0lo=x - xerr, x0hi=x + xerr, x1lo=y - yerr, x1hi=y + yerr, y=z, staterror=zerr)
            else:
                raise ValueError("Set xerr and yerr, or set neither!")

        return data

    def make_simfit(self, numdata):
        """
        This makes a single datasets into a simdatafit at allow fitting of multiple models by copying the single dataset!

        Parameter:

        numdata: int
            the number of times you want to copy the dataset i.e if you want 2 datasets total you put 1!
        """

        self.data = DataSimulFit("wrapped_data", [self.data for _ in xrange(numdata)])
        self.ndata = numdata + 1


class ConvertedModel(SherpaWrapper):
    """
    This  wraps the model convertion to sherpa models and from astropy models and back!

    Parameters:
        models: model : `~astropy.modeling.FittableModel` (or list of)

        tie_list: list (optional)
            a list of parameter pairs which will be tied accross models
            e.g. [(modelB.y, modelA.x)] will mean that y in modelB will be tied to x of modelA
    """

    def __init__(self, models, tie_list=None):
        self.model_dict = OrderedDict()
        try:
            models.parameters  # does it quack
            self.sherpa_model = self._astropy_to_sherpa_model(models)
            self.model_dict[models] = self.sherpa_model
        except AttributeError:
            for mod in models:
                self.model_dict[mod] = self._astropy_to_sherpa_model(mod)
                if tie_list is not None:
                    for par1, par2 in tie_list:
                        getattr(self.model_dict[par1._model], par1.name).link = getattr(self.model_dict[par2._model], par2.name)
            self.sherpa_model = SimulFitModel("wrapped_fit_model", self.model_dict.values())

    @staticmethod
    def _astropy_to_sherpa_model(model):
        """
        Converts the model using sherpa's usermodel suppling the parameter detail to sherpa
        then using a decorator to allow the call method to act like the calc method
        """
        def _calc2call(func):
            """This decorator makes call and calc work together."""
            def _converter(inp, *x):
                if func.n_inputs == 1:
                    retvals = func.evaluate(x[0], *inp)
                else:
                    retvals = func.evaluate(x[0], x[1], *inp)
                return retvals
            return _converter

        if len(model.ineqcons) > 0 or len(model.eqcons) > 0:
            AstropyUserWarning('In/eqcons are not supported by sherpa these will be ignored!')

        pars = []
        linkedpars = []
        for pname in model.param_names:
            param = getattr(model, pname)
            vals = [param.name, param.value, param.min, param.max, param.min,
                    param.max, None, param.fixed, False]
            attrnames = ["name", "val", "min", "max", "hard_min", "hard_max",
                         "units", "frozen", "alwaysfrozen"]
            if model.name is None:
                model._name = ""

            pars.append(Parameter(modelname="wrap_" + model.name, **dict([(atr, val) for atr, val in zip(attrnames, vals) if val is not None])))
            if param.tied is not False:
                linkedpars.append(pname)

        smodel = UserModel(model.name, pars)
        smodel.calc = _calc2call(model)

        for pname in linkedpars:
            param = getattr(model, pname)
            sparam = getattr(smodel, pname)
            sparam.link = param.tied(smodel)

        return smodel

    def get_astropy_model(self):
        """Returns an astropy model based on the sherpa model"""
        return_models = []

        for apymod, shmod in self.model_dict.items():
            return_models.append(apymod.copy())
            for pname, pval in map(lambda p: (p.name, p.val), shmod.pars):
                getattr(return_models[-1], pname.split(".")[-1]).value = pval

        if len(return_models) > 1:
            return return_models
        else:
            return return_models[0]
