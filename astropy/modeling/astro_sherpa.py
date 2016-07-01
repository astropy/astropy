import numpy as np
from collections import OrderedDict
from sherpa.fit import Fit
from sherpa.data import Data1D, Data1DInt, Data2D, Data2DInt, DataSimulFit, BaseData
from sherpa.models import UserModel, Parameter, SimulFitModel
from .fitting import Fitter
from ..utils.exceptions import AstropyUserWarning
from sherpa.stats import Chi2, Chi2ConstVar, Chi2DataVar, Chi2Gehrels, Chi2ModVar, Chi2XspecVar, LeastSq
from sherpa.stats import CStat, WStat, Cash
from sherpa.optmethods import GridSearch, LevMar, MonCar, NelderMead
from sherpa.estmethods import Confidence, Covariance, Projection
# from astropy.modeling

__all__ = ('SherpaFitter',)


class SherpaWrapper(object):
    value = None

    def __init__(self, value=None):
        if value is not None:
            self.set(value)

    def set(self, value):
        try:
            self.value = self._sherpa_values[value.lower()]
        except KeyError:
            UserWarning("Value not found")  # todo handle


class Stat(SherpaWrapper):

    """
    A wrapper for the fit statistics of sherpa

    Parameter:
        value: String
            the name of a sherpa statistics.
    """

    _sherpa_values = {'cash': Cash, 'wstat': WStat, 'cstat': CStat,
                      'chi2': Chi2, 'chi2constvar': Chi2ConstVar,
                      'chi2datavar': Chi2DataVar,
                      'chi2gehrels': Chi2Gehrels,
                      'chi2modvar': Chi2ModVar,
                      'chi2xspecvar': Chi2XspecVar,
                      'leastsq': LeastSq}


class OptMethod(SherpaWrapper):

    """
    A wrapper for the optimization methods of sherpa

    Parameter:
        value: String
            the name of a sherpa optimization method.
    """
    _sherpa_values = {'simplex': GridSearch, 'levmar': LevMar,
                      'moncar': MonCar, 'neldermead': NelderMead}


class EstMethod(SherpaWrapper):

    """
    A wrapper for the error estimation methods of sherpa

    Parameter:
        value: String
            the name of a sherpa statistics.
    """

    _sherpa_values = {'confidence': Confidence, 'covariance': Covariance,
                      'projection': Projection}


class SherpaFitter(Fitter):
    __doc__ = """
    Sherpa Fitter for astropy models. Yay :)

    Parameters
        ----------
        optimizer : string
            the name of a sherpa optimizer.
            posible options include:
                {opt}
        statistic : string
            the name of a sherpa statistic.
            posible options include:
                {stat}
        estmethod : string
            the name of a sherpa estmethod.
            posible options include:
                {est}
    """.format(opt=", ".join(OptMethod._sherpa_values.keys()),
               stat=", ".join(Stat._sherpa_values.keys()),
               est=", ".join(EstMethod._sherpa_values.keys()))  # is this evil?

    def __init__(self, optimizer="levmar", statistic="leastsq", estmethod="covariance"):
        try:
            optimizer = optimizer.value
        except AttributeError:
            optimizer = OptMethod(optimizer).value

        try:
            statistic = statistic.value
        except AttributeError:
            statistic = Stat(statistic).value

        super(SherpaFitter, self).__init__(optimizer=optimizer, statistic=statistic)

        try:
            self._est_method = estmethod.value
        except AttributeError:
            self._est_method = EstMethod(estmethod).value

        self.fit_info = {}
        self._fitter = None  # a handle for sherpa fit function
        self._fitmodel = None  # a handle for sherpa fit model
        self._data = None  # a handle for sherpa dataset

    def __call__(self, models, x, y, z=None, xerr=None, yerr=None, zerr=None, bkg=None, bkg_scale=1, **kwargs):
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
           self._data = Dataset(models[0].n_inputs, x, y, z, xerr, yerr, zerr, bkg, bkg_scale)
        except TypeError:
            self._data = Dataset(models.n_inputs, x, y, z, xerr, yerr, zerr, bkg, bkg_scale)

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

    def __init__(self, n_dim, x, y, z=None, xerr=None, yerr=None, zerr=None, bkg=None, bkg_scale=1):

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
            try:
                iter(bkg_scale)
            except TypeError:
                bkg_scale = len(x) * [bkg_scale]


            for nn, (xx, yy, zz, xxe, yye, zze, bkg, bkg_scale) in enumerate(zip(x, y, z, xerr, yerr, zerr, bkg, bkg_scale)):
                data.append(self._make_dataset(n_dim, x=xx, y=yy, z=zz, xerr=xxe, yerr=yye, zerr=zze, bkg=bkg, bkg_scale=bkg_scale, n=nn))
            self.data = DataSimulFit("wrapped_data", data)
            self.ndata = nn + 1
        else:
            self.data = self._make_dataset(n_dim, x=x, y=y, z=z, xerr=xerr, yerr=yerr, zerr=zerr, bkg=bkg, bkg_scale=bkg_scale)
            self.ndata = 1

    @staticmethod
    def _make_dataset(n_dim, x, y, z=None, xerr=None, yerr=None, zerr=None, bkg=None, bkg_scale=1, n=0):
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
                    if bkg is None:
                        data = Data1D("wrapped_data", x=x, y=y)
                    else:
                        data = Data1DBkg("wrapped_data", x=x, y=y, bkg=bkg, bkg_scale=bkg_scale)
                else:
                    if bkg is None:
                        data = Data1D("wrapped_data", x=x, y=y, staterror=yerr)
                    else:
                        data = Data1DBkg("wrapped_data", x=x, y=y, staterror=yerr, bkg=bkg, bkg_scale=bkg_scale)
            else:
                if yerr is None:
                    if bkg is None:
                        data = Data1DInt("wrapped_data", xlo=x - xerr, xhi=x + xerr, y=y)
                    else:
                        data = Data1DIntBkg("wrapped_data", xlo=x - xerr, xhi=x + xerr, y=y, bkg=bkg, bkg_scale=bkg_scale)
                else:
                    if bkg is None:
                        data = Data1DInt("wrapped_data", xlo=x - xerr, xhi=x + xerr, y=y, staterror=yerr)
                    else:
                        data = Data1DIntBkg("wrapped_data", xlo=x - xerr, xhi=x + xerr, y=y, staterror=yerr, bkg=bkg, bkg_scale=bkg_scale)
        else:
            if xerr is None and yerr is None:
                if zerr is None:
                    if bkg is None:
                        data = Data2D("wrapped_data", x0=x, x1=y, y=z)
                    else:
                        data = Data2DBkg("wrapped_data", x0=x, x1=y, y=z, bkg=bkg, bkg_scale=bkg_scale)
                else:
                    if bkg is None:
                        data = Data2D("wrapped_data", x0=x, x1=y, y=z, staterror=zerr)
                    else:
                        data = Data2DBkg("wrapped_data", x0=x, x1=y, y=z, staterror=zerr, bkg=bkg, bkg_scale=bkg_scale)
            elif xerr is not None and yerr is not None:
                if zerr is None:
                    if bkg is None:
                        data = Data2DInt("wrapped_data", x0lo=x - xerr, x0hi=x + xerr, x1lo=y - yerr, x1hi=y + yerr, y=z)
                    else:
                        data = Data2DIntBkg("wrapped_data", x0lo=x - xerr, x0hi=x + xerr, x1lo=y - yerr, x1hi=y + yerr, y=z, bkg=bkg, bkg_scale=bkg_scale)
                else:
                    if bkg is None:
                        data = Data2DInt("wrapped_data", x0lo=x - xerr, x0hi=x + xerr, x1lo=y - yerr, x1hi=y + yerr, y=z, staterror=zerr)
                    else:
                        data = Data2DIntBkg("wrapped_data", x0lo=x - xerr, x0hi=x + xerr, x1lo=y - yerr, x1hi=y + yerr, y=z, staterror=zerr, bkg=bkg, bkg_scale=bkg_scale)
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


class ConvertedModel(object):

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


class Data1DIntBkg(Data1DInt):

    _response_ids = [0]
    _background_ids = [0]

    @property
    def response_ids(self):
        return self._response_ids

    @property
    def background_ids(self):
        return self._background_ids

    @property
    def backscal(self):
        return self._bkg_scale

    def get_background(self, index):
        return self._backgrounds[index]

    def __init__(self, name, xlo, xhi, y, bkg, staterror=None, bkg_scale=1):
        self._bkg = np.asanyarray(bkg)
        self._bkg_scale = bkg_scale
        self.exposure = 1

        self.subtracted = False

        self._backgrounds = [BkgDataset(bkg, bkg_scale)]
        BaseData.__init__(self)

        self.xlo = xlo
        self.xhi = xhi
        self.y = y
        self.staterror = staterror


class Data1DBkg(Data1D):

    _response_ids = [0]
    _background_ids = [0]

    @property
    def response_ids(self):
        return self._response_ids

    @property
    def background_ids(self):
        return self._background_ids

    @property
    def backscal(self):
        return self._bkg_scale

    def get_background(self, index):
        return self._backgrounds[index]

    def __init__(self, name, x, y, bkg, staterror=None, bkg_scale=1):
        self._bkg = np.asanyarray(bkg)
        self._bkg_scale = bkg_scale
        self.exposure = 1
        self.subtracted = False

        self._backgrounds = [BkgDataset(bkg, bkg_scale)]
        BaseData.__init__(self)

        self.x = x
        self.y = y
        self.staterror = staterror


class Data2DIntBkg(Data2DInt):

    _response_ids = [0]
    _background_ids = [0]

    @property
    def response_ids(self):
        return self._response_ids

    @property
    def background_ids(self):
        return self._background_ids

    @property
    def backscal(self):
        return self._bkg_scale

    def get_background(self, index):
        return self._backgrounds[index]

    def __init__(self, name, xlo, xhi, ylo, yhi, z, bkg, staterror=None, bkg_scale=1):
        self._bkg = np.asanyarray(bkg)
        self._bkg_scale = bkg_scale
        self.exposure = 1

        self.subtracted = False

        self._backgrounds = [BkgDataset(bkg, bkg_scale)]
        BaseData.__init__(self)

        self.xlo = xlo
        self.xhi = xhi
        self.ylo = ylo
        self.yhi = yhi
        self.z = z
        self.staterror = staterror


class Data2DBkg(Data2D):

    _response_ids = [0]
    _background_ids = [0]

    @property
    def response_ids(self):
        return self._response_ids

    @property
    def background_ids(self):
        return self._background_ids

    @property
    def backscal(self):
        return self._bkg_scale

    def get_background(self, index):
        return self._backgrounds[index]

    def __init__(self, name, x, y, z, bkg, staterror=None, bkg_scale=1):
        self._bkg = np.asanyarray(bkg)
        self._bkg_scale = bkg_scale
        self.exposure = 1
        self.subtracted = False

        self._backgrounds = [BkgDataset(bkg, bkg_scale)]
        BaseData.__init__(self)

        self.x = x
        self.y = y
        self.z = z
        self.staterror = staterror


class BkgDataset(object):

    def __init__(self, bkg, bkg_scale):
        self._bkg = bkg
        self._bkg_scale = bkg_scale
        self.exposure = 1

    def get_dep(self, flag):
        return self._bkg

    @property
    def backscal(self):
        return self._bkg_scale
