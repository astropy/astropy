import inspect
import numpy as np
import operator as op
from collections import OrderedDict
from pyparsing import nestedExpr
from sherpa.stats import Stat, LeastSq
from sherpa.fit import Fit
from sherpa.optmethods import OptMethod, LevMar
from sherpa.estmethods import EstMethod, Covariance
from sherpa.data import Data1D, Data1DInt, Data2D, Data2DInt, DataSimulFit
from sherpa.models import UserModel, Parameter, SimulFitModel
from .fitting import Fitter
from .core import FittableModel, _CompoundModel
# from astropy.modeling

__all__ = ('SherpaFitter',)

# need to check aviable operators
OPERATORS = OrderedDict([
    ('^', op.pow),
    ('**', op.pow),
    ('*', op.mul),
    ('/', op.div),
    ('//', op.floordiv),
    ('+', op.add),
    ('-', op.sub)
])


def _expression_to_list(expr):
    """
    This takes an expression i.e '[0]+[1]*([2]+[3])' and returns the list ['[0]+[1]*',['[2]+[3]']]

    Parameters
    ----------
        expr: parseString
            this expression you wish to pars
        returns: list
            a list of parsed expression
    """
    return nestedExpr().parseString("(" + expr + ")").asList()[0]


def _parse_expression(exp_list, models):
    """
    This takes an expression list and a list of models it coverts the models
    then constucts a compound model based on the expression list.

    Parameters:
        exp_list: list
            a list expression like the output of _expression_to_list
        models: list
            a list of astropy models which are to be coverted
    """
    to_return = []
    if isinstance(exp_list, list):
        for part in exp_list:
            if isinstance(part, list):
                to_return.extend(_parse_expression(part, models))
            else:
                if part in OPERATORS.keys():
                    to_return.append(part)
                else:
                    to_return.append(models[part.replace("[", "").replace("]", "")])

        for txtop, oper in OPERATORS.items():
            while txtop in to_return:
                opindex = to_return.index(txtop)
                to_return[opindex - 1] = oper(to_return[opindex - 1], to_return[opindex + 1])
                to_return.pop(opindex + 1)
                to_return.pop(opindex)
        return to_return
    else:
        return exp_list


def _astropy_to_sherpa_model(model, deconstruct=False):
    """
    This takes astropy models and converts them to sherpa models.
    If deconstruct is true then the submodels are converted then recombined.

    Parameters
    ----------
        model : `~astropy.modeling.FittableModel`
            the model you wish to convert
        diconstruct: bool
            if True each submodel is converted then re-assembled based on the
            model._format_expression()
    """

    def _convert_to_sherpa_model(model):
        '''
        Converts the model using sherpa's usermodel suppling the parameter detail to sherpa
        then using a decorator to allow the call method to act like the calc method
        '''
        def _calc2call(func):
            '''This decorator makes call and calc work together.'''
            def _converter(inp, *x):
                if func.n_inputs == 1:
                    # startp = func._parameters[:]
                    # func._parameters = np.array(inp)
                    # retvals = func.__call__(*x)
                    retvals = func.evaluate(x[0], *inp)
                    # func._parameters = startp
                else:
                    retvals = func.evaluate(x[0], x[1], *inp)  # I need to test the 2d case!
                return retvals
            return _converter

        pars = []
        linkedpars = []
        for pname in model.param_names:
            param = model.__getattribute__(pname)
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
            param = model.__getattribute__(pname)
            sparam = smodel.__getattribute__(pname)
            sparam.link = param.tied(smodel)

        return smodel

    if deconstruct:
        models = [_convert_to_sherpa_model(m) for m in model.submodels]
        expression = model._format_expression()

        if isinstance(models, list) or isinstance(models, tuple):
            models = {str(num): mod for num, mod in enumerate(models)}

            expression_list = _expression_to_list(expression)
            return [_parse_expression(elist, models)[0] for elist in expression_list]
        else:
            expression_list = _expression_to_list(expression)
            return _parse_expression(expression_list, models)[0]
    else:
        return _convert_to_sherpa_model(model)


def _make_dataset(n_dim, x, y, z=None, xerr=None, yerr=None, zerr=None, n=0):
    '''
    Parameters
        ----------
        n_dim: int
            Used to veirfy required number of dimentions.
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
        n  : int
            used in error reporting

    returns:
        _data: a sherpa dataset
    '''

    if (z is None and n_dim > 1) or (z is not None and n_dim == 1):
        raise Exception("Model and data dimentions don't match in dataset %i" % n)

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

    '''
    So the generall method is to dstack(and transpose each row) so
        [x1,x2,x3],[y1,y2,y3],[z1,z2,z3],...
    becomes
        [x1,y1,,z1],[x2,y2,,z2],[x3,y3,,z3],...
        the and put it into a data set! this assumes that all columns are set for all datasets!
        athough a none can be put in place of missing columns
    '''
    if z is None:
        if xerr is None:
            if yerr is None:
                _data = Data1D("wrapped_data", x=x, y=y)
            else:
                _data = Data1D("wrapped_data", x=x, y=y, staterror=yerr)
        else:
            if yerr is None:
                _data = Data1DInt("wrapped_data", xlo=x - xerr, xhi=x + xerr, y=y)
            else:
                _data = Data1DInt("wrapped_data", xlo=x - xerr, xhi=x + xerr, y=y, staterror=yerr)
    else:
        if xerr is None and yerr is None:
            if zerr is None:
                _data = Data2D("wrapped_data", x0=x, x1=y, y=z)
            else:
                _data = Data2D("wrapped_data", x0=x, x1=y, y=z, staterror=zerr)
        elif xerr is not None and yerr is not None:
            if zerr is None:
                _data = Data2DInt("wrapped_data", x0lo=x - xerr, x0hi=x + xerr, x1lo=y - yerr, x1hi=y + yerr, y=z)
            else:
                _data = Data2DInt("wrapped_data", x0lo=x - xerr, x0hi=x + xerr, x1lo=y - yerr, x1hi=y + yerr, y=z, staterror=zerr)
        else:
            raise Exception("Set xerr and yerr, or set neither!")

    return _data


def make_datasets(n_dim, x, y, z=None, xerr=None, yerr=None, zerr=None):
    '''
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
    '''

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

        for nn, (xx, yy, zz, xxe, yye, zze) in enumerate(zip(x, y, z, xerr, yerr, zerr)):
            data.append(_make_dataset(n_dim, x=xx, y=yy, z=zz, xerr=xxe, yerr=yye, zerr=zze, n=nn))
        return DataSimulFit("wrapped_data", data)
    else:
        return _make_dataset(n_dim, x=x, y=y, z=z, xerr=xerr, yerr=yerr, zerr=zerr)


def _make_simfit_model(models, tie_list, model_dict):
    '''
    makes a SimulFitModel for the models suplied it also links and parameters within tie_list
    model_dict updated with astropy models and their converted conterparts this is used for
    repopulating the model later
    '''
    for mod in models:
        model_dict[mod] = _astropy_to_sherpa_model(mod)

    if tie_list is not None:
        for par1, par2 in tie_list:
            getattr(model_dict[par1._model], par1.name).link = getattr(model_dict[par2._model], par2.name)
    return SimulFitModel("wrapped_fit_model", model_dict.values())


class SherpaFitter(Fitter):

    """
    Sherpa Fitter for astropy models. Yay :)

    Parameters
        ----------
        optimizer : (callable)
                    A sherpa optimizer method.
        statistic : (callable)
                    A sherpa statistic function.
        estmethod : (callable)
                    A shepra error estimation method.
    """

    def __init__(self, optimizer=LevMar, statistic=LeastSq, estmethod=Covariance):
        '''
        check the types of optimizer, statistic and estmethod are sherpa based classes
         so we can use them
        '''
        if inspect.isclass(optimizer) and not issubclass(optimizer, OptMethod):
            raise TypeError("optimizer must be sherpa OptMethod class or isinstance")
        elif inspect.isfunction(optimizer) and not issubclass(optimizer.__class__, OptMethod):
            raise TypeError("optimizer must be sherpa OptMethod class or isinstance")

        if inspect.isclass(statistic) and not issubclass(statistic, Stat):
            raise TypeError("statistic must be sherpa Stat class or isinstance")
        elif inspect.isfunction(statistic) and not issubclass(statistic.__class__, Stat):
            raise TypeError("statistic must be sherpa Stat class or isinstance")

        if inspect.isclass(estmethod) and not issubclass(estmethod, EstMethod):
            raise TypeError("estmethod must be sherpa EstMethod class or isinstance")
        elif inspect.isfunction(estmethod) and not issubclass(estmethod.__class__, EstMethod):
            raise TypeError("estmethod must be sherpa EstMethod class or isinstance")

        super(SherpaFitter, self).__init__(optimizer=optimizer, statistic=statistic)
        if inspect.isclass(estmethod):
            self._est_method = estmethod()
        else:
            self._est_method = estmethod

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

        model_dict = OrderedDict()
        tie_list = []

        if isinstance(models, (FittableModel, _CompoundModel)):
            models = [models]

        self._data = make_datasets(models[0].n_inputs, x, y, z, xerr, yerr, zerr)

        if isinstance(self._data, DataSimulFit):
            _ndata = len(self._data.datasets)
            if len(models) == 1:
                self._fitmodel = _make_simfit_model([models[0].copy() for _ in xrange(_ndata)], tie_list, model_dict)
                # Copy the model so each data set has the same model!
            elif len(models) == _ndata:
                self._fitmodel = _make_simfit_model(models, tie_list, model_dict)
            else:
                raise Exception("Don't know how to handle multiple models unless there is one foreach dataset")
        else:
            if len(models) > 1:
                self._data = DataSimulFit("wrapped_data", [self._data for _ in xrange(len(models))])
                self._fitmodel = _make_simfit_model(models, tie_list, model_dict)
            else:
                self._fitmodel = _astropy_to_sherpa_model(models[0])

        self._fitter = Fit(self._data, self._fitmodel, self._stat_method, self._opt_method, self._est_method, **kwargs)
        self.fit_info = self._fitter.fit()

        if isinstance(self._fitmodel, SimulFitModel):
            return_models = []
            for apymod, shmod in model_dict.items():
                return_models.append(apymod.copy())
                for pname, pval in map(lambda p: (p.name, p.val), shmod.pars):
                    getattr(return_models[-1], pname.split(".")[-1]).value = pval
        else:
            return_models = models[0].copy()
            for pname, pval in map(lambda p: (p.name, p.val), self._fitmodel.pars):
                getattr(return_models, pname.split(".")[-1]).value = pval
        return return_models

    def est_errors(self, sigma=None, maxiters=None, numcores=1, methoddict=None, parlist=None):
        '''
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
        '''
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
