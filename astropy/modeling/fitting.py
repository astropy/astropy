# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module provides wrappers, called Fitters, around some Numpy and Scipy
fitting functions. All Fitters take an instance of `~astropy.modeling.core.ParametricModel`
as input and define a ``__call__`` method which fits the model to the data and changes the
model's parameters attribute. The idea is to make this extensible and allow
users to easily add other fitters.

Linear fitting is done using Numpy's `~numpy.linalg.lstsq` function.
There are currently two non-linear fitters which use `~scipy.optimize.leastsq` and
`~scipy.optimize.slsqp` functions in scipy.optimize.\
"""

from __future__ import division

import abc
import numbers
import warnings

from functools import reduce

import numpy as np

from ..logger import log
from .utils import poly_map_domain
from ..utils.exceptions import AstropyUserWarning



__all__ = ['LinearLSQFitter', 'NonLinearLSQFitter', 'SLSQPFitter',
           'JointFitter', 'Fitter']



DEFAULT_MAXITER = 100
DEFAULT_EPS = np.sqrt(np.finfo(float).eps)
DEFAULT_MAX_BOUND = 10 ** 12
DEFAULT_MIN_BOUND = -10 ** 12


def _convert_input(x, y, z=None):
    x = np.asarray(x)
    y = np.asarray(y)
    if x.shape[0] != y.shape[0]:
        raise ValueError("x and y should have the same shape")
    if z is None:
        farg = (x, y)
    else:
        z = np.asarray(z)
        if x.shape != z.shape:
            raise ValueError("x, y and z should have the same shape")
        farg = (x, y, z)
    return farg


class ModelsError(Exception):
    """Base class for model exceptions"""


class ModelLinearityError(ModelsError):
    """
    Raised when a linear model is passed to a non-linear fitter and vice versa.
    """


class UnsupportedConstraintError(ModelsError, ValueError):
    """
    Raised when a fitter does not support a type of constraint.
    """


class Fitter(object):
    """
    Base class for all fitters.

    The purpose of this class is to manage constraints.
    """

    __metaclass__ = abc.ABCMeta

    # The base Fitter does not support any constraints by default; individual
    # fitters should explicitly set this list to the specific constraints
    # it supports
    # May contain any of 'bounds', 'eqcons', 'ineqcons', 'fixed', or 'tied'
    supported_constraints = []

    def __init__(self):
        self._weights = None

    def _validate_constraints(self, model):
        message = '{0} cannot handle {{0}} constraints.'.format(
                self.__class__.__name__)

        if (any(model.fixed.values()) and
                'fixed' not in self.supported_constraints):
            raise UnsupportedConstraintError(
                    message.format('fixed parameter'))

        if (any(model.tied.values()) and
                'tied' not in self.supported_constraints):
            raise UnsupportedConstraintError(
                    message.format('tied parameter'))

        if (any([tuple(b) != (None, None) for b in model.bounds.values()]) and
                'bounds' not in self.supported_constraints):
            raise UnsupportedConstraintError(
                    message.format('bound parameter'))

        if model.eqcons and 'eqcons' not in self.supported_constraints:
            raise UnsupportedConstraintError(message.format('equality'))

        if (model.ineqcons and
                'ineqcons' not in self.supported_constraints):
            raise UnsupportedConstraintError(message.format('inequality'))

    # TODO
    # @property
    # def covar(self):
    #     return None

    @abc.abstractmethod
    def __call__(self):
        """
        This method performs the actual fitting and modifies the parameter list
        of a model.

        Fitter subclasses should implement this method.
        """

        raise NotImplementedError("Subclasses should implement this")

    def _fitter_to_model_params(self, model, fps):
        _fit_params, _fit_param_indices = model._model_to_fit_params()
        if any(model.fixed.values()) or any(model.tied.values()):
            model.parameters[_fit_param_indices] = fps
            for idx, name in enumerate(model.param_names):
                if model.tied[name] != False:
                    value = model.tied[name](model)
                    slice_ = model._param_metrics[name][0]
                    model.parameters[slice_] = value
        elif any([tuple(b) != (None, None) for b in model.bounds.values()]):
            for name, par in zip(model.param_names, _fit_params):
                if model.bounds[name] != (None, None):
                    b = model.bounds[name]
                    if b[0] is not None:
                        par = max(par, model.bounds[name][0])
                    if b[1] is not None:
                        par = min(par, model.bounds[name][1])
                    setattr(model, name, par)
        else:
            model.parameters = fps


class LinearLSQFitter(Fitter):
    """
    A class performing a linear least square fitting.

    Uses `numpy.linalg.lstsq` to do the fitting.
    Given a model and data, fits the model to the data and changes the
    model's parameters. Keeps a dictionary of auxiliary fitting information.

    Parameters
    ----------
    model : an instance of `~astropy.modeling.core.ParametricModel`

    Raises
    ------
    ModelLinearityError
        A nonlinear model is passed to a linear fitter
    """

    supported_constraints = ['fixed']

    def __init__(self):
        super(LinearLSQFitter, self).__init__()
        self.fit_info = {'residuals': None,
                         'rank': None,
                         'singular_values': None,
                         'params': None
                         }

    @staticmethod
    def _deriv_with_constraints(model, param_indices, x=None, y=None):
        if y is None:
            d = np.array(model.deriv(x=x))
        else:
            d = np.array(model.deriv(x=x, y=y))

        if model.col_deriv:
            return d[param_indices]
        else:
            return d[:, param_indices]

    def _map_domain_window(self, model, x, y=None):
        """
        Maps domain into window for a polynomial model which has these
        attributes.
        """

        if y is None:
            if hasattr(model, 'domain') and model.domain is None:
                model.domain = [x.min(), x.max()]
            if hasattr(model, 'window') and model.window is None:
                model.window = [-1, 1]
            return poly_map_domain(x, model.domain, model.window)
        else:
            if hasattr(model, 'x_domain') and model.x_domain is None:
                model.x_domain = [x.min(), x.max()]
            if hasattr(model, 'y_domain') and model.y_domain is None:
                model.y_domain = [y.min(), y.max()]
            if hasattr(model, 'x_window') and model.x_window is None:
                model.x_window = [-1., 1.]
            if hasattr(model, 'y_window') and model.y_window is None:
                model.y_window = [-1., 1.]

            xnew = poly_map_domain(x, model.x_domain, model.x_window)
            ynew = poly_map_domain(y, model.y_domain, model.y_window)
            return xnew, ynew

    def __call__(self, model, x, y, z=None, weights=None, rcond=None):
        """
        Fit data to this model.

        Parameters
        ----------
        model : `ParametricModel`
            model to fit to x, y, z
        x : array
            input coordinates
        y : array
            input coordinates
        z : array (optional)
            input coordinates
        weights : array (optional)
            weights
        rcond :  float, optional
            Cut-off ratio for small singular values of `a`.
            Singular values are set to zero if they are smaller than `rcond`
            times the largest singular value of `a`.

        Returns
        ------
        model_copy : `ParametricModel`
            a copy of the input model with parameters set by the fitter
        """
        if not model.fittable:
            raise ValueError("Model must be a subclass of ParametricModel")
        if not model.linear:
            raise ModelLinearityError('Model is not linear in parameters, '
                                      'linear fit methods should not be used.')

        self._validate_constraints(model)
        multiple = False
        model_copy = model.copy()
        _, fitparam_indices = model_copy._model_to_fit_params()
        self._weights = weights
        if model_copy.n_inputs == 2 and z is None:
            raise ValueError("Expected x, y and z for a 2 dimensional model.")

        farg = _convert_input(x, y, z)

        if len(farg) == 2:
            x, y = farg
            if y.ndim == 2:
                assert y.shape[1] == model_copy.param_dim, (
                    "Number of data sets (Y array is expected to equal "
                    "the number of parameter sets")
            # map domain into window
            if hasattr(model_copy, 'domain'):
                x = self._map_domain_window(model_copy, x)
            if any(model_copy.fixed.values()):
                lhs = self._deriv_with_constraints(model_copy,
                                                   fitparam_indices,
                                                   x=x)
            else:
                lhs = model_copy.deriv(x=x)
            if len(y.shape) == 2:
                rhs = y
                multiple = y.shape[1]
            else:
                rhs = y
        else:
            x, y, z = farg
            if x.shape[-1] != z.shape[-1]:
                raise ValueError("x and z should have equal last dimensions")

            # map domain into window
            if hasattr(model_copy, 'x_domain'):
                x, y = self._map_domain_window(model_copy, x, y)

            if any(model_copy.fixed.values()):
                lhs = self._deriv_with_constraints(model_copy,
                                                   fitparam_indices, x=x, y=y)
            else:
                lhs = model_copy.deriv(x=x, y=y)
            if len(z.shape) == 3:
                rhs = np.array([i.flatten() for i in z]).T
                multiple = z.shape[0]
            else:
                rhs = z.flatten()

        if weights is not None:
            weights = np.asarray(weights, dtype=np.float)
            if len(x) != len(weights):
                raise ValueError("x and weights should have the same length")
            if rhs.ndim == 2:
                lhs *= weights[:, np.newaxis]
                rhs *= weights[:, np.newaxis]
            else:
                lhs *= weights[:, np.newaxis]
                rhs *= weights

        if not multiple and model_copy.param_dim > 1:
            raise ValueError("Attempting to fit a 1D data set to a model "
                             "with multiple parameter sets")
        if rcond is None:
            rcond = len(x) * np.finfo(x.dtype).eps

        scl = (lhs * lhs).sum(0)
        lacoef, resids, rank, sval = np.linalg.lstsq(lhs / scl, rhs, rcond)

        self.fit_info['residuals'] = resids
        self.fit_info['rank'] = rank
        self.fit_info['singular_values'] = sval

        # If y.n_inputs > model.n_inputs we are doing a simultanious 1D fitting
        # of several 1D arrays. Otherwise the model is 2D.
        # if y.n_inputs > self.model.n_inputs:
        if multiple and model_copy.param_dim != multiple:
            model_copy.param_dim = multiple
        # TODO: Changing the model's param_dim needs to be handled more
        # carefully; for now it's not actually allowed
        lacoef = (lacoef.T / scl).T
        self.fit_info['params'] = lacoef
        # TODO: Only Polynomial models currently have an _order attribute;
        # maybe change this to read isinstance(model, PolynomialBase)
        if hasattr(model_copy, '_order') and rank != model_copy._order:
            warnings.warn("The fit may be poorly conditioned\n",
                          AstropyUserWarning)
        self._fitter_to_model_params(model_copy, lacoef.flatten())
        return model_copy


class NonLinearLSQFitter(Fitter):
    """
    A class performing non-linear least squares fitting using the
    Levenberg-Marquardt algorithm implemented in `scipy.optimize.leastsq`.

    Parameters
    ----------
    model : a fittable `~astropy.modeling.core.ParametricModel`
        model to fit to data

    Raises
    ------
    ModelLinearityError
        A linear model is passed to a nonlinear fitter
    """

    supported_constraints = ['fixed', 'tied', 'bounds']

    def __init__(self):

        self.fit_info = {'nfev': None,
                         'fvec': None,
                         'fjac': None,
                         'ipvt': None,
                         'qtf': None,
                         'message': None,
                         'ierr': None,
                         'status': None}

        super(NonLinearLSQFitter, self).__init__()

    def errorfunc(self, fps, *args):
        model = args[0]
        self._fitter_to_model_params(model, fps)
        meas = args[-1]
        if self._weights is None:
            return np.ravel(model(*args[1 : -1]) - meas)
        else:
            return np.ravel(self._weights * (model(*args[1 : -1]) - meas))

    # @property
    # def covar(self):
    #     """
    #     Calculate the covariance matrix (doesn't take into account
    #     constraints).
    #     """

        # NOTE: Not currently used; won't work with Fitter implementation where
        # models are not tied to fitters

        # n = len(self.model.parameters)
        #  construct the permutation matrix
        # P = np.take(np.eye(n), self.fit_info['ipvt'] - 1, 0)
        #  construct the R matrix as in JP = QR
        # r = np.triu(self.fit_info['fjac'].T[:n, :])
        # r_pt = np.dot(r, P.T)
        # p_rt = np.dot(P, r.T)
        # try:
        #     return np.dual.inv(np.dot(p_rt, r_pt))
        # except:
        #     log.info("Could not construct a covariance matrix")
        #     return None

    def __call__(self, model, x, y, z=None, weights=None,
                 maxiter=DEFAULT_MAXITER,
                 epsilon=DEFAULT_EPS, estimate_jacobian=False):
        """
        Fit data to this model.

        Parameters
        ----------
        model : `ParametricModel`
            model to fit to x, y, z
        x : array
           input coordinates
        y : array
           input coordinates
        z : array (optional)
           input coordinates
        weights : array (optional
           weights
        maxiter : int
            maximum number of iterations
        epsilon : float
            A suitable step length for the forward-difference
            approximation of the Jacobian (if model.fjac=None). If
            epsfcn is less than the machine precision, it is
            assumed that the relative errors in the functions are
            of the order of the machine precision.
        estimate_jacobian : bool
            If False (default) and if the model has a deriv method,
            it will be used. Otherwise the Jacobian will be estimated.
            If True, the Jacobian will be estimated in any case.

        Returns
        ------
        model_copy : `ParametricModel`
            a copy of the input model with parameters set by the fitter
        """
        if not model.fittable:
            raise ValueError("Model must be a subclass of ParametricModel")
        self._validate_constraints(model)
        from scipy import optimize
        model_copy = model.copy()
        farg = (model_copy, ) + _convert_input(x, y, z)
        self._weights = weights
        if model_copy.param_dim != 1:
            # only single data sets ca be fitted
            raise ValueError("NonLinearLSQFitter can only fit one "
                             "data set at a time")
        if model_copy.deriv is None or estimate_jacobian:
            dfunc = None
        else:
            dfunc = self._wrap_deriv
        init_values, _ = model_copy._model_to_fit_params()
        fitparams, status, dinfo, mess, ierr = optimize.leastsq(
            self.errorfunc, init_values, args=farg, Dfun=dfunc,
            col_deriv=model_copy.col_deriv, maxfev=maxiter, epsfcn=epsilon,
            full_output=True)
        self._fitter_to_model_params(model_copy, fitparams)
        self.fit_info.update(dinfo)
        self.fit_info['status'] = status
        self.fit_info['message'] = mess
        self.fit_info['ierr'] = ierr
        if ierr not in [1, 2, 3, 4]:
            warnings.warn("The fit may be unsuccessful; check "
                          "fit_info['message'] for more information.",
                          AstropyUserWarning)
        return model_copy

    @staticmethod
    def _wrap_deriv(params, model, x, y, z=None):
        """
        Wraps the method calculating the Jacobian of the function to account
        for model constraints.

        Currently the only fitter that uses a derivative is the
        `NonLinearLSQFitter`. This wrapper may need to be revised when other
        fitters using function derivative are added or when the statistic is
        separated from the fitting routines.

        `~scipy.optimize.leastsq` expects the function derivative to have the
        above signature (parlist, (argtuple)). In order to accomodate model
        constraints, instead of using p directly, we set the parameter list in
        this function.
        """
        if any(model.fixed.values()) or any(model.tied.values()):

            if z is None:
                full_deriv = np.array(model.deriv(x, *model.parameters))
            else:
                full_deriv = np.array(model.deriv(x, y, *model.parameters))

            pars = [getattr(model, name) for name in model.param_names]
            fixed = [par.fixed for par in pars]
            tied = [par.tied for par in pars]
            tied = list(np.where([par.tied != False for par in pars], True, tied))
            fix_and_tie = np.logical_or(fixed, tied)
            ind = np.logical_not(fix_and_tie)

            if not model.col_deriv:
                full_deriv = np.asarray(full_deriv).T
                residues = np.asarray(full_deriv[np.nonzero(ind)])
            else:
                residues = full_deriv[np.nonzero(ind)]

            return [np.ravel(_) for _ in residues]
        else:
            if z is None:
                return model.deriv(x, *params)
            else:
                return [np.ravel(_) for _ in model.deriv(x, y, *params)]


class SLSQPFitter(Fitter):
    """
    Sequential Least Squares Programming optimization algorithm.

    The algorithm is described in [1]_. It supports tied and fixed
    parameters, as well as bounded constraints. Uses
    `scipy.optimize.fmin_slsqp`.

    Raises
    ------
    ModelLinearityError
        A linear model is passed to a nonlinear fitter

    References
    ----------
    .. [1] http://www.netlib.org/toms/733
    """

    supported_constraints = ['bounds', 'eqcons', 'ineqcons', 'fixed', 'tied']

    def __init__(self):
        super(SLSQPFitter, self).__init__()
        self.fit_info = {
            'final_func_val': None,
            'numiter': None,
            'exit_mode': None,
            'message': None
        }

    def errorfunc(self, fps, *args):
        """
        Compute the sum of the squared residuals

        Parameters
        ----------
        fps : list
            parameters returned by the fitter
        args : list
            input coordinates
        """
        model = args[0]
        meas = args[-1]
        self._fitter_to_model_params(model, fps)
        res = model(*args[1:-1]) - meas

        if self._weights is None:
            return np.sum(res ** 2)
        else:
            return np.sum(self._weights * res ** 2)

    def __call__(self, model, x, y, z=None, weights=None, verblevel=0,
                 maxiter=DEFAULT_MAXITER, epsilon=DEFAULT_EPS):
        """
        Fit data to this model.

        Parameters
        ----------
        model : `ParametricModel`
            model to fit to x, y, z
        x : array
            input coordinates
        y : array
            input coordinates
        z : array (optional)
            input coordinates
        weights : array (optional)
            weights
        verblevel : int
            0-silent
            1-print summary upon completion,
            2-print summary after each iteration
        maxiter : int
            maximum number of iterations
        epsilon : float
            the step size for finite-difference derivative estimates

        Returns
        ------
        model_copy : `ParametricModel`
            a copy of the input model with parameters set by the fitter
        """
        if not model.fittable:
            raise ValueError("Model must be a subclass of ParametricModel")
        if model.linear:
            warnings.warn('Model is linear in parameters; '
                          'consider using linear fitting methods.',
                          AstropyUserWarning)
        self._validate_constraints(model)
        model_copy = model.copy()
        from scipy import optimize

        farg = _convert_input(x, y, z)
        farg = (model_copy, ) + farg
        self._weights = weights
        if model_copy.param_dim != 1:
            # for now only single data sets ca be fitted
            raise ValueError("NonLinearLSQFitter can only fit "
                             "one data set at a time")

        p0, param_indices = model_copy._model_to_fit_params()
        pars = [getattr(model_copy, name) for name in model_copy.param_names]
        bounds = [par.bounds for par in pars if par.fixed != True and par.tied == False]
        
        bounds = np.asarray(bounds)
        for i in bounds:
            if i[0] is None:
                i[0] = DEFAULT_MIN_BOUND
            if i[1] is None:
                i[1] = DEFAULT_MAX_BOUND
        # older versions of scipy require this array to be float
        bounds = np.asarray(bounds, dtype=np.float)
        eqcons = np.array(model_copy.eqcons)
        ineqcons = np.array(model_copy.ineqcons) 
        fitparams, final_func_val, numiter, exit_mode, mess = \
            optimize.fmin_slsqp(
            self.errorfunc, p0, args=farg, disp=verblevel, full_output=1,
            bounds=bounds, eqcons=eqcons, ieqcons=ineqcons, iter=maxiter,
            acc=1.E-6, epsilon=DEFAULT_EPS)

        self._fitter_to_model_params(model_copy, fitparams)
        self.fit_info['final_func_val'] = final_func_val
        self.fit_info['numiter'] = numiter
        self.fit_info['exit_mode'] = exit_mode
        self.fit_info['message'] = mess

        if exit_mode != 0:
            warnings.warn("The fit may be unsuccessful; check "
                          "fit_info['message'] for more information.",
                          AstropyUserWarning)
        return model_copy


class JointFitter(object):
    """
    Fit models which share a parameter.

    For example, fit two gaussians to two data sets but keep
    the FWHM the same.

    Parameters
    ----------
    models : list
        a list of model instances
    jointparameters : list
        a list of joint parameters
    initvals : list
        a list of initial values
    """

    def __init__(self, models, jointparameters, initvals):
        self.models = list(models)
        self.initvals = list(initvals)
        self.jointparams = jointparameters
        self._verify_input()
        for m in self.jointparams.keys():
            m.set_joint_parameters(self.jointparams[m])
        self.fitparams = self._model_to_fit_params()

        # a list of model.n_inputs
        self.modeldims = [m.n_inputs for m in self.models]
        # sum all model dimensions
        self.ndim = np.sum(self.modeldims)

    def _model_to_fit_params(self):
        fparams = []
        fparams.extend(self.initvals)
        for model in self.models:
            params = [p.flatten() for p in model.parameters]
            for pname in model.joint:
                slc = model._param_metrics[pname][0]
                del params[slc]
            fparams.extend(params)
        return fparams

    def errorfunc(self, fps, *args):
        """
        fps : list
            the fitted parameters - result of an one iteration of the
            fitting algorithm
        args : dict
            tuple of measured and input coordinates
            args is always passed as a tuple from optimize.leastsq
        """

        lstsqargs = list(args)
        fitted = []
        fitparams = list(fps)
        numjp = len(self.initvals)
        # make a separate list of the joint fitted parameters
        jointfitparams = fitparams[:numjp]
        del fitparams[:numjp]

        for model in self.models:
            margs = lstsqargs[:model.n_inputs + 1]
            del lstsqargs[:model.n_inputs + 1]
            # separate each model separately fitted parameters
            numfp = len(model._parameters) - len(model.joint)
            mfparams = fitparams[:numfp]

            del fitparams[:numfp]
            # recreate the model parameters
            mparams = []
            for pname in model.param_names:
                if pname in model.joint:
                    index = model.joint.index(pname)
                    # should do this with slices in case the
                    # parameter is not a number
                    mparams.extend([jointfitparams[index]])
                else:
                    slc = model._param_metrics[pname][0]
                    plen = slc.stop - slc.start
                    mparams.extend(mfparams[:plen])
                    del mfparams[:plen]
            modelfit = model.eval(margs[:-1], *mparams)
            fitted.extend(modelfit - margs[-1])
        return np.ravel(fitted)

    def _verify_input(self):
        assert(len(self.models) > 1)
        assert(len(self.jointparams.keys()) >= 2)
        for j in self.jointparams.keys():
            assert(len(self.jointparams[j]) == len(self.initvals))

    def __call__(self, *args):
        """
        Fit data to these models keeping some of the pramaters common to the
        two models.
        """

        from scipy import optimize

        assert(len(args) == reduce(lambda x, y: x + 1 + y + 1, self.modeldims))
        self.fitparams[:], _ = optimize.leastsq(self.errorfunc, self.fitparams,
                                              args=args)

        fparams = self.fitparams[:]
        numjp = len(self.initvals)
        # make a separate list of the joint fitted parameters
        jointfitparams = fparams[:numjp]
        del fparams[:numjp]

        for model in self.models:
            # extract each model's fitted parameters
            numfp = len(model._parameters) - len(model.joint)
            mfparams = fparams[:numfp]

            del fparams[:numfp]
            # recreate the model parameters
            mparams = []
            for pname in model.param_names:
                if pname in model.joint:
                    index = model.joint.index(pname)
                    # should do this with slices in case the parameter
                    # is not a number
                    mparams.extend([jointfitparams[index]])
                else:
                    slc = model._param_metrics[pname][0]
                    plen = slc.stop - slc.start
                    mparams.extend(mfparams[:plen])
                    del mfparams[:plen]
            model.parameters = np.array(mparams)
