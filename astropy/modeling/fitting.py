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


MAXITER = 100
EPS = np.sqrt(np.finfo(float).eps)


# supported constraints
# TODO: Register these automatically through a metaclass
constraintsdef = {
    'NonLinearLSQFitter': ['fixed', 'tied', 'bounds'],
    'SLSQPFitter': ['bounds', 'eqcons', 'ineqcons', 'fixed', 'tied'],
    'LinearLSQFitter': ['fixed']}


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


class UnsupportedConstraintError(ModelsError):
    """
    Raised when a fitter does not support a type of constraint.
    """


class Fitter(object):
    """
    Base class for all fitters.

    The purpose of this class is to manage constraints.
    """

    __metaclass__ = abc.ABCMeta

    def __init__(self, model):
        if not model.fittable:
            raise ValueError("Model must be a subclass of ParametricModel")

        self._model = model
        self.bounds = []
        self.fixed = []
        self.tied = []
        self._set_constraints()
        self._fitparams = self._model_to_fit_params()
        self._validate_constraints()
        self._weights = None

    def _set_constraints(self):
        pars = [getattr(self.model, name) for name in self.model.param_names]
        self.fixed = [par.fixed for par in pars]
        self.tied = [par.tied for par in pars]
        min_values = [par.min for par in pars]
        max_values = [par.max for par in pars]
        b = []
        for i, j in zip(min_values, max_values):
            if i is None:
                i = -10 ** 12
            if j is None:
                j = 10 ** 12
            b.append((i, j))
        self.bounds = b[:]

    @property
    def model(self):
        """The model being fitted."""

        return self._model

    @model.setter
    def model(self, val):
        self._model = val

    @property
    def fitparams(self):
        return self._fitparams

    @fitparams.setter
    def fitparams(self, fps):
        """
        Update ``.model.parameters`` from fitparams in the presence of
        constraints.

        This is the opposite of ``_model_to_fit_params``.

        Parameters
        ----------
        fps : list
            list of parameters, fitted in a succesive iteration of the
            fitting algorithm
        """

        self._fitparams[:] = fps
        if any(self.fixed) or any(self.tied):
            fitparams = list(fps[:])
            mpars = []
            for i, name in enumerate(self.model.param_names):
                if self.fixed[i] is True:
                    par = getattr(self.model, name)
                    if len(par.shape) == 0:
                        mpars.extend([par.value])
                    else:
                        mpars.extend(par.value)
                elif self.tied[i] is not False:
                    val = self.tied[i](self.model)
                    if isinstance(val, numbers.Number):
                        mpars.append(val)
                    else:
                        mpars.extend(val)
                else:
                    sl = self.model._parameters.parinfo[name][0]
                    plen = sl.stop - sl.start
                    mpars.extend(fitparams[:plen])
                    del fitparams[:plen]
            self.model.parameters = mpars
        elif any([b != (-1E12, 1E12) for b in self.bounds]):
            self._set_bounds(fps)
        else:
            self.model.parameters[:] = fps

    def _model_to_fit_params(self):
        """
        Create a set of parameters to be fitted.

        These may be a subset of the model parameters, if some of them are held
        constant or tied.
        """

        if any(self.model.fixed.values()) or any(self.model.tied.values()):
            pars = self.model._parameters[:]
            for item in self.model.param_names[::-1]:
                if self.model.fixed[item] or self.model.tied[item]:
                    sl = self.model._parameters.parinfo[item][0]
                    del pars[sl]
            return pars
        else:
            return self.model._parameters

    def _set_bounds(self, pars):
        """
        Diferent fitting algorithms deal with bounds in a different way.  For
        example, the SLSQP algorithm accepts bounds as input while the leastsq
        algorithm does not handle bounds at all and they are dealt with in a
        separate method.

        This method is to be implemented by subcclasses of Fitter if necessary.
        """

        raise NotImplementedError("Subclasses should implement this")

    def _wrap_deriv(self, p, x, y, z=None):
        """
        Wraps the method calculating the Jacobian of the function to account
        for model constraints.

        Currently the only fitter that uses a derivative is the
        `NonLinearLSQFitter`. This wrapper may neeed to be revised when other
        fitters using function derivative are added or when the statistic is
        separated from the fitting routines.

        `~scipy.optimize.leastsq` expects the function derivative to have the
        above signature (parlist, (argtuple)). In order to accomodate model
        constraints, instead of using p directly, we set the parameter list in
        this function.
        """

        if any(self.fixed) or any(self.tied):
            params = self.model._parameters

            if z is None:
                full_deriv = np.array(self.model.deriv(x, *params))
            else:
                full_deriv = np.array(self.model.deriv(x, y, *params))

            ind = range(len(self.model.param_names))
            fix_tie_ind = list(np.nonzero(self.fixed)[0])
            fix_tie_ind.extend(np.nonzero(self.tied)[0])

            for index in fix_tie_ind:
                ind.remove(index)

            res = np.empty((full_deriv.shape[0],
                            full_deriv.shape[1] - len(ind)))
            res = full_deriv[ind, :]
            return [np.ravel(_) for _ in res]
        else:
            params = p[:]
            if z is None:
                return self.model.deriv(x, *params)
            else:
                return [np.ravel(_) for _ in self.model.deriv(x, y, *params)]

    def _validate_constraints(self):
        fname = self.__class__.__name__
        try:
            c = constraintsdef[fname]
        except KeyError:
            raise UnsupportedConstraintError("{0} does not support fitting",
                                             "with constraints".format(fname))
        if any(self.fixed) and 'fixed' not in c:
            raise ValueError("{0} cannot handle fixed parameter",
                             "constraints .".format(fname))
        if any(self.tied) and 'tied' not in c:
            raise ValueError("{0} cannot handle tied parameter",
                             "constraints ".format(fname))
        if any([b != (-1E12, 1E12) for b in self.bounds]) and 'bounds' not in c:
            raise ValueError("{0} cannot handle bound parameter",
                             "constraints".format(fname))
        if self.model.eqcons and 'eqcons' not in c:
            raise ValueError("{0} cannot handle equality constraints but ",
                             "eqcons given".format(fname))
        if self.model.ineqcons and 'ineqcons' not in c:
            raise ValueError("{0} cannot handle inequality constraints but ",
                             "ineqcons given".format(fname))

    @property
    def covar(self):
        return None

    @property
    def weights(self):
        """Fitting weights."""

        return self._weights

    @weights.setter
    def weights(self, val):
        """Set fitting weights."""

        self._weights = val

    def _update_constraints(self):
        self._set_constraints()
        self._fitparams = self._model_to_fit_params()

    @abc.abstractmethod
    def __call__(self):
        """
        This method performs the actual fitting and modifies the parameter list
        of a model.

        Fitter subclasses should implement this method.
        """

        raise NotImplementedError("Subclasses should implement this")


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

    def __init__(self, model):
        super(LinearLSQFitter, self).__init__(model)
        if not self.model.linear:
            raise ModelLinearityError('Model is not linear in parameters, '
                                      'linear fit methods should not be used.')
        self.fit_info = {'residuals': None,
                         'rank': None,
                         'singular_values': None,
                         'params': None
                         }

    def _deriv_with_constraints(self, params=None, x=None, y=None):
        if y is None:
            d = self.model.deriv(x=x)
        else:
            d = self.model.deriv(x=x, y=y)
        fixed = [name for name in self.model.fixed if
                 self.model.fixed[name]]
        ind = range(len(self.model.param_names))
        for name in fixed:
            index = self.model.param_names.index(name)
            ind.remove(index)
        res = d[:, ind]
        return res

    def _map_domain_window(self, x, y=None):
        """
        Maps domain into window for a polynomial model which has these
        attributes.
        """

        if y is None:
            if hasattr(self.model, 'domain') and self.model.domain is None:
                self.model.domain = [x.min(), x.max()]
            if hasattr(self.model, 'window') and self.model.window is None:
                self.model.window = [-1, 1]
            return poly_map_domain(x, self.model.domain, self.model.window)
        else:
            if hasattr(self.model, 'x_domain') and self.model.x_domain is None:
                self.model.x_domain = [x.min(), x.max()]
            if hasattr(self.model, 'y_domain') and self.model.y_domain is None:
                self.model.y_domain = [y.min(), y.max()]
            if hasattr(self.model, 'x_window') and self.model.x_window is None:
                self.model.x_window = [-1., 1.]
            if hasattr(self.model, 'y_window') and self.model.y_window is None:
                self.model.y_window = [-1., 1.]

            xnew = poly_map_domain(x, self.model.x_domain, self.model.x_window)
            ynew = poly_map_domain(y, self.model.y_domain, self.model.y_window)
            return xnew, ynew

    def __call__(self, x, y, z=None, weights=None, rcond=None):
        """
        Fit data to this model.

        Parameters
        ----------
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
        """

        super(LinearLSQFitter, self)._update_constraints()

        multiple = False

        if self.model.n_inputs == 2 and z is None:
            raise ValueError("Expected x, y and z for a 2 dimensional model.")

        farg = _convert_input(x, y, z)

        if len(farg) == 2:
            x, y = farg
            if y.ndim == 2:
                assert y.shape[1] == self.model._parameters.dim, (
                    "Number of data sets (Y array is expected to equal "
                    "the number of parameter sets")
            # map domain into window
            if hasattr(self.model, 'domain'):
                x = self._map_domain_window(x)
            if any(self.model.fixed.values()):
                lhs = self._deriv_with_constraints(x=x)
            else:
                lhs = self.model.deriv(x=x)
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
            if hasattr(self.model, 'x_domain'):
                x, y = self._map_domain_window(x, y)

            if any(self.model.fixed.values()):
                lhs = self._deriv_with_constraints(x=x, y=y)
            else:
                lhs = self.model.deriv(x=x, y=y)
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

        if not multiple and self.model._parameters.dim > 1:
            raise ValueError("Attempting to fit a 1D data set to a model "
                             "with multiple parameter sets")
        if rcond is None:
            rcond = len(x) * np.finfo(x.dtype).eps

        scl = (lhs * lhs).sum(0)
        lacoef, resids, rank, sval = np.linalg.lstsq(lhs / scl, rhs, rcond)

        self.fit_info['residuals'] = resids
        self.fit_info['rank'] = rank
        self.fit_info['singular_values'] = sval

        self.model._parameters._changed = True
        # If y.n_inputs > model.n_inputs we are doing a simultanious 1D fitting
        # of several 1D arrays. Otherwise the model is 2D.
        # if y.n_inputs > self.model.n_inputs:
        if multiple:
            self.model._parameters.dim = multiple
        lacoef = (lacoef.T / scl).T
        self.fit_info['params'] = lacoef
        if rank != self.model._order:
            warnings.warn("The fit may be poorly conditioned\n",
                          AstropyUserWarning)
        self.fitparams = lacoef.flatten()[:]


class NonLinearLSQFitter(Fitter):
    """
    A class performing non-linear least squares fitting using the
    Levenberg-Marquardt algorithm implemented in `scipy.optimize.leastsq`.

    Parameters
    ----------
    model : a fittable :class: `~astropy.modeling.core.ParametricModel`
        model to fit to data

    Raises
    ------
    ModelLinearityError
        A linear model is passed to a nonlinear fitter
    """

    def __init__(self, model):

        self.fit_info = {'nfev': None,
                         'fvec': None,
                         'fjac': None,
                         'ipvt': None,
                         'qtf': None,
                         'message': None,
                         'ierr': None,
                         'status': None}

        super(NonLinearLSQFitter, self).__init__(model)
        if self.model.linear:
            warnings.warn('Model is linear in parameters, '
                          'consider using linear fitting methods.', AstropyUserWarning)

    def errorfunc(self, fps, *args):
        self.fitparams = fps
        meas = args[-1]
        if self.weights is None:
            return np.ravel(self.model(*args[: -1]) - meas)
        else:
            return np.ravel(self.weights * (self.model(*args[: -1]) - meas))

    def _set_bounds(self, fitparams):
        for c in self.model.bounds.values():
            if c != (-1E12, 1E12):
                for name, par, b in zip(self.model.param_names, fitparams,
                                        self.bounds):
                    par = max(par, b[0])
                    par = min(par, b[1])
                    setattr(self.model, name, par)

    @property
    def covar(self):
        """
        Calculate the covariance matrix (doesn't take into account
        constraints).
        """

        n = len(self.model.parameters)
        # construct the permutation matrix
        P = np.take(np.eye(n), self.fit_info['ipvt'] - 1, 0)
        # construct the R matrix as in JP = QR
        r = np.triu(self.fit_info['fjac'].T[:n, :])
        r_pt = np.dot(r, P.T)
        p_rt = np.dot(P, r.T)
        try:
            return np.dual.inv(np.dot(p_rt, r_pt))
        except:
            log.info("Could not construct a covariance matrix")
            return None

    def __call__(self, x, y, z=None, weights=None, maxiter=MAXITER,
                 epsilon=EPS, estimate_jacobian=False):
        """
        Fit data to this model.

        Parameters
        ----------
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
        """

        from scipy import optimize

        self._update_constraints()
        farg = _convert_input(x, y, z)
        self.weights = weights
        if self.model._parameters.dim != 1:
            # for now only single data sets ca be fitted
            raise ValueError("NonLinearLSQFitter can only fit one "
                             "data set at a time")
        if self._model.deriv is None or estimate_jacobian:
            self.dfunc = None
        else:
            self.dfunc = self._wrap_deriv

        self.fitparams, status, dinfo, mess, ierr = optimize.leastsq(
            self.errorfunc, self.fitparams, args=farg, Dfun=self.dfunc,
            col_deriv=self.model.col_deriv, maxfev=maxiter, epsfcn=epsilon,
            full_output=True)

        self.fit_info.update(dinfo)
        self.fit_info['status'] = status
        self.fit_info['message'] = mess
        self.fit_info['ierr'] = ierr
        if ierr not in [1, 2, 3, 4]:
            warnings.warn("The fit may be unsuccessful; check "
                          "fit_info['message'] for more information.",
                          AstropyUserWarning)


class SLSQPFitter(Fitter):
    """
    Sequential Least Squares Programming optimization algorithm.

    The algorithm is described in [1]_. It supports tied and fixed
    parameters, as well as bounded constraints. Uses
    `scipy.optimize.slsqp`.

    Parameters
    ----------
    model : a fittable :class: `models.ParametricModel`
        model to fit to data

    Raises
    ------
    ModelLinearityError
        A linear model is passed to a nonlinear fitter

    References
    ----------
    .. [1] http://www.netlib.org/toms/733
    """

    def __init__(self, model):
        super(SLSQPFitter, self).__init__(model)
        if self.model.linear:
            warnings.warn('Model is linear in parameters, '
                          'consider using linear fitting methods.', AstropyUserWarning)

        self.fit_info = {'final_func_val': None,
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

        meas = args[-1]
        self.fitparams = fps
        res = self.model(*args[:-1]) - meas

        if self.weights is None:
            return np.sum(res ** 2)
        else:
            return np.sum(self.weights * res ** 2)

    def _set_bounds(self, fitparams):
        """
        Set this as a dummy method because the SLSQP fitter
        handles bounds internally.
        """
        self._fitparams[:] = fitparams
        self.model.parameters = fitparams

    def __call__(self, x, y, z=None, weights=None, verblevel=0,
                 maxiter=MAXITER, epsilon=EPS):
        """
        Fit data to this model.

        Parameters
        ----------
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
        """

        from scipy import optimize

        # update constraints to pick up changes to parameters
        self._update_constraints()
        farg = _convert_input(x, y, z)
        self.weights = weights
        if self.model._parameters.dim != 1:
            # for now only single data sets ca be fitted
            raise ValueError("NonLinearLSQFitter can only fit "
                             "one data set at a time")

        p0 = self.model._parameters[:]

        self.fitparams, final_func_val, numiter, exit_mode, mess = \
            optimize.fmin_slsqp(
                self.errorfunc, p0, args=farg, disp=verblevel, full_output=1,
                bounds=self.bounds, eqcons=self.model.eqcons,
                ieqcons=self.model.ineqcons, iter=maxiter, acc=1.E-6,
                epsilon=EPS)

        self.fit_info['final_func_val'] = final_func_val
        self.fit_info['numiter'] = numiter
        self.fit_info['exit_mode'] = exit_mode
        self.fit_info['message'] = mess

        if exit_mode != 0:
            warnings.warn("The fit may be unsuccessful; check "
                          "fit_info['message'] for more information.",
                          AstropyUserWarning)


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
            params = model._parameters[:]
            for pname in model.joint:
                sl = model._parameters.parinfo[pname][0]
                del params[sl]
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
                    sl = model._parameters.parinfo[pname][0]
                    plen = sl.stop - sl.start
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
                    sl = model._parameters.parinfo[pname][0]
                    plen = sl.stop - sl.start
                    mparams.extend(mfparams[:plen])
                    del mfparams[:plen]
            model._parameters[:] = np.array(mparams)
