# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module implements classes (called Fitters) which combine optimization
algorithms (typically from `scipy.optimize`) with statistic functions to perform
fitting. Fitters are implemented as callable classes. In addition to the data
to fit, the ``__call__`` method takes an instance of
`~astropy.modeling.core.FittableModel` as input, and returns a copy of the
model with its parameters determined by the optimizer.

Optimization algorithms, called "optimizers" are implemented in
`~astropy.modeling.optimizers` and statistic functions are in
`~astropy.modeling.statistic`. The goal is to provide an easy to extend
framework and allow users to easily create new fitters by combining statistics
with optimizers.

There are two exceptions to the above scheme.
`~astropy.modeling.fitting.LinearLSQFitter` uses Numpy's `~numpy.linalg.lstsq`
function.  `~astropy.modeling.fitting.LevMarLSQFitter` uses
`~scipy.optimize.leastsq` which combines optimization and statistic in one
implementation.
"""

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import abc
import inspect
import operator
import warnings

from functools import reduce, wraps

import numpy as np

from .utils import poly_map_domain, _combine_equivalency_dict
from ..units import Quantity
from ..utils.exceptions import AstropyUserWarning
from ..extern import six
from ..extern.six.moves import range
from .optimizers import (SLSQP, Simplex)
from .statistic import (leastsquare)

# Check pkg_resources exists
try:
    from pkg_resources import iter_entry_points
    HAS_PKG = True
except ImportError:
    HAS_PKG = False


__all__ = ['LinearLSQFitter', 'LevMarLSQFitter', 'FittingWithOutlierRemoval',
           'SLSQPLSQFitter', 'SimplexLSQFitter', 'JointFitter', 'Fitter']


# Statistic functions implemented in `astropy.modeling.statistic.py
STATISTICS = [leastsquare]

# Optimizers implemented in `astropy.modeling.optimizers.py
OPTIMIZERS = [Simplex, SLSQP]

from .optimizers import (DEFAULT_MAXITER, DEFAULT_EPS, DEFAULT_ACC)


class ModelsError(Exception):
    """Base class for model exceptions"""


class ModelLinearityError(ModelsError):
    """ Raised when a non-linear model is passed to a linear fitter."""


class UnsupportedConstraintError(ModelsError, ValueError):
    """
    Raised when a fitter does not support a type of constraint.
    """


class _FitterMeta(abc.ABCMeta):
    """
    Currently just provides a registry for all Fitter classes.
    """

    registry = set()

    def __new__(mcls, name, bases, members):
        cls = super(_FitterMeta, mcls).__new__(mcls, name, bases, members)

        if not inspect.isabstract(cls) and not name.startswith('_'):
            mcls.registry.add(cls)

        return cls


def fitter_unit_support(func):
    """
    This is a decorator that can be used to add support for dealing with
    quantities to any __call__ method on a fitter which may not support
    quantities itself. This is done by temporarily removing units from all
    parameters then adding them back once the fitting has completed.
    """
    @wraps(func)
    def wrapper(self, model, x, y, z=None, **kwargs):
        equivalencies = kwargs.pop('equivalencies', None)

        data_has_units = (isinstance(x, Quantity) or
                          isinstance(y, Quantity) or
                          isinstance(z, Quantity))

        model_has_units = model._has_units

        if data_has_units or model_has_units:

            if model._supports_unit_fitting:

                # We now combine any instance-level input equivalencies with user
                # specified ones at call-time.


                input_units_equivalencies = _combine_equivalency_dict(
                    model.inputs, equivalencies, model.input_units_equivalencies)

                # If input_units is defined, we transform the input data into those
                # expected by the model. We hard-code the input names 'x', and 'y'
                # here since FittableModel instances have input names ('x',) or
                # ('x', 'y')

                if model.input_units is not None:
                    if isinstance(x, Quantity):
                        x = x.to(model.input_units['x'], equivalencies=input_units_equivalencies['x'])
                    if isinstance(y, Quantity) and z is not None:
                        y = y.to(model.input_units['y'], equivalencies=input_units_equivalencies['y'])

                # We now strip away the units from the parameters, taking care to
                # first convert any parameters to the units that correspond to the
                # input units (to make sure that initial guesses on the parameters)
                # are in the right unit system

                model = model.without_units_for_data(x=x, y=y, z=z)

                # We strip away the units from the input itself

                add_back_units = False

                if isinstance(x, Quantity):
                    add_back_units = True
                    xdata = x.value
                else:
                    xdata = np.asarray(x)

                if isinstance(y, Quantity):
                    add_back_units = True
                    ydata = y.value
                else:
                    ydata = np.asarray(y)

                if z is not None:
                    if isinstance(y, Quantity):
                        add_back_units = True
                        zdata = z.value
                    else:
                        zdata = np.asarray(z)

                # We run the fitting
                if z is None:
                    model_new = func(self, model, xdata, ydata, **kwargs)
                else:
                    model_new = func(self, model, xdata, ydata, zdata, **kwargs)

                # And finally we add back units to the parameters
                if add_back_units:
                    model_new = model_new.with_units_from_data(x=x, y=y, z=z)

                return model_new

            else:

                raise NotImplementedError("This model does not support being fit to data with units")

        else:

            return func(self, model, x, y, z=z, **kwargs)

    return wrapper


@six.add_metaclass(_FitterMeta)
class Fitter(object):
    """
    Base class for all fitters.

    Parameters
    ----------
    optimizer : callable
        A callable implementing an optimization algorithm
    statistic : callable
        Statistic function
    """

    def __init__(self, optimizer, statistic):
        if optimizer is None:
            raise ValueError("Expected an optimizer.")
        if statistic is None:
            raise ValueError("Expected a statistic function.")
        if inspect.isclass(optimizer):
            # a callable class
            self._opt_method = optimizer()
        elif inspect.isfunction(optimizer):
            self._opt_method = optimizer
        else:
            raise ValueError("Expected optimizer to be a callable class or a function.")
        if inspect.isclass(statistic):
            self._stat_method = statistic()
        else:
            self._stat_method = statistic

    def objective_function(self, fps, *args):
        """
        Function to minimize.

        Parameters
        ----------
        fps : list
            parameters returned by the fitter
        args : list
            [model, [other_args], [input coordinates]]
            other_args may include weights or any other quantities specific for
            a statistic

        Notes
        -----
        The list of arguments (args) is set in the `__call__` method.
        Fitters may overwrite this method, e.g. when statistic functions
        require other arguments.

        """
        model = args[0]
        meas = args[-1]
        _fitter_to_model_params(model, fps)
        res = self._stat_method(meas, model, *args[1:-1])
        return res

    @abc.abstractmethod
    def __call__(self):
        """
        This method performs the actual fitting and modifies the parameter list
        of a model.

        Fitter subclasses should implement this method.
        """

        raise NotImplementedError("Subclasses should implement this method.")


# TODO: I have ongoing branch elsewhere that's refactoring this module so that
# all the fitter classes in here are Fitter subclasses.  In the meantime we
# need to specify that _FitterMeta is its metaclass.
@six.add_metaclass(_FitterMeta)
class LinearLSQFitter(object):
    """
    A class performing a linear least square fitting.

    Uses `numpy.linalg.lstsq` to do the fitting.
    Given a model and data, fits the model to the data and changes the
    model's parameters. Keeps a dictionary of auxiliary fitting information.
    """

    supported_constraints = ['fixed']

    def __init__(self):
        self.fit_info = {'residuals': None,
                         'rank': None,
                         'singular_values': None,
                         'params': None
                         }

    @staticmethod
    def _deriv_with_constraints(model, param_indices, x=None, y=None):
        if y is None:
            d = np.array(model.fit_deriv(x, *model.parameters))
        else:
            d = np.array(model.fit_deriv(x, y, *model.parameters))

        if model.col_fit_deriv:
            return d[param_indices]
        else:
            return d[..., param_indices]

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

    @fitter_unit_support
    def __call__(self, model, x, y, z=None, weights=None, rcond=None):
        """
        Fit data to this model.

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
        weights : array (optional)
            Weights for fitting.
            For data with Gaussian uncertainties, the weights should be
            1/sigma.
        rcond :  float, optional
            Cut-off ratio for small singular values of ``a``.
            Singular values are set to zero if they are smaller than ``rcond``
            times the largest singular value of ``a``.
        equivalencies : list or None, optional and keyword-only argument
            List of *additional* equivalencies that are should be applied in
            case x, y and/or z have units. Default is None.

        Returns
        -------
        model_copy : `~astropy.modeling.FittableModel`
            a copy of the input model with parameters set by the fitter
        """

        if not model.fittable:
            raise ValueError("Model must be a subclass of FittableModel")

        if not model.linear:
            raise ModelLinearityError('Model is not linear in parameters, '
                                      'linear fit methods should not be used.')

        _validate_constraints(self.supported_constraints, model)

        model_copy = model.copy()
        _, fitparam_indices = _model_to_fit_params(model_copy)

        if model_copy.n_inputs == 2 and z is None:
            raise ValueError("Expected x, y and z for a 2 dimensional model.")

        farg = _convert_input(x, y, z, n_models=len(model_copy),
                              model_set_axis=model_copy.model_set_axis)

        has_fixed = any(model_copy.fixed.values())

        if has_fixed:

            # The list of fixed params is the complement of those being fitted:
            fixparam_indices = [idx for idx in
                                range(len(model_copy.param_names))
                                if idx not in fitparam_indices]

            # Construct matrix of user-fixed parameters that can be dotted with
            # the corresponding fit_deriv() terms, to evaluate corrections to
            # the dependent variable in order to fit only the remaining terms:
            fixparams = np.asarray([getattr(model_copy,
                                            model_copy.param_names[idx]).value
                                    for idx in fixparam_indices])

        if len(farg) == 2:
            x, y = farg

            # map domain into window
            if hasattr(model_copy, 'domain'):
                x = self._map_domain_window(model_copy, x)
            if has_fixed:
                lhs = self._deriv_with_constraints(model_copy,
                                                   fitparam_indices,
                                                   x=x)
                fixderivs = self._deriv_with_constraints(model_copy,
                                                         fixparam_indices,
                                                         x=x)
            else:
                lhs = model_copy.fit_deriv(x, *model_copy.parameters)
            sum_of_implicit_terms = model_copy.sum_of_implicit_terms(x)
            rhs = y
        else:
            x, y, z = farg

            # map domain into window
            if hasattr(model_copy, 'x_domain'):
                x, y = self._map_domain_window(model_copy, x, y)

            if has_fixed:
                lhs = self._deriv_with_constraints(model_copy,
                                                   fitparam_indices, x=x, y=y)
                fixderivs = self._deriv_with_constraints(model_copy,
                                                    fixparam_indices, x=x, y=y)
            else:
                lhs = model_copy.fit_deriv(x, y, *model_copy.parameters)
            sum_of_implicit_terms = model_copy.sum_of_implicit_terms(x, y)

            if len(model_copy) > 1:
                if z.ndim > 2:
                    # Basically this code here is making the assumption that if
                    # z has 3 dimensions it represents multiple models where
                    # the value of z is one plane per model.  It's then
                    # flattening each plane and transposing so that the model
                    # axis is *last*.  That's fine, but this could be
                    # generalized for other dimensionalities of z.
                    # TODO: See above comment
                    rhs = np.array([i.flatten() for i in z]).T
                else:
                    rhs = z.T
            else:
                rhs = z.flatten()

        # If the derivative is defined along rows (as with non-linear models)
        if model_copy.col_fit_deriv:
            lhs = np.asarray(lhs).T

        # Subtract any terms fixed by the user from (a copy of) the RHS, in
        # order to fit the remaining terms correctly:
        if has_fixed:
            if model_copy.col_fit_deriv:
                fixderivs = np.asarray(fixderivs).T  # as for lhs above
            rhs = rhs - fixderivs.dot(fixparams)  # evaluate user-fixed terms

        # Subtract any terms implicit in the model from the RHS, which, like
        # user-fixed terms, affect the dependent variable but are not fitted:
        if sum_of_implicit_terms is not None:
            # If we have a model set, the extra axis must be added to
            # sum_of_implicit_terms as its innermost dimension, to match the
            # dimensionality of rhs after _convert_input "rolls" it as needed
            # by np.linalg.lstsq. The vector then gets broadcast to the right
            # number of sets (columns). This assumes all the models share the
            # same input co-ordinates, as is currently the case.
            if len(model_copy) > 1:
                sum_of_implicit_terms = sum_of_implicit_terms[..., np.newaxis]
            rhs = rhs - sum_of_implicit_terms

        if weights is not None:
            weights = np.asarray(weights, dtype=np.float)
            if len(x) != len(weights):
                raise ValueError("x and weights should have the same length")
            if rhs.ndim == 2:
                lhs *= weights[:, np.newaxis]
                # Don't modify in-place in case rhs was the original dependent
                # variable array
                rhs = rhs * weights[:, np.newaxis]
            else:
                lhs *= weights[:, np.newaxis]
                rhs = rhs * weights

        if rcond is None:
            rcond = len(x) * np.finfo(x.dtype).eps

        scl = (lhs * lhs).sum(0)
        lacoef, resids, rank, sval = np.linalg.lstsq(lhs / scl, rhs, rcond)

        self.fit_info['residuals'] = resids
        self.fit_info['rank'] = rank
        self.fit_info['singular_values'] = sval

        lacoef = (lacoef.T / scl).T
        self.fit_info['params'] = lacoef

        # TODO: Only Polynomial models currently have an _order attribute;
        # maybe change this to read isinstance(model, PolynomialBase)
        if hasattr(model_copy, '_order') and rank != model_copy._order:
            warnings.warn("The fit may be poorly conditioned\n",
                          AstropyUserWarning)

        _fitter_to_model_params(model_copy, lacoef.flatten())
        return model_copy


class FittingWithOutlierRemoval(object):
    """
    This class combines an outlier removal technique with a fitting procedure.
    Basically, given a number of iterations ``niter``, outliers are removed
    and fitting is performed for each iteration.

    Parameters
    ----------
    fitter : An Astropy fitter
        An instance of any Astropy fitter, i.e., LinearLSQFitter,
        LevMarLSQFitter, SLSQPLSQFitter, SimplexLSQFitter, JointFitter.
    outlier_func : function
        A function for outlier removal.
    niter : int (optional)
        Number of iterations.
    outlier_kwargs : dict (optional)
        Keyword arguments for outlier_func.
    """

    def __init__(self, fitter, outlier_func, niter=3, **outlier_kwargs):
        self.fitter = fitter
        self.outlier_func = outlier_func
        self.niter = niter
        self.outlier_kwargs = outlier_kwargs

    def __str__(self):
        return ("Fitter: {0}\nOutlier function: {1}\nNum. of iterations: {2}" +
                ("\nOutlier func. args.: {3}"))\
                .format(self.fitter__class__.__name__,
                        self.outlier_func.__name__, self.niter,
                        self.outlier_kwargs)

    def __repr__(self):
        return ("{0}(fitter: {1}, outlier_func: {2}," +
                " niter: {3}, outlier_kwargs: {4})")\
                 .format(self.__class__.__name__,
                         self.fitter.__class__.__name__,
                         self.outlier_func.__name__, self.niter,
                         self.outlier_kwargs)

    def __call__(self, model, x, y, z=None, weights=None, **kwargs):
        """
        Parameters
        ----------
        model : `~astropy.modeling.FittableModel`
            An analytic model which will be fit to the provided data.
            This also contains the initial guess for an optimization
            algorithm.
        x : array-like
            Input coordinates.
        y : array-like
            Data measurements (1D case) or input coordinates (2D case).
        z : array-like (optional)
            Data measurements (2D case).
        weights : array-like (optional)
            Weights to be passed to the fitter.
        kwargs : dict (optional)
            Keyword arguments to be passed to the fitter.

        Returns
        -------
        filtered_data : numpy.ma.core.MaskedArray
            Data used to perform the fitting after outlier removal.
        fitted_model : `~astropy.modeling.FittableModel`
            Fitted model after outlier removal.
        """

        fitted_model = self.fitter(model, x, y, z, weights=weights, **kwargs)
        filtered_weights = weights
        if z is None:
            filtered_data = y
            for n in range(self.niter):
                filtered_data = self.outlier_func(filtered_data - fitted_model(x),
                                                  **self.outlier_kwargs)
                filtered_data += fitted_model(x)
                if weights is not None:
                    filtered_weights = weights[~filtered_data.mask]
                fitted_model = self.fitter(fitted_model,
                               x[~filtered_data.mask],
                               filtered_data.data[~filtered_data.mask],
                               weights=filtered_weights,
                               **kwargs)
        else:
            filtered_data = z
            for n in range(self.niter):
                filtered_data = self.outlier_func(filtered_data - fitted_model(x, y),
                                                  **self.outlier_kwargs)
                filtered_data += fitted_model(x, y)
                if weights is not None:
                    filtered_weights = weights[~filtered_data.mask]
                fitted_model = self.fitter(fitted_model,
                               x[~filtered_data.mask],
                               y[~filtered_data.mask],
                               filtered_data.data[~filtered_data.mask],
                               weights=filtered_weights,
                               **kwargs)
        return filtered_data, fitted_model


@six.add_metaclass(_FitterMeta)
class LevMarLSQFitter(object):
    """
    Levenberg-Marquardt algorithm and least squares statistic.

    Attributes
    ----------
    fit_info : dict
        The `scipy.optimize.leastsq` result for the most recent fit (see
        notes).

    Notes
    -----
    The ``fit_info`` dictionary contains the values returned by
    `scipy.optimize.leastsq` for the most recent fit, including the values from
    the ``infodict`` dictionary it returns. See the `scipy.optimize.leastsq`
    documentation for details on the meaning of these values. Note that the
    ``x`` return value is *not* included (as it is instead the parameter values
    of the returned model).

    Additionally, one additional element of ``fit_info`` is computed whenever a
    model is fit, with the key 'param_cov'. The corresponding value is the
    covariance matrix of the parameters as a 2D numpy array.  The order of the
    matrix elements matches the order of the parameters in the fitted model
    (i.e., the same order as ``model.param_names``).
    """

    supported_constraints = ['fixed', 'tied', 'bounds']
    """
    The constraint types supported by this fitter type.
    """

    def __init__(self):
        self.fit_info = {'nfev': None,
                         'fvec': None,
                         'fjac': None,
                         'ipvt': None,
                         'qtf': None,
                         'message': None,
                         'ierr': None,
                         'param_jac': None,
                         'param_cov': None}

        super(LevMarLSQFitter, self).__init__()

    def objective_function(self, fps, *args):
        """
        Function to minimize.

        Parameters
        ----------
        fps : list
            parameters returned by the fitter
        args : list
            [model, [weights], [input coordinates]]
        """

        model = args[0]
        weights = args[1]
        _fitter_to_model_params(model, fps)
        meas = args[-1]
        if weights is None:
            return np.ravel(model(*args[2: -1]) - meas)
        else:
            return np.ravel(weights * (model(*args[2: -1]) - meas))

    @fitter_unit_support
    def __call__(self, model, x, y, z=None, weights=None,
                 maxiter=DEFAULT_MAXITER, acc=DEFAULT_ACC,
                 epsilon=DEFAULT_EPS, estimate_jacobian=False):
        """
        Fit data to this model.

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
        weights : array (optional)
            Weights for fitting.
            For data with Gaussian uncertainties, the weights should be
            1/sigma.
        maxiter : int
            maximum number of iterations
        acc : float
            Relative error desired in the approximate solution
        epsilon : float
            A suitable step length for the forward-difference
            approximation of the Jacobian (if model.fjac=None). If
            epsfcn is less than the machine precision, it is
            assumed that the relative errors in the functions are
            of the order of the machine precision.
        estimate_jacobian : bool
            If False (default) and if the model has a fit_deriv method,
            it will be used. Otherwise the Jacobian will be estimated.
            If True, the Jacobian will be estimated in any case.
        equivalencies : list or None, optional and keyword-only argument
            List of *additional* equivalencies that are should be applied in
            case x, y and/or z have units. Default is None.

        Returns
        -------
        model_copy : `~astropy.modeling.FittableModel`
            a copy of the input model with parameters set by the fitter
        """

        from scipy import optimize

        model_copy = _validate_model(model, self.supported_constraints)
        farg = (model_copy, weights, ) + _convert_input(x, y, z)

        if model_copy.fit_deriv is None or estimate_jacobian:
            dfunc = None
        else:
            dfunc = self._wrap_deriv
        init_values, _ = _model_to_fit_params(model_copy)
        fitparams, cov_x, dinfo, mess, ierr = optimize.leastsq(
            self.objective_function, init_values, args=farg, Dfun=dfunc,
            col_deriv=model_copy.col_fit_deriv, maxfev=maxiter, epsfcn=epsilon,
            xtol=acc, full_output=True)
        _fitter_to_model_params(model_copy, fitparams)
        self.fit_info.update(dinfo)
        self.fit_info['cov_x'] = cov_x
        self.fit_info['message'] = mess
        self.fit_info['ierr'] = ierr
        if ierr not in [1, 2, 3, 4]:
            warnings.warn("The fit may be unsuccessful; check "
                          "fit_info['message'] for more information.",
                          AstropyUserWarning)

        # now try to compute the true covariance matrix
        if (len(y) > len(init_values)) and cov_x is not None:
            sum_sqrs = np.sum(self.objective_function(fitparams, *farg)**2)
            dof = len(y) - len(init_values)
            self.fit_info['param_cov'] = cov_x * sum_sqrs / dof
        else:
            self.fit_info['param_cov'] = None

        return model_copy

    @staticmethod
    def _wrap_deriv(params, model, weights, x, y, z=None):
        """
        Wraps the method calculating the Jacobian of the function to account
        for model constraints.

        `scipy.optimize.leastsq` expects the function derivative to have the
        above signature (parlist, (argtuple)). In order to accommodate model
        constraints, instead of using p directly, we set the parameter list in
        this function.
        """

        if weights is None:
            weights = 1.0

        if any(model.fixed.values()) or any(model.tied.values()):
            # update the parameters with the current values from the fitter
            _fitter_to_model_params(model, params)
            if z is None:
                full = np.array(model.fit_deriv(x, *model.parameters))
                if not model.col_fit_deriv:
                    full_deriv = np.ravel(weights) * full.T
                else:
                    full_deriv = np.ravel(weights) * full
            else:
                full = np.array([np.ravel(_) for _ in model.fit_deriv(x, y, *model.parameters)])
                if not model.col_fit_deriv:
                    full_deriv = np.ravel(weights) * full.T
                else:
                    full_deriv = np.ravel(weights) * full

            pars = [getattr(model, name) for name in model.param_names]
            fixed = [par.fixed for par in pars]
            tied = [par.tied for par in pars]
            tied = list(np.where([par.tied is not False for par in pars],
                                 True, tied))
            fix_and_tie = np.logical_or(fixed, tied)
            ind = np.logical_not(fix_and_tie)

            if not model.col_fit_deriv:
                residues = np.asarray(full_deriv[np.nonzero(ind)]).T
            else:
                residues = full_deriv[np.nonzero(ind)]

            return [np.ravel(_) for _ in residues]
        else:
            if z is None:
                return [np.ravel(_) for _ in np.ravel(weights) * np.array(model.fit_deriv(x, *params))]
            else:
                if not model.col_fit_deriv:
                    return [np.ravel(_) for _ in (
                        np.ravel(weights) * np.array(model.fit_deriv(x, y, *params)).T).T]
                else:
                    return [np.ravel(_) for _ in (weights * np.array(model.fit_deriv(x, y, *params)))]


class SLSQPLSQFitter(Fitter):
    """
    SLSQP optimization algorithm and least squares statistic.


    Raises
    ------
    ModelLinearityError
        A linear model is passed to a nonlinear fitter

    """

    supported_constraints = SLSQP.supported_constraints

    def __init__(self):
        super(SLSQPLSQFitter, self).__init__(optimizer=SLSQP, statistic=leastsquare)
        self.fit_info = {}

    @fitter_unit_support
    def __call__(self, model, x, y, z=None, weights=None, **kwargs):
        """
        Fit data to this model.

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
        weights : array (optional)
            Weights for fitting.
            For data with Gaussian uncertainties, the weights should be
            1/sigma.
        kwargs : dict
            optional keyword arguments to be passed to the optimizer or the statistic

        verblevel : int
            0-silent
            1-print summary upon completion,
            2-print summary after each iteration
        maxiter : int
            maximum number of iterations
        epsilon : float
            the step size for finite-difference derivative estimates
        acc : float
            Requested accuracy
        equivalencies : list or None, optional and keyword-only argument
            List of *additional* equivalencies that are should be applied in
            case x, y and/or z have units. Default is None.

        Returns
        -------
        model_copy : `~astropy.modeling.FittableModel`
            a copy of the input model with parameters set by the fitter
        """

        model_copy = _validate_model(model, self._opt_method.supported_constraints)
        farg = _convert_input(x, y, z)
        farg = (model_copy, weights, ) + farg
        p0, _ = _model_to_fit_params(model_copy)
        fitparams, self.fit_info = self._opt_method(
            self.objective_function, p0, farg, **kwargs)
        _fitter_to_model_params(model_copy, fitparams)

        return model_copy


class SimplexLSQFitter(Fitter):
    """

    Simplex algorithm and least squares statistic.

    Raises
    ------
    ModelLinearityError
        A linear model is passed to a nonlinear fitter

    """

    supported_constraints = Simplex.supported_constraints

    def __init__(self):
        super(SimplexLSQFitter, self).__init__(optimizer=Simplex,
                                               statistic=leastsquare)
        self.fit_info = {}

    @fitter_unit_support
    def __call__(self, model, x, y, z=None, weights=None, **kwargs):
        """
        Fit data to this model.

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
        weights : array (optional)
            Weights for fitting.
            For data with Gaussian uncertainties, the weights should be
            1/sigma.
        kwargs : dict
            optional keyword arguments to be passed to the optimizer or the statistic

        maxiter : int
            maximum number of iterations
        acc : float
            Relative error in approximate solution
        equivalencies : list or None, optional and keyword-only argument
            List of *additional* equivalencies that are should be applied in
            case x, y and/or z have units. Default is None.

        Returns
        -------
        model_copy : `~astropy.modeling.FittableModel`
            a copy of the input model with parameters set by the fitter
        """

        model_copy = _validate_model(model,
                                     self._opt_method.supported_constraints)
        farg = _convert_input(x, y, z)
        farg = (model_copy, weights, ) + farg

        p0, _ = _model_to_fit_params(model_copy)

        fitparams, self.fit_info = self._opt_method(
            self.objective_function, p0, farg, **kwargs)
        _fitter_to_model_params(model_copy, fitparams)
        return model_copy


@six.add_metaclass(_FitterMeta)
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
            joint_params = self.jointparams[model]
            param_metrics = model._param_metrics
            for param_name in joint_params:
                slice_ = param_metrics[param_name]['slice']
                del params[slice_]
            fparams.extend(params)
        return fparams

    def objective_function(self, fps, *args):
        """
        Function to minimize.

        Parameters
        ----------
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
            joint_params = self.jointparams[model]
            margs = lstsqargs[:model.n_inputs + 1]
            del lstsqargs[:model.n_inputs + 1]
            # separate each model separately fitted parameters
            numfp = len(model._parameters) - len(joint_params)
            mfparams = fitparams[:numfp]

            del fitparams[:numfp]
            # recreate the model parameters
            mparams = []
            param_metrics = model._param_metrics
            for param_name in model.param_names:
                if param_name in joint_params:
                    index = joint_params.index(param_name)
                    # should do this with slices in case the
                    # parameter is not a number
                    mparams.extend([jointfitparams[index]])
                else:
                    slice_ = param_metrics[param_name]['slice']
                    plen = slice_.stop - slice_.start
                    mparams.extend(mfparams[:plen])
                    del mfparams[:plen]
            modelfit = model.evaluate(margs[:-1], *mparams)
            fitted.extend(modelfit - margs[-1])
        return np.ravel(fitted)

    def _verify_input(self):
        if len(self.models) <= 1:
            raise TypeError("Expected >1 models, {} is given".format(
                    len(self.models)))
        if len(self.jointparams.keys()) < 2:
            raise TypeError("At least two parameters are expected, "
                            "{} is given".format(len(self.jointparams.keys())))
        for j in self.jointparams.keys():
            if len(self.jointparams[j]) != len(self.initvals):
                raise TypeError("{} parameter(s) provided but {} expected".format(
                        len(self.jointparams[j]), len(self.initvals)))

    def __call__(self, *args):
        """
        Fit data to these models keeping some of the parameters common to the
        two models.
        """

        from scipy import optimize

        if len(args) != reduce(lambda x, y: x + 1 + y + 1, self.modeldims):
            raise ValueError("Expected {} coordinates in args but {} provided"
                             .format(reduce(lambda x, y: x + 1 + y + 1,
                                            self.modeldims), len(args)))

        self.fitparams[:], _ = optimize.leastsq(self.objective_function,
                                                self.fitparams, args=args)

        fparams = self.fitparams[:]
        numjp = len(self.initvals)
        # make a separate list of the joint fitted parameters
        jointfitparams = fparams[:numjp]
        del fparams[:numjp]

        for model in self.models:
            # extract each model's fitted parameters
            joint_params = self.jointparams[model]
            numfp = len(model._parameters) - len(joint_params)
            mfparams = fparams[:numfp]

            del fparams[:numfp]
            # recreate the model parameters
            mparams = []
            param_metrics = model._param_metrics
            for param_name in model.param_names:
                if param_name in joint_params:
                    index = joint_params.index(param_name)
                    # should do this with slices in case the parameter
                    # is not a number
                    mparams.extend([jointfitparams[index]])
                else:
                    slice_ = param_metrics[param_name]['slice']
                    plen = slice_.stop - slice_.start
                    mparams.extend(mfparams[:plen])
                    del mfparams[:plen]
            model.parameters = np.array(mparams)


def _convert_input(x, y, z=None, n_models=1, model_set_axis=0):
    """Convert inputs to float arrays."""

    x = np.asarray(x, dtype=np.float)
    y = np.asarray(y, dtype=np.float)

    if z is not None:
        z = np.asarray(z, dtype=np.float)

    # For compatibility with how the linear fitter code currently expects to
    # work, shift the dependent variable's axes to the expected locations
    if n_models > 1:
        if z is None:
            if y.shape[model_set_axis] != n_models:
                raise ValueError(
                    "Number of data sets (y array is expected to equal "
                    "the number of parameter sets)")
            # For a 1-D model the y coordinate's model-set-axis is expected to
            # be last, so that its first dimension is the same length as the x
            # coordinates.  This is in line with the expectations of
            # numpy.linalg.lstsq:
            # http://docs.scipy.org/doc/numpy/reference/generated/numpy.linalg.lstsq.html
            # That is, each model should be represented by a column.  TODO:
            # Obviously this is a detail of np.linalg.lstsq and should be
            # handled specifically by any fitters that use it...
            y = np.rollaxis(y, model_set_axis, y.ndim)
        else:
            # Shape of z excluding model_set_axis
            z_shape = z.shape[:model_set_axis] + z.shape[model_set_axis + 1:]

            if not (x.shape == y.shape == z_shape):
                raise ValueError("x, y and z should have the same shape")

    if z is None:
        farg = (x, y)
    else:
        farg = (x, y, z)
    return farg


# TODO: These utility functions are really particular to handling
# bounds/tied/fixed constraints for scipy.optimize optimizers that do not
# support them inherently; this needs to be reworked to be clear about this
# distinction (and the fact that these are not necessarily applicable to any
# arbitrary fitter--as evidenced for example by the fact that JointFitter has
# its own versions of these)
# TODO: Most of this code should be entirely rewritten; it should not be as
# inefficient as it is.
def _fitter_to_model_params(model, fps):
    """
    Constructs the full list of model parameters from the fitted and
    constrained parameters.
    """

    _, fit_param_indices = _model_to_fit_params(model)

    has_tied = any(model.tied.values())
    has_fixed = any(model.fixed.values())
    has_bound = any(b != (None, None) for b in model.bounds.values())

    if not (has_tied or has_fixed or has_bound):
        # We can just assign directly
        model.parameters = fps
        return

    fit_param_indices = set(fit_param_indices)
    offset = 0
    param_metrics = model._param_metrics
    for idx, name in enumerate(model.param_names):
        if idx not in fit_param_indices:
            continue

        slice_ = param_metrics[name]['slice']
        shape = param_metrics[name]['shape']
        # This is determining which range of fps (the fitted parameters) maps
        # to parameters of the model
        size = reduce(operator.mul, shape, 1)

        values = fps[offset:offset + size]

        # Check bounds constraints
        if model.bounds[name] != (None, None):
            _min, _max = model.bounds[name]
            if _min is not None:
                values = np.fmax(values, _min)
            if _max is not None:
                values = np.fmin(values, _max)

        model.parameters[slice_] = values
        offset += size

    # This has to be done in a separate loop due to how tied parameters are
    # currently evaluated (the fitted parameters need to actually be *set* on
    # the model first, for use in evaluating the "tied" expression--it might be
    # better to change this at some point
    if has_tied:
        for idx, name in enumerate(model.param_names):
            if model.tied[name]:
                value = model.tied[name](model)
                slice_ = param_metrics[name]['slice']
                model.parameters[slice_] = value


def _model_to_fit_params(model):
    """
    Convert a model instance's parameter array to an array that can be used
    with a fitter that doesn't natively support fixed or tied parameters.
    In particular, it removes fixed/tied parameters from the parameter
    array.

    These may be a subset of the model parameters, if some of them are held
    constant or tied.
    """

    fitparam_indices = list(range(len(model.param_names)))
    if any(model.fixed.values()) or any(model.tied.values()):
        params = list(model.parameters)
        param_metrics = model._param_metrics
        for idx, name in list(enumerate(model.param_names))[::-1]:
            if model.fixed[name] or model.tied[name]:
                slice_ = param_metrics[name]['slice']
                del params[slice_]
                del fitparam_indices[idx]
        return (np.array(params), fitparam_indices)
    else:
        return (model.parameters, fitparam_indices)


def _validate_constraints(supported_constraints, model):
    """Make sure model constraints are supported by the current fitter."""

    message = 'Optimizer cannot handle {0} constraints.'

    if (any(six.itervalues(model.fixed)) and
            'fixed' not in supported_constraints):
        raise UnsupportedConstraintError(
                message.format('fixed parameter'))

    if any(six.itervalues(model.tied)) and 'tied' not in supported_constraints:
        raise UnsupportedConstraintError(
                message.format('tied parameter'))

    if (any(tuple(b) != (None, None) for b in six.itervalues(model.bounds)) and
            'bounds' not in supported_constraints):
        raise UnsupportedConstraintError(
                message.format('bound parameter'))

    if model.eqcons and 'eqcons' not in supported_constraints:
        raise UnsupportedConstraintError(message.format('equality'))

    if model.ineqcons and 'ineqcons' not in supported_constraints:
        raise UnsupportedConstraintError(message.format('inequality'))


def _validate_model(model, supported_constraints):
    """
    Check that model and fitter are compatible and return a copy of the model.
    """

    if not model.fittable:
        raise ValueError("Model does not appear to be fittable.")
    if model.linear:
        warnings.warn('Model is linear in parameters; '
                      'consider using linear fitting methods.',
                      AstropyUserWarning)
    elif len(model) != 1:
        # for now only single data sets ca be fitted
        raise ValueError("Non-linear fitters can only fit "
                         "one data set at a time.")
    _validate_constraints(supported_constraints, model)

    model_copy = model.copy()
    return model_copy


def populate_entry_points(entry_points):
    """
    This injects entry points into the `astropy.modeling.fitting` namespace.
    This provides a means of inserting a fitting routine without requirement
    of it being merged into astropy's core.

    Parameters
    ----------

    entry_points : a list of `~pkg_resources.EntryPoint`
                  entry_points are objects which encapsulate
                  importable objects and are defined on the
                  installation of a package.
    Notes
    -----
    An explanation of entry points can be found `here <http://setuptools.readthedocs.io/en/latest/setuptools.html#dynamic-discovery-of-services-and-plugins>`

    """

    for entry_point in entry_points:
        name = entry_point.name
        try:
            entry_point = entry_point.load()
        except Exception as e:
            # This stops the fitting from choking if an entry_point produces an error.
            warnings.warn(AstropyUserWarning('{type} error occurred in entry '
                                             'point {name}.' .format(type=type(e).__name__, name=name)))
        else:
            if not inspect.isclass(entry_point):
                warnings.warn(AstropyUserWarning(
                    'Modeling entry point {0} expected to be a '
                    'Class.' .format(name)))
            else:
                if issubclass(entry_point, Fitter):
                    name = entry_point.__name__
                    globals()[name] = entry_point
                    __all__.append(name)
                else:
                    warnings.warn(AstropyUserWarning(
                        'Modeling entry point {0} expected to extend '
                        'astropy.modeling.Fitter' .format(name)))


# this is so fitting doesn't choke if pkg_resources doesn't exist
if HAS_PKG:
    populate_entry_points(iter_entry_points(group='astropy.modeling', name=None))
