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
# pylint: disable=invalid-name

import abc
import inspect
import operator
import warnings
from functools import reduce, wraps
from importlib.metadata import entry_points

import numpy as np

from astropy.units import Quantity
from astropy.utils.decorators import deprecated
from astropy.utils.exceptions import AstropyUserWarning

from .optimizers import DEFAULT_ACC, DEFAULT_EPS, DEFAULT_MAXITER, SLSQP, Simplex
from .spline import (
    SplineExactKnotsFitter,
    SplineInterpolateFitter,
    SplineSmoothingFitter,
    SplineSplrepFitter,
)
from .statistic import leastsquare
from .utils import _combine_equivalency_dict, poly_map_domain

__all__ = [
    "LinearLSQFitter",
    "LevMarLSQFitter",
    "TRFLSQFitter",
    "DogBoxLSQFitter",
    "LMLSQFitter",
    "FittingWithOutlierRemoval",
    "SLSQPLSQFitter",
    "SimplexLSQFitter",
    "JointFitter",
    "Fitter",
    "ModelLinearityError",
    "ModelsError",
    "SplineExactKnotsFitter",
    "SplineInterpolateFitter",
    "SplineSmoothingFitter",
    "SplineSplrepFitter",
]


# Statistic functions implemented in `astropy.modeling.statistic.py
STATISTICS = [leastsquare]

# Optimizers implemented in `astropy.modeling.optimizers.py
OPTIMIZERS = [Simplex, SLSQP]


class NonFiniteValueError(RuntimeError):
    """
    Error raised when attempting to a non-finite value.
    """


class Covariance:
    """Class for covariance matrix calculated by fitter."""

    def __init__(self, cov_matrix, param_names):
        self.cov_matrix = cov_matrix
        self.param_names = param_names

    def pprint(self, max_lines, round_val):
        # Print and label lower triangle of covariance matrix
        # Print rows for params up to `max_lines`, round floats to 'round_val'
        longest_name = max(len(x) for x in self.param_names)
        ret_str = "parameter variances / covariances \n"
        fstring = f'{"": <{longest_name}}| {{0}}\n'
        for i, row in enumerate(self.cov_matrix):
            if i <= max_lines - 1:
                param = self.param_names[i]
                ret_str += fstring.replace(" " * len(param), param, 1).format(
                    repr(np.round(row[: i + 1], round_val))[7:-2]
                )
            else:
                ret_str += "..."
        return ret_str.rstrip()

    def __repr__(self):
        return self.pprint(max_lines=10, round_val=3)

    def __getitem__(self, params):
        # index covariance matrix by parameter names or indices
        if len(params) != 2:
            raise ValueError("Covariance must be indexed by two values.")
        if all(isinstance(item, str) for item in params):
            i1, i2 = self.param_names.index(params[0]), self.param_names.index(
                params[1]
            )
        elif all(isinstance(item, int) for item in params):
            i1, i2 = params
        else:
            raise TypeError(
                "Covariance can be indexed by two parameter names or integer indices."
            )
        return self.cov_matrix[i1][i2]


class StandardDeviations:
    """Class for fitting uncertainties."""

    def __init__(self, cov_matrix, param_names):
        self.param_names = param_names
        self.stds = self._calc_stds(cov_matrix)

    def _calc_stds(self, cov_matrix):
        # sometimes scipy lstsq returns a non-sensical negative vals in the
        # diagonals of the cov_x it computes.
        stds = [np.sqrt(x) if x > 0 else None for x in np.diag(cov_matrix)]
        return stds

    def pprint(self, max_lines, round_val):
        longest_name = max(len(x) for x in self.param_names)
        ret_str = "standard deviations\n"
        for i, std in enumerate(self.stds):
            if i <= max_lines - 1:
                param = self.param_names[i]
                ret_str += (
                    f"{param}{' ' * (longest_name - len(param))}| "
                    f"{np.round(std, round_val)}\n"
                )
            else:
                ret_str += "..."
        return ret_str.rstrip()

    def __repr__(self):
        return self.pprint(max_lines=10, round_val=3)

    def __getitem__(self, param):
        if isinstance(param, str):
            i = self.param_names.index(param)
        elif isinstance(param, int):
            i = param
        else:
            raise TypeError(
                "Standard deviation can be indexed by parameter name or integer."
            )
        return self.stds[i]


class ModelsError(Exception):
    """Base class for model exceptions."""


class ModelLinearityError(ModelsError):
    """Raised when a non-linear model is passed to a linear fitter."""


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
        cls = super().__new__(mcls, name, bases, members)

        if not inspect.isabstract(cls) and not name.startswith("_"):
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
        equivalencies = kwargs.pop("equivalencies", None)

        data_has_units = (
            isinstance(x, Quantity)
            or isinstance(y, Quantity)
            or isinstance(z, Quantity)
        )

        model_has_units = model._has_units

        if data_has_units or model_has_units:
            if model._supports_unit_fitting:
                # We now combine any instance-level input equivalencies with user
                # specified ones at call-time.

                input_units_equivalencies = _combine_equivalency_dict(
                    model.inputs, equivalencies, model.input_units_equivalencies
                )

                # If input_units is defined, we transform the input data into those
                # expected by the model. We hard-code the input names 'x', and 'y'
                # here since FittableModel instances have input names ('x',) or
                # ('x', 'y')

                if model.input_units is not None:
                    if isinstance(x, Quantity):
                        x = x.to(
                            model.input_units[model.inputs[0]],
                            equivalencies=input_units_equivalencies[model.inputs[0]],
                        )
                    if isinstance(y, Quantity) and z is not None:
                        y = y.to(
                            model.input_units[model.inputs[1]],
                            equivalencies=input_units_equivalencies[model.inputs[1]],
                        )

                # Create a dictionary mapping the real model inputs and outputs
                # names to the data. This remapping of names must be done here, after
                # the input data is converted to the correct units.
                rename_data = {model.inputs[0]: x}
                if z is not None:
                    rename_data[model.outputs[0]] = z
                    rename_data[model.inputs[1]] = y
                else:
                    rename_data[model.outputs[0]] = y
                    rename_data["z"] = None

                # We now strip away the units from the parameters, taking care to
                # first convert any parameters to the units that correspond to the
                # input units (to make sure that initial guesses on the parameters)
                # are in the right unit system
                model = model.without_units_for_data(**rename_data)
                if isinstance(model, tuple):
                    rename_data["_left_kwargs"] = model[1]
                    rename_data["_right_kwargs"] = model[2]
                    model = model[0]

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
                    if isinstance(z, Quantity):
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
                    model_new = model_new.with_units_from_data(**rename_data)
                return model_new

            else:
                raise NotImplementedError(
                    "This model does not support being fit to data with units."
                )

        else:
            return func(self, model, x, y, z=z, **kwargs)

    return wrapper


class Fitter(metaclass=_FitterMeta):
    """
    Base class for all fitters.

    Parameters
    ----------
    optimizer : callable
        A callable implementing an optimization algorithm
    statistic : callable
        Statistic function

    """

    supported_constraints = []

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
        fitter_to_model_params(model, fps)
        res = self._stat_method(meas, model, *args[1:-1])
        return res

    @staticmethod
    def _add_fitting_uncertainties(*args):
        """
        When available, calculate and sets the parameter covariance matrix
        (model.cov_matrix) and standard deviations (model.stds).
        """
        return None

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
class LinearLSQFitter(metaclass=_FitterMeta):
    """
    A class performing a linear least square fitting.
    Uses `numpy.linalg.lstsq` to do the fitting.
    Given a model and data, fits the model to the data and changes the
    model's parameters. Keeps a dictionary of auxiliary fitting information.

    Notes
    -----
    Note that currently LinearLSQFitter does not support compound models.
    """

    supported_constraints = ["fixed"]
    supports_masked_input = True

    def __init__(self, calc_uncertainties=False):
        self.fit_info = {
            "residuals": None,
            "rank": None,
            "singular_values": None,
            "params": None,
        }
        self._calc_uncertainties = calc_uncertainties

    @staticmethod
    def _is_invertible(m):
        """Check if inverse of matrix can be obtained."""
        if m.shape[0] != m.shape[1]:
            return False
        if np.linalg.matrix_rank(m) < m.shape[0]:
            return False
        return True

    def _add_fitting_uncertainties(self, model, a, n_coeff, x, y, z=None, resids=None):
        """
        Calculate and parameter covariance matrix and standard deviations
        and set `cov_matrix` and `stds` attributes.
        """
        x_dot_x_prime = np.dot(a.T, a)
        masked = False or hasattr(y, "mask")

        # check if invertible. if not, can't calc covariance.
        if not self._is_invertible(x_dot_x_prime):
            return model
        inv_x_dot_x_prime = np.linalg.inv(x_dot_x_prime)

        if z is None:  # 1D models
            if len(model) == 1:  # single model
                mask = None
                if masked:
                    mask = y.mask
                xx = np.ma.array(x, mask=mask)
                RSS = [(1 / (xx.count() - n_coeff)) * resids]

            if len(model) > 1:  # model sets
                RSS = []  # collect sum residuals squared for each model in set
                for j in range(len(model)):
                    mask = None
                    if masked:
                        mask = y.mask[..., j].flatten()
                    xx = np.ma.array(x, mask=mask)
                    eval_y = model(xx, model_set_axis=False)
                    eval_y = np.rollaxis(eval_y, model.model_set_axis)[j]
                    RSS.append(
                        (1 / (xx.count() - n_coeff)) * np.sum((y[..., j] - eval_y) ** 2)
                    )

        else:  # 2D model
            if len(model) == 1:
                mask = None
                if masked:
                    warnings.warn(
                        "Calculation of fitting uncertainties "
                        "for 2D models with masked values not "
                        "currently supported.\n",
                        AstropyUserWarning,
                    )
                    return
                xx, _ = np.ma.array(x, mask=mask), np.ma.array(y, mask=mask)
                # len(xx) instead of xx.count. this will break if values are masked?
                RSS = [(1 / (len(xx) - n_coeff)) * resids]
            else:
                RSS = []
                for j in range(len(model)):
                    eval_z = model(x, y, model_set_axis=False)
                    mask = None  # need to figure out how to deal w/ masking here.
                    if model.model_set_axis == 1:
                        # model_set_axis passed when evaluating only refers to input shapes
                        # so output must be reshaped for model_set_axis=1.
                        eval_z = np.rollaxis(eval_z, 1)
                    eval_z = eval_z[j]
                    RSS.append(
                        [(1 / (len(x) - n_coeff)) * np.sum((z[j] - eval_z) ** 2)]
                    )

        covs = [inv_x_dot_x_prime * r for r in RSS]
        free_param_names = [
            x
            for x in model.fixed
            if (model.fixed[x] is False) and (model.tied[x] is False)
        ]

        if len(covs) == 1:
            model.cov_matrix = Covariance(covs[0], model.param_names)
            model.stds = StandardDeviations(covs[0], free_param_names)
        else:
            model.cov_matrix = [Covariance(cov, model.param_names) for cov in covs]
            model.stds = [StandardDeviations(cov, free_param_names) for cov in covs]

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
            if hasattr(model, "domain") and model.domain is None:
                model.domain = [x.min(), x.max()]
            if hasattr(model, "window") and model.window is None:
                model.window = [-1, 1]
            return poly_map_domain(x, model.domain, model.window)
        else:
            if hasattr(model, "x_domain") and model.x_domain is None:
                model.x_domain = [x.min(), x.max()]
            if hasattr(model, "y_domain") and model.y_domain is None:
                model.y_domain = [y.min(), y.max()]
            if hasattr(model, "x_window") and model.x_window is None:
                model.x_window = [-1.0, 1.0]
            if hasattr(model, "y_window") and model.y_window is None:
                model.y_window = [-1.0, 1.0]

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
            Input coordinates
        y : array-like
            Input coordinates
        z : array-like, optional
            Input coordinates.
            If the dependent (``y`` or ``z``) coordinate values are provided
            as a `numpy.ma.MaskedArray`, any masked points are ignored when
            fitting. Note that model set fitting is significantly slower when
            there are masked points (not just an empty mask), as the matrix
            equation has to be solved for each model separately when their
            coordinate grids differ.
        weights : array, optional
            Weights for fitting.
            For data with Gaussian uncertainties, the weights should be
            1/sigma.
        rcond :  float, optional
            Cut-off ratio for small singular values of ``a``.
            Singular values are set to zero if they are smaller than ``rcond``
            times the largest singular value of ``a``.
        equivalencies : list or None, optional, keyword-only
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
            raise ModelLinearityError(
                "Model is not linear in parameters, "
                "linear fit methods should not be used."
            )

        if hasattr(model, "submodel_names"):
            raise ValueError("Model must be simple, not compound")

        _validate_constraints(self.supported_constraints, model)

        model_copy = model.copy()
        model_copy.sync_constraints = False
        _, fitparam_indices, _ = model_to_fit_params(model_copy)

        if model_copy.n_inputs == 2 and z is None:
            raise ValueError("Expected x, y and z for a 2 dimensional model.")

        farg = _convert_input(
            x, y, z, n_models=len(model_copy), model_set_axis=model_copy.model_set_axis
        )

        n_fixed = sum(model_copy.fixed.values())

        # This is also done by _convert_inputs, but we need it here to allow
        # checking the array dimensionality before that gets called:
        if weights is not None:
            weights = np.asarray(weights, dtype=float)

        if n_fixed:
            # The list of fixed params is the complement of those being fitted:
            fixparam_indices = [
                idx
                for idx in range(len(model_copy.param_names))
                if idx not in fitparam_indices
            ]

            # Construct matrix of user-fixed parameters that can be dotted with
            # the corresponding fit_deriv() terms, to evaluate corrections to
            # the dependent variable in order to fit only the remaining terms:
            fixparams = np.asarray(
                [
                    getattr(model_copy, model_copy.param_names[idx]).value
                    for idx in fixparam_indices
                ]
            )

        if len(farg) == 2:
            x, y = farg

            if weights is not None:
                # If we have separate weights for each model, apply the same
                # conversion as for the data, otherwise check common weights
                # as if for a single model:
                _, weights = _convert_input(
                    x,
                    weights,
                    n_models=len(model_copy) if weights.ndim == y.ndim else 1,
                    model_set_axis=model_copy.model_set_axis,
                )

            # map domain into window
            if hasattr(model_copy, "domain"):
                x = self._map_domain_window(model_copy, x)
            if n_fixed:
                lhs = np.asarray(
                    self._deriv_with_constraints(model_copy, fitparam_indices, x=x)
                )
                fixderivs = self._deriv_with_constraints(
                    model_copy, fixparam_indices, x=x
                )
            else:
                lhs = np.asarray(model_copy.fit_deriv(x, *model_copy.parameters))
            sum_of_implicit_terms = model_copy.sum_of_implicit_terms(x)
            rhs = y
        else:
            x, y, z = farg

            if weights is not None:
                # If we have separate weights for each model, apply the same
                # conversion as for the data, otherwise check common weights
                # as if for a single model:
                _, _, weights = _convert_input(
                    x,
                    y,
                    weights,
                    n_models=len(model_copy) if weights.ndim == z.ndim else 1,
                    model_set_axis=model_copy.model_set_axis,
                )

            # map domain into window
            if hasattr(model_copy, "x_domain"):
                x, y = self._map_domain_window(model_copy, x, y)

            if n_fixed:
                lhs = np.asarray(
                    self._deriv_with_constraints(model_copy, fitparam_indices, x=x, y=y)
                )
                fixderivs = self._deriv_with_constraints(
                    model_copy, fixparam_indices, x=x, y=y
                )
            else:
                lhs = np.asanyarray(model_copy.fit_deriv(x, y, *model_copy.parameters))
            sum_of_implicit_terms = model_copy.sum_of_implicit_terms(x, y)

            if len(model_copy) > 1:
                # Just to be explicit (rather than baking in False == 0):
                model_axis = model_copy.model_set_axis or 0

                if z.ndim > 2:
                    # For higher-dimensional z, flatten all the axes except the
                    # dimension along which models are stacked and transpose so
                    # the model axis is *last* (I think this resolves Erik's
                    # pending generalization from 80a6f25a):
                    rhs = np.rollaxis(z, model_axis, z.ndim)
                    rhs = rhs.reshape(-1, rhs.shape[-1])
                else:
                    # This "else" seems to handle the corner case where the
                    # user has already flattened x/y before attempting a 2D fit
                    # but z has a second axis for the model set. NB. This is
                    # ~5-10x faster than using rollaxis.
                    rhs = z.T if model_axis == 0 else z

                if weights is not None:
                    # Same for weights
                    if weights.ndim > 2:
                        # Separate 2D weights for each model:
                        weights = np.rollaxis(weights, model_axis, weights.ndim)
                        weights = weights.reshape(-1, weights.shape[-1])
                    elif weights.ndim == z.ndim:
                        # Separate, flattened weights for each model:
                        weights = weights.T if model_axis == 0 else weights
                    else:
                        # Common weights for all the models:
                        weights = weights.flatten()
            else:
                rhs = z.flatten()
                if weights is not None:
                    weights = weights.flatten()

        # If the derivative is defined along rows (as with non-linear models)
        if model_copy.col_fit_deriv:
            lhs = np.asarray(lhs).T

        # Some models (eg. Polynomial1D) don't flatten multi-dimensional inputs
        # when constructing their Vandermonde matrix, which can lead to obscure
        # failures below. Ultimately, np.linalg.lstsq can't handle >2D matrices,
        # so just raise a slightly more informative error when this happens:
        if np.asanyarray(lhs).ndim > 2:
            raise ValueError(
                f"{type(model_copy).__name__} gives unsupported >2D "
                "derivative matrix for this x/y"
            )

        # Subtract any terms fixed by the user from (a copy of) the RHS, in
        # order to fit the remaining terms correctly:
        if n_fixed:
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
            # same input coordinates, as is currently the case.
            if len(model_copy) > 1:
                sum_of_implicit_terms = sum_of_implicit_terms[..., np.newaxis]
            rhs = rhs - sum_of_implicit_terms

        if weights is not None:
            if rhs.ndim == 2:
                if weights.shape == rhs.shape:
                    # separate weights for multiple models case: broadcast
                    # lhs to have more dimension (for each model)
                    lhs = lhs[..., np.newaxis] * weights[:, np.newaxis]
                    rhs = rhs * weights
                else:
                    lhs *= weights[:, np.newaxis]
                    # Don't modify in-place in case rhs was the original
                    # dependent variable array
                    rhs = rhs * weights[:, np.newaxis]
            else:
                lhs *= weights[:, np.newaxis]
                rhs = rhs * weights

        scl = (lhs * lhs).sum(0)
        lhs /= scl

        masked = np.any(np.ma.getmask(rhs))
        if weights is not None and not masked and np.any(np.isnan(lhs)):
            raise ValueError(
                "Found NaNs in the coefficient matrix, which "
                "should not happen and would crash the lapack "
                "routine. Maybe check that weights are not null."
            )

        a = None  # need for calculating covarience

        if (masked and len(model_copy) > 1) or (
            weights is not None and weights.ndim > 1
        ):
            # Separate masks or weights for multiple models case: Numpy's
            # lstsq supports multiple dimensions only for rhs, so we need to
            # loop manually on the models. This may be fixed in the future
            # with https://github.com/numpy/numpy/pull/15777.

            # Initialize empty array of coefficients and populate it one model
            # at a time. The shape matches the number of coefficients from the
            # Vandermonde matrix and the number of models from the RHS:
            lacoef = np.zeros(lhs.shape[1:2] + rhs.shape[-1:], dtype=rhs.dtype)

            # Arrange the lhs as a stack of 2D matrices that we can iterate
            # over to get the correctly-orientated lhs for each model:
            if lhs.ndim > 2:
                lhs_stack = np.rollaxis(lhs, -1, 0)
            else:
                lhs_stack = np.broadcast_to(lhs, rhs.shape[-1:] + lhs.shape)

            # Loop over the models and solve for each one. By this point, the
            # model set axis is the second of two. Transpose rather than using,
            # say, np.moveaxis(array, -1, 0), since it's slightly faster and
            # lstsq can't handle >2D arrays anyway. This could perhaps be
            # optimized by collecting together models with identical masks
            # (eg. those with no rejected points) into one operation, though it
            # will still be relatively slow when calling lstsq repeatedly.
            for model_lhs, model_rhs, model_lacoef in zip(lhs_stack, rhs.T, lacoef.T):
                # Cull masked points on both sides of the matrix equation:
                good = ~model_rhs.mask if masked else slice(None)
                model_lhs = model_lhs[good]
                model_rhs = model_rhs[good][..., np.newaxis]
                a = model_lhs

                # Solve for this model:
                t_coef, resids, rank, sval = np.linalg.lstsq(
                    model_lhs, model_rhs, rcond
                )
                model_lacoef[:] = t_coef.T

        else:
            # If we're fitting one or more models over a common set of points,
            # we only have to solve a single matrix equation, which is an order
            # of magnitude faster than calling lstsq() once per model below:

            good = ~rhs.mask if masked else slice(None)  # latter is a no-op
            a = lhs[good]
            # Solve for one or more models:
            lacoef, resids, rank, sval = np.linalg.lstsq(lhs[good], rhs[good], rcond)

        self.fit_info["residuals"] = resids
        self.fit_info["rank"] = rank
        self.fit_info["singular_values"] = sval

        lacoef /= scl[:, np.newaxis] if scl.ndim < rhs.ndim else scl
        self.fit_info["params"] = lacoef

        fitter_to_model_params(model_copy, lacoef.flatten())

        # TODO: Only Polynomial models currently have an _order attribute;
        # maybe change this to read isinstance(model, PolynomialBase)
        if (
            hasattr(model_copy, "_order")
            and len(model_copy) == 1
            and rank < (model_copy._order - n_fixed)
        ):
            warnings.warn("The fit may be poorly conditioned\n", AstropyUserWarning)

        # calculate and set covariance matrix and standard devs. on model
        if self._calc_uncertainties:
            if len(y) > len(lacoef):
                self._add_fitting_uncertainties(
                    model_copy, a * scl, len(lacoef), x, y, z, resids
                )
        model_copy.sync_constraints = True
        return model_copy


class FittingWithOutlierRemoval:
    """
    This class combines an outlier removal technique with a fitting procedure.
    Basically, given a maximum number of iterations ``niter``, outliers are
    removed and fitting is performed for each iteration, until no new outliers
    are found or ``niter`` is reached.

    Parameters
    ----------
    fitter : `Fitter`
        An instance of any Astropy fitter, i.e., LinearLSQFitter,
        LevMarLSQFitter, SLSQPLSQFitter, SimplexLSQFitter, JointFitter. For
        model set fitting, this must understand masked input data (as
        indicated by the fitter class attribute ``supports_masked_input``).
    outlier_func : callable
        A function for outlier removal.
        If this accepts an ``axis`` parameter like the `numpy` functions, the
        appropriate value will be supplied automatically when fitting model
        sets (unless overridden in ``outlier_kwargs``), to find outliers for
        each model separately; otherwise, the same filtering must be performed
        in a loop over models, which is almost an order of magnitude slower.
    niter : int, optional
        Maximum number of iterations.
    outlier_kwargs : dict, optional
        Keyword arguments for outlier_func.

    Attributes
    ----------
    fit_info : dict
        The ``fit_info`` (if any) from the last iteration of the wrapped
        ``fitter`` during the most recent fit. An entry is also added with the
        keyword ``niter`` that records the actual number of fitting iterations
        performed (as opposed to the user-specified maximum).
    """

    def __init__(self, fitter, outlier_func, niter=3, **outlier_kwargs):
        self.fitter = fitter
        self.outlier_func = outlier_func
        self.niter = niter
        self.outlier_kwargs = outlier_kwargs
        self.fit_info = {"niter": None}

    def __str__(self):
        return (
            f"Fitter: {self.fitter.__class__.__name__}\n"
            f"Outlier function: {self.outlier_func.__name__}\n"
            f"Num. of iterations: {self.niter}\n"
            f"Outlier func. args.: {self.outlier_kwargs}"
        )

    def __repr__(self):
        return (
            f"{self.__class__.__name__}(fitter: {self.fitter.__class__.__name__}, "
            f"outlier_func: {self.outlier_func.__name__},"
            f" niter: {self.niter}, outlier_kwargs: {self.outlier_kwargs})"
        )

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
        z : array-like, optional
            Data measurements (2D case).
        weights : array-like, optional
            Weights to be passed to the fitter.
        kwargs : dict, optional
            Keyword arguments to be passed to the fitter.

        Returns
        -------
        fitted_model : `~astropy.modeling.FittableModel`
            Fitted model after outlier removal.
        mask : `numpy.ndarray`
            Boolean mask array, identifying which points were used in the final
            fitting iteration (False) and which were found to be outliers or
            were masked in the input (True).
        """
        # For single models, the data get filtered here at each iteration and
        # then passed to the fitter, which is the historical behavior and
        # works even for fitters that don't understand masked arrays. For model
        # sets, the fitter must be able to filter masked data internally,
        # because fitters require a single set of x/y coordinates whereas the
        # eliminated points can vary between models. To avoid this limitation,
        # we could fall back to looping over individual model fits, but it
        # would likely be fiddly and involve even more overhead (and the
        # non-linear fitters don't work with model sets anyway, as of writing).

        if len(model) == 1:
            model_set_axis = None
        else:
            if (
                not hasattr(self.fitter, "supports_masked_input")
                or self.fitter.supports_masked_input is not True
            ):
                raise ValueError(
                    f"{type(self.fitter).__name__} cannot fit model sets with masked "
                    "values"
                )

            # Fitters use their input model's model_set_axis to determine how
            # their input data are stacked:
            model_set_axis = model.model_set_axis
        # Construct input coordinate tuples for fitters & models that are
        # appropriate for the dimensionality being fitted:
        if z is None:
            coords = (x,)
            data = y
        else:
            coords = x, y
            data = z

        # For model sets, construct a numpy-standard "axis" tuple for the
        # outlier function, to treat each model separately (if supported):
        if model_set_axis is not None:
            if model_set_axis < 0:
                model_set_axis += data.ndim

            if "axis" not in self.outlier_kwargs:  # allow user override
                # This also works for False (like model instantiation):
                self.outlier_kwargs["axis"] = tuple(
                    n for n in range(data.ndim) if n != model_set_axis
                )

        loop = False

        # Starting fit, prior to any iteration and masking:
        fitted_model = self.fitter(model, x, y, z, weights=weights, **kwargs)
        filtered_data = np.ma.masked_array(data)
        if filtered_data.mask is np.ma.nomask:
            filtered_data.mask = False
        filtered_weights = weights
        last_n_masked = filtered_data.mask.sum()
        n = 0  # (allow recording no. of iterations when 0)

        # Perform the iterative fitting:
        for n in range(1, self.niter + 1):
            # (Re-)evaluate the last model:
            model_vals = fitted_model(*coords, model_set_axis=False)

            # Determine the outliers:
            if not loop:
                # Pass axis parameter if outlier_func accepts it, otherwise
                # prepare for looping over models:
                try:
                    filtered_data = self.outlier_func(
                        filtered_data - model_vals, **self.outlier_kwargs
                    )
                # If this happens to catch an error with a parameter other
                # than axis, the next attempt will fail accordingly:
                except TypeError:
                    if model_set_axis is None:
                        raise
                    else:
                        self.outlier_kwargs.pop("axis", None)
                        loop = True

                        # Construct MaskedArray to hold filtered values:
                        filtered_data = np.ma.masked_array(
                            filtered_data,
                            dtype=np.result_type(filtered_data, model_vals),
                            copy=True,
                        )
                        # Make sure the mask is an array, not just nomask:
                        if filtered_data.mask is np.ma.nomask:
                            filtered_data.mask = False

                        # Get views transposed appropriately for iteration
                        # over the set (handling data & mask separately due to
                        # NumPy issue #8506):
                        data_T = np.rollaxis(filtered_data, model_set_axis, 0)
                        mask_T = np.rollaxis(filtered_data.mask, model_set_axis, 0)

            if loop:
                model_vals_T = np.rollaxis(model_vals, model_set_axis, 0)
                for row_data, row_mask, row_mod_vals in zip(
                    data_T, mask_T, model_vals_T
                ):
                    masked_residuals = self.outlier_func(
                        row_data - row_mod_vals, **self.outlier_kwargs
                    )
                    row_data.data[:] = masked_residuals.data
                    row_mask[:] = masked_residuals.mask

                # Issue speed warning after the fact, so it only shows up when
                # the TypeError is genuinely due to the axis argument.
                warnings.warn(
                    "outlier_func did not accept axis argument; "
                    "reverted to slow loop over models.",
                    AstropyUserWarning,
                )

            # Recombine newly-masked residuals with model to get masked values:
            filtered_data += model_vals

            # Re-fit the data after filtering, passing masked/unmasked values
            # for single models / sets, respectively:
            if model_set_axis is None:
                good = ~filtered_data.mask

                if weights is not None:
                    filtered_weights = weights[good]

                fitted_model = self.fitter(
                    fitted_model,
                    *(c[good] for c in coords),
                    filtered_data.data[good],
                    weights=filtered_weights,
                    **kwargs,
                )
            else:
                fitted_model = self.fitter(
                    fitted_model,
                    *coords,
                    filtered_data,
                    weights=filtered_weights,
                    **kwargs,
                )

            # Stop iteration if the masked points are no longer changing (with
            # cumulative rejection we only need to compare how many there are):
            this_n_masked = filtered_data.mask.sum()  # (minimal overhead)
            if this_n_masked == last_n_masked:
                break
            last_n_masked = this_n_masked

        self.fit_info = {"niter": n}
        self.fit_info.update(getattr(self.fitter, "fit_info", {}))

        return fitted_model, filtered_data.mask


class _NonLinearLSQFitter(metaclass=_FitterMeta):
    """
    Base class for Non-Linear least-squares fitters.

    Parameters
    ----------
    calc_uncertainties : bool
        If the covarience matrix should be computed and set in the fit_info.
        Default: False
    use_min_max_bounds : bool
        If the set parameter bounds for a model will be enforced each given
        parameter while fitting via a simple min/max condition.
        Default: True
    """

    supported_constraints = ["fixed", "tied", "bounds"]
    """
    The constraint types supported by this fitter type.
    """

    def __init__(self, calc_uncertainties=False, use_min_max_bounds=True):
        self.fit_info = None
        self._calc_uncertainties = calc_uncertainties
        self._use_min_max_bounds = use_min_max_bounds
        super().__init__()

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
        fitter_to_model_params(model, fps, self._use_min_max_bounds)
        meas = args[-1]

        if weights is None:
            value = np.ravel(model(*args[2:-1]) - meas)
        else:
            value = np.ravel(weights * (model(*args[2:-1]) - meas))

        if not np.all(np.isfinite(value)):
            raise NonFiniteValueError(
                "Objective function has encountered a non-finite value, "
                "this will cause the fit to fail!\n"
                "Please remove non-finite values from your input data before "
                "fitting to avoid this error."
            )

        return value

    @staticmethod
    def _add_fitting_uncertainties(model, cov_matrix):
        """
        Set ``cov_matrix`` and ``stds`` attributes on model with parameter
        covariance matrix returned by ``optimize.leastsq``.
        """
        free_param_names = [
            x
            for x in model.fixed
            if (model.fixed[x] is False) and (model.tied[x] is False)
        ]

        model.cov_matrix = Covariance(cov_matrix, free_param_names)
        model.stds = StandardDeviations(cov_matrix, free_param_names)

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
            fitter_to_model_params(model, params)
            if z is None:
                full = np.array(model.fit_deriv(x, *model.parameters))
                if not model.col_fit_deriv:
                    full_deriv = np.ravel(weights) * full.T
                else:
                    full_deriv = np.ravel(weights) * full
            else:
                full = np.array(
                    [np.ravel(_) for _ in model.fit_deriv(x, y, *model.parameters)]
                )
                if not model.col_fit_deriv:
                    full_deriv = np.ravel(weights) * full.T
                else:
                    full_deriv = np.ravel(weights) * full

            pars = [getattr(model, name) for name in model.param_names]
            fixed = [par.fixed for par in pars]
            tied = [par.tied for par in pars]
            tied = list(np.where([par.tied is not False for par in pars], True, tied))
            fix_and_tie = np.logical_or(fixed, tied)
            ind = np.logical_not(fix_and_tie)

            if not model.col_fit_deriv:
                residues = np.asarray(full_deriv[np.nonzero(ind)]).T
            else:
                residues = full_deriv[np.nonzero(ind)]

            return [np.ravel(_) for _ in residues]
        else:
            if z is None:
                fit_deriv = np.array(model.fit_deriv(x, *params))
                try:
                    output = np.array(
                        [np.ravel(_) for _ in np.array(weights) * fit_deriv]
                    )
                    if output.shape != fit_deriv.shape:
                        output = np.array(
                            [np.ravel(_) for _ in np.atleast_2d(weights).T * fit_deriv]
                        )
                    return output
                except ValueError:
                    return np.array(
                        [
                            np.ravel(_)
                            for _ in np.array(weights) * np.moveaxis(fit_deriv, -1, 0)
                        ]
                    ).transpose()
            else:
                if not model.col_fit_deriv:
                    return [
                        np.ravel(_)
                        for _ in (
                            np.ravel(weights)
                            * np.array(model.fit_deriv(x, y, *params)).T
                        ).T
                    ]
                return [
                    np.ravel(_)
                    for _ in weights * np.array(model.fit_deriv(x, y, *params))
                ]

    def _compute_param_cov(
        self, model, y, init_values, cov_x, fitparams, farg, weights=None
    ):
        # now try to compute the true covariance matrix
        if (len(y) > len(init_values)) and cov_x is not None:
            self.fit_info["param_cov"] = cov_x
            if weights is None:
                # if there are no measurement uncertainties given in `weights`,
                # fall back on the default behavior in scipy.optimize.curve_fit
                # when `absolute_sigma == False`. If there are uncertainties,
                # assume they are "absolute" and not "relative".
                # For details, see curve_fit:
                #   https://github.com/scipy/scipy/blob/
                #   c1ed5ece8ffbf05356a22a8106affcd11bd3aee0/scipy/
                #   optimize/_minpack_py.py#L591-L602
                sum_sqrs = np.sum(self.objective_function(fitparams, *farg) ** 2)
                dof = len(y) - len(init_values)
                self.fit_info["param_cov"] *= sum_sqrs / dof
        else:
            self.fit_info["param_cov"] = None
        if self._calc_uncertainties is True:
            if self.fit_info["param_cov"] is not None:
                self._add_fitting_uncertainties(model, self.fit_info["param_cov"])

    def _run_fitter(self, model, farg, maxiter, acc, epsilon, estimate_jacobian):
        return None, None, None

    def _filter_non_finite(self, x, y, z=None, weights=None):
        """
        Filter out non-finite values in x, y, z.

        Returns
        -------
        x, y, z : ndarrays
            x, y, and z with non-finite values filtered out.
        """
        MESSAGE = "Non-Finite input data has been removed by the fitter."

        mask = np.ones_like(x, dtype=bool) if weights is None else np.isfinite(weights)
        mask &= np.isfinite(y) if z is None else np.isfinite(z)

        if not np.all(mask):
            warnings.warn(MESSAGE, AstropyUserWarning)

        return (
            x[mask],
            y[mask],
            None if z is None else z[mask],
            None if weights is None else weights[mask],
        )

    @fitter_unit_support
    def __call__(
        self,
        model,
        x,
        y,
        z=None,
        weights=None,
        maxiter=DEFAULT_MAXITER,
        acc=DEFAULT_ACC,
        epsilon=DEFAULT_EPS,
        estimate_jacobian=False,
        filter_non_finite=False,
    ):
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
        z : array, optional
           input coordinates
        weights : array, optional
            Weights for fitting. For data with Gaussian uncertainties, the weights
            should be 1/sigma.

            .. versionchanged:: 5.3
                Calculate parameter covariances while accounting for ``weights``
                as "absolute" inverse uncertainties. To recover the old behavior,
                choose ``weights=None``.

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
        equivalencies : list or None, optional, keyword-only
            List of *additional* equivalencies that are should be applied in
            case x, y and/or z have units. Default is None.
        filter_non_finite : bool, optional
            Whether or not to filter data with non-finite values. Default is False

        Returns
        -------
        model_copy : `~astropy.modeling.FittableModel`
            a copy of the input model with parameters set by the fitter

        """
        model_copy = _validate_model(model, self.supported_constraints)
        model_copy.sync_constraints = False

        if filter_non_finite:
            x, y, z, weights = self._filter_non_finite(x, y, z, weights)
        farg = (
            model_copy,
            weights,
        ) + _convert_input(x, y, z)

        init_values, fitparams, cov_x = self._run_fitter(
            model_copy, farg, maxiter, acc, epsilon, estimate_jacobian
        )

        self._compute_param_cov(
            model_copy, y, init_values, cov_x, fitparams, farg, weights
        )

        model.sync_constraints = True
        return model_copy


class LevMarLSQFitter(_NonLinearLSQFitter):
    """
    Levenberg-Marquardt algorithm and least squares statistic.

    Parameters
    ----------
    calc_uncertainties : bool
        If the covarience matrix should be computed and set in the fit_info.
        Default: False

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

    def __init__(self, calc_uncertainties=False):
        super().__init__(calc_uncertainties)
        self.fit_info = {
            "nfev": None,
            "fvec": None,
            "fjac": None,
            "ipvt": None,
            "qtf": None,
            "message": None,
            "ierr": None,
            "param_jac": None,
            "param_cov": None,
        }

    def _run_fitter(self, model, farg, maxiter, acc, epsilon, estimate_jacobian):
        from scipy import optimize

        if model.fit_deriv is None or estimate_jacobian:
            dfunc = None
        else:
            dfunc = self._wrap_deriv
        init_values, _, _ = model_to_fit_params(model)
        fitparams, cov_x, dinfo, mess, ierr = optimize.leastsq(
            self.objective_function,
            init_values,
            args=farg,
            Dfun=dfunc,
            col_deriv=model.col_fit_deriv,
            maxfev=maxiter,
            epsfcn=epsilon,
            xtol=acc,
            full_output=True,
        )
        fitter_to_model_params(model, fitparams)
        self.fit_info.update(dinfo)
        self.fit_info["cov_x"] = cov_x
        self.fit_info["message"] = mess
        self.fit_info["ierr"] = ierr
        if ierr not in [1, 2, 3, 4]:
            warnings.warn(
                "The fit may be unsuccessful; check "
                "fit_info['message'] for more information.",
                AstropyUserWarning,
            )

        return init_values, fitparams, cov_x


class _NLLSQFitter(_NonLinearLSQFitter):
    """
    Wrapper class for `scipy.optimize.least_squares` method, which provides:
        - Trust Region Reflective
        - dogbox
        - Levenberg-Marqueardt
    algorithms using the least squares statistic.

    Parameters
    ----------
    method : str
        ‘trf’ :  Trust Region Reflective algorithm, particularly suitable
            for large sparse problems with bounds. Generally robust method.
        ‘dogbox’ : dogleg algorithm with rectangular trust regions, typical
            use case is small problems with bounds. Not recommended for
            problems with rank-deficient Jacobian.
        ‘lm’ : Levenberg-Marquardt algorithm as implemented in MINPACK.
            Doesn’t handle bounds and sparse Jacobians. Usually the most
            efficient method for small unconstrained problems.
    calc_uncertainties : bool
        If the covarience matrix should be computed and set in the fit_info.
        Default: False
    use_min_max_bounds: bool
        If the set parameter bounds for a model will be enforced each given
        parameter while fitting via a simple min/max condition. A True setting
        will replicate how LevMarLSQFitter enforces bounds.
        Default: False

    Attributes
    ----------
    fit_info :
        A `scipy.optimize.OptimizeResult` class which contains all of
        the most recent fit information
    """

    def __init__(self, method, calc_uncertainties=False, use_min_max_bounds=False):
        super().__init__(calc_uncertainties, use_min_max_bounds)
        self._method = method

    def _run_fitter(self, model, farg, maxiter, acc, epsilon, estimate_jacobian):
        from scipy import optimize
        from scipy.linalg import svd

        if model.fit_deriv is None or estimate_jacobian:
            dfunc = "2-point"
        else:

            def _dfunc(params, model, weights, x, y, z=None):
                if model.col_fit_deriv:
                    return np.transpose(
                        self._wrap_deriv(params, model, weights, x, y, z)
                    )
                else:
                    return self._wrap_deriv(params, model, weights, x, y, z)

            dfunc = _dfunc

        init_values, _, bounds = model_to_fit_params(model)

        # Note, if use_min_max_bounds is True we are defaulting to enforcing bounds
        # using the old method employed by LevMarLSQFitter, this is different
        # from the method that optimize.least_squares employs to enforce bounds
        # thus we override the bounds being passed to optimize.least_squares so
        # that it will not enforce any bounding.
        if self._use_min_max_bounds:
            bounds = (-np.inf, np.inf)

        self.fit_info = optimize.least_squares(
            self.objective_function,
            init_values,
            args=farg,
            jac=dfunc,
            max_nfev=maxiter,
            diff_step=np.sqrt(epsilon),
            xtol=acc,
            method=self._method,
            bounds=bounds,
        )

        # Adapted from ~scipy.optimize.minpack, see:
        # https://github.com/scipy/scipy/blob/47bb6febaa10658c72962b9615d5d5aa2513fa3a/scipy/optimize/minpack.py#L795-L816
        # Do Moore-Penrose inverse discarding zero singular values.
        _, s, VT = svd(self.fit_info.jac, full_matrices=False)
        threshold = np.finfo(float).eps * max(self.fit_info.jac.shape) * s[0]
        s = s[s > threshold]
        VT = VT[: s.size]
        cov_x = np.dot(VT.T / s**2, VT)

        fitter_to_model_params(model, self.fit_info.x, False)
        if not self.fit_info.success:
            warnings.warn(
                f"The fit may be unsuccessful; check: \n    {self.fit_info.message}",
                AstropyUserWarning,
            )

        return init_values, self.fit_info.x, cov_x


class TRFLSQFitter(_NLLSQFitter):
    """
    Trust Region Reflective algorithm and least squares statistic.

    Parameters
    ----------
    calc_uncertainties : bool
        If the covarience matrix should be computed and set in the fit_info.
        Default: False
    use_min_max_bounds: bool
        If the set parameter bounds for a model will be enforced each given
        parameter while fitting via a simple min/max condition. A True setting
        will replicate how LevMarLSQFitter enforces bounds.
        Default: False

    Attributes
    ----------
    fit_info :
        A `scipy.optimize.OptimizeResult` class which contains all of
        the most recent fit information
    """

    def __init__(self, calc_uncertainties=False, use_min_max_bounds=False):
        super().__init__("trf", calc_uncertainties, use_min_max_bounds)


class DogBoxLSQFitter(_NLLSQFitter):
    """
    DogBox algorithm and least squares statistic.

    Parameters
    ----------
    calc_uncertainties : bool
        If the covarience matrix should be computed and set in the fit_info.
        Default: False
    use_min_max_bounds: bool
        If the set parameter bounds for a model will be enforced each given
        parameter while fitting via a simple min/max condition. A True setting
        will replicate how LevMarLSQFitter enforces bounds.
        Default: False

    Attributes
    ----------
    fit_info :
        A `scipy.optimize.OptimizeResult` class which contains all of
        the most recent fit information
    """

    def __init__(self, calc_uncertainties=False, use_min_max_bounds=False):
        super().__init__("dogbox", calc_uncertainties, use_min_max_bounds)


class LMLSQFitter(_NLLSQFitter):
    """
    `scipy.optimize.least_squares` Levenberg-Marquardt algorithm and least squares statistic.

    Parameters
    ----------
    calc_uncertainties : bool
        If the covarience matrix should be computed and set in the fit_info.
        Default: False

    Attributes
    ----------
    fit_info :
        A `scipy.optimize.OptimizeResult` class which contains all of
        the most recent fit information
    """

    def __init__(self, calc_uncertainties=False):
        super().__init__("lm", calc_uncertainties, True)


class SLSQPLSQFitter(Fitter):
    """
    Sequential Least Squares Programming (SLSQP) optimization algorithm and
    least squares statistic.

    Raises
    ------
    ModelLinearityError
        A linear model is passed to a nonlinear fitter

    Notes
    -----
    See also the `~astropy.modeling.optimizers.SLSQP` optimizer.

    """

    supported_constraints = SLSQP.supported_constraints

    def __init__(self):
        super().__init__(optimizer=SLSQP, statistic=leastsquare)
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
        z : array, optional
            input coordinates
        weights : array, optional
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
        equivalencies : list or None, optional, keyword-only
            List of *additional* equivalencies that are should be applied in
            case x, y and/or z have units. Default is None.

        Returns
        -------
        model_copy : `~astropy.modeling.FittableModel`
            a copy of the input model with parameters set by the fitter

        """
        model_copy = _validate_model(model, self._opt_method.supported_constraints)
        model_copy.sync_constraints = False
        farg = _convert_input(x, y, z)
        farg = (
            model_copy,
            weights,
        ) + farg
        init_values, _, _ = model_to_fit_params(model_copy)
        fitparams, self.fit_info = self._opt_method(
            self.objective_function, init_values, farg, **kwargs
        )
        fitter_to_model_params(model_copy, fitparams)

        model_copy.sync_constraints = True
        return model_copy


class SimplexLSQFitter(Fitter):
    """
    Simplex algorithm and least squares statistic.

    Raises
    ------
    `ModelLinearityError`
        A linear model is passed to a nonlinear fitter

    """

    supported_constraints = Simplex.supported_constraints

    def __init__(self):
        super().__init__(optimizer=Simplex, statistic=leastsquare)
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
        z : array, optional
            input coordinates
        weights : array, optional
            Weights for fitting.
            For data with Gaussian uncertainties, the weights should be
            1/sigma.
        kwargs : dict
            optional keyword arguments to be passed to the optimizer or the statistic
        maxiter : int
            maximum number of iterations
        acc : float
            Relative error in approximate solution
        equivalencies : list or None, optional, keyword-only
            List of *additional* equivalencies that are should be applied in
            case x, y and/or z have units. Default is None.

        Returns
        -------
        model_copy : `~astropy.modeling.FittableModel`
            a copy of the input model with parameters set by the fitter

        """
        model_copy = _validate_model(model, self._opt_method.supported_constraints)
        model_copy.sync_constraints = False
        farg = _convert_input(x, y, z)
        farg = (
            model_copy,
            weights,
        ) + farg

        init_values, _, _ = model_to_fit_params(model_copy)

        fitparams, self.fit_info = self._opt_method(
            self.objective_function, init_values, farg, **kwargs
        )
        fitter_to_model_params(model_copy, fitparams)
        model_copy.sync_constraints = True
        return model_copy


class JointFitter(metaclass=_FitterMeta):
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
        self.fitparams = self.model_to_fit_params()

        # a list of model.n_inputs
        self.modeldims = [m.n_inputs for m in self.models]
        # sum all model dimensions
        self.ndim = np.sum(self.modeldims)

    def model_to_fit_params(self):
        fparams = []
        fparams.extend(self.initvals)
        for model in self.models:
            params = model.parameters.tolist()
            joint_params = self.jointparams[model]
            param_metrics = model._param_metrics
            for param_name in joint_params:
                slice_ = param_metrics[param_name]["slice"]
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
            margs = lstsqargs[: model.n_inputs + 1]
            del lstsqargs[: model.n_inputs + 1]
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
                    slice_ = param_metrics[param_name]["slice"]
                    plen = slice_.stop - slice_.start
                    mparams.extend(mfparams[:plen])
                    del mfparams[:plen]
            modelfit = model.evaluate(margs[:-1], *mparams)
            fitted.extend(modelfit - margs[-1])
        return np.ravel(fitted)

    def _verify_input(self):
        if len(self.models) <= 1:
            raise TypeError(f"Expected >1 models, {len(self.models)} is given")
        if len(self.jointparams.keys()) < 2:
            raise TypeError(
                "At least two parameters are expected, "
                f"{len(self.jointparams.keys())} is given"
            )
        for j in self.jointparams.keys():
            if len(self.jointparams[j]) != len(self.initvals):
                raise TypeError(
                    f"{len(self.jointparams[j])} parameter(s) "
                    f"provided but {len(self.initvals)} expected"
                )

    def __call__(self, *args):
        """
        Fit data to these models keeping some of the parameters common to the
        two models.
        """
        from scipy import optimize

        if len(args) != reduce(lambda x, y: x + 1 + y + 1, self.modeldims):
            raise ValueError(
                f"Expected {reduce(lambda x, y: x + 1 + y + 1, self.modeldims)} "
                f"coordinates in args but {len(args)} provided"
            )

        self.fitparams[:], _ = optimize.leastsq(
            self.objective_function, self.fitparams, args=args
        )

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
                    slice_ = param_metrics[param_name]["slice"]
                    plen = slice_.stop - slice_.start
                    mparams.extend(mfparams[:plen])
                    del mfparams[:plen]
            model.parameters = np.array(mparams)


def _convert_input(x, y, z=None, n_models=1, model_set_axis=0):
    """Convert inputs to float arrays."""
    x = np.asanyarray(x, dtype=float)
    y = np.asanyarray(y, dtype=float)

    if z is not None:
        z = np.asanyarray(z, dtype=float)
        data_ndim, data_shape = z.ndim, z.shape
    else:
        data_ndim, data_shape = y.ndim, y.shape

    # For compatibility with how the linear fitter code currently expects to
    # work, shift the dependent variable's axes to the expected locations
    if n_models > 1 or data_ndim > x.ndim:
        if (model_set_axis or 0) >= data_ndim:
            raise ValueError("model_set_axis out of range")
        if data_shape[model_set_axis] != n_models:
            raise ValueError(
                "Number of data sets (y or z array) is expected to equal "
                "the number of parameter sets"
            )
        if z is None:
            # For a 1-D model the y coordinate's model-set-axis is expected to
            # be last, so that its first dimension is the same length as the x
            # coordinates.  This is in line with the expectations of
            # numpy.linalg.lstsq:
            # https://numpy.org/doc/stable/reference/generated/numpy.linalg.lstsq.html
            # That is, each model should be represented by a column.  TODO:
            # Obviously this is a detail of np.linalg.lstsq and should be
            # handled specifically by any fitters that use it...
            y = np.rollaxis(y, model_set_axis, y.ndim)
            data_shape = y.shape[:-1]
        else:
            # Shape of z excluding model_set_axis
            data_shape = z.shape[:model_set_axis] + z.shape[model_set_axis + 1 :]

    if z is None:
        if data_shape != x.shape:
            raise ValueError("x and y should have the same shape")
        farg = (x, y)
    else:
        if not (x.shape == y.shape == data_shape):
            raise ValueError("x, y and z should have the same shape")
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
def fitter_to_model_params(model, fps, use_min_max_bounds=True):
    """
    Constructs the full list of model parameters from the fitted and
    constrained parameters.

    Parameters
    ----------
    model :
        The model being fit
    fps :
        The fit parameter values to be assigned
    use_min_max_bounds: bool
        If the set parameter bounds for model will be enforced on each
        parameter with bounds.
        Default: True
    """
    _, fit_param_indices, _ = model_to_fit_params(model)

    has_tied = any(model.tied.values())
    has_fixed = any(model.fixed.values())
    has_bound = any(b != (None, None) for b in model.bounds.values())
    parameters = model.parameters

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

        slice_ = param_metrics[name]["slice"]
        shape = param_metrics[name]["shape"]
        # This is determining which range of fps (the fitted parameters) maps
        # to parameters of the model
        size = reduce(operator.mul, shape, 1)

        values = fps[offset : offset + size]

        # Check bounds constraints
        if model.bounds[name] != (None, None) and use_min_max_bounds:
            _min, _max = model.bounds[name]
            if _min is not None:
                values = np.fmax(values, _min)
            if _max is not None:
                values = np.fmin(values, _max)

        parameters[slice_] = values
        offset += size

    # Update model parameters before calling ``tied`` constraints.
    model._array_to_parameters()

    # This has to be done in a separate loop due to how tied parameters are
    # currently evaluated (the fitted parameters need to actually be *set* on
    # the model first, for use in evaluating the "tied" expression--it might be
    # better to change this at some point
    if has_tied:
        for idx, name in enumerate(model.param_names):
            if model.tied[name]:
                value = model.tied[name](model)
                slice_ = param_metrics[name]["slice"]

                # To handle multiple tied constraints, model parameters
                # need to be updated after each iteration.
                parameters[slice_] = value
                model._array_to_parameters()


@deprecated(
    since="5.1",
    message="private method: _fitter_to_model_params has been made public now",
)
def _fitter_to_model_params(model, fps):
    return fitter_to_model_params(model, fps)


def model_to_fit_params(model):
    """
    Convert a model instance's parameter array to an array that can be used
    with a fitter that doesn't natively support fixed or tied parameters.
    In particular, it removes fixed/tied parameters from the parameter
    array.
    These may be a subset of the model parameters, if some of them are held
    constant or tied.
    """
    fitparam_indices = list(range(len(model.param_names)))
    model_params = model.parameters
    model_bounds = list(model.bounds.values())
    if any(model.fixed.values()) or any(model.tied.values()):
        params = list(model_params)
        param_metrics = model._param_metrics
        for idx, name in list(enumerate(model.param_names))[::-1]:
            if model.fixed[name] or model.tied[name]:
                slice_ = param_metrics[name]["slice"]
                del params[slice_]
                del model_bounds[slice_]
                del fitparam_indices[idx]
        model_params = np.array(params)

    for idx, bound in enumerate(model_bounds):
        if bound[0] is None:
            lower = -np.inf
        else:
            lower = bound[0]

        if bound[1] is None:
            upper = np.inf
        else:
            upper = bound[1]

        model_bounds[idx] = (lower, upper)
    model_bounds = tuple(zip(*model_bounds))
    return model_params, fitparam_indices, model_bounds


@deprecated(
    since="5.1",
    message="private method: _model_to_fit_params has been made public now",
)
def _model_to_fit_params(model):
    return model_to_fit_params(model)


def _validate_constraints(supported_constraints, model):
    """Make sure model constraints are supported by the current fitter."""
    message = "Optimizer cannot handle {0} constraints."

    if any(model.fixed.values()) and "fixed" not in supported_constraints:
        raise UnsupportedConstraintError(message.format("fixed parameter"))

    if any(model.tied.values()) and "tied" not in supported_constraints:
        raise UnsupportedConstraintError(message.format("tied parameter"))

    if (
        any(tuple(b) != (None, None) for b in model.bounds.values())
        and "bounds" not in supported_constraints
    ):
        raise UnsupportedConstraintError(message.format("bound parameter"))

    if model.eqcons and "eqcons" not in supported_constraints:
        raise UnsupportedConstraintError(message.format("equality"))

    if model.ineqcons and "ineqcons" not in supported_constraints:
        raise UnsupportedConstraintError(message.format("inequality"))


def _validate_model(model, supported_constraints):
    """
    Check that model and fitter are compatible and return a copy of the model.
    """
    if not model.fittable:
        raise ValueError("Model does not appear to be fittable.")
    if model.linear:
        warnings.warn(
            "Model is linear in parameters; consider using linear fitting methods.",
            AstropyUserWarning,
        )
    elif len(model) != 1:
        # for now only single data sets ca be fitted
        raise ValueError("Non-linear fitters can only fit one data set at a time.")
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
    entry_points : list of `~importlib.metadata.EntryPoint`
        entry_points are objects which encapsulate importable objects and
        are defined on the installation of a package.

    Notes
    -----
    An explanation of entry points can be found `here
    <http://setuptools.readthedocs.io/en/latest/setuptools.html#dynamic-discovery-of-services-and-plugins>`_
    """
    for entry_point in entry_points:
        name = entry_point.name
        try:
            entry_point = entry_point.load()
        except Exception as e:
            # This stops the fitting from choking if an entry_point produces an error.
            warnings.warn(
                AstropyUserWarning(
                    f"{type(e).__name__} error occurred in entry point {name}."
                )
            )
        else:
            if not inspect.isclass(entry_point):
                warnings.warn(
                    AstropyUserWarning(
                        f"Modeling entry point {name} expected to be a Class."
                    )
                )
            else:
                if issubclass(entry_point, Fitter):
                    name = entry_point.__name__
                    globals()[name] = entry_point
                    __all__.append(name)
                else:
                    warnings.warn(
                        AstropyUserWarning(
                            f"Modeling entry point {name} expected to extend "
                            "astropy.modeling.Fitter"
                        )
                    )


def _populate_ep():
    # TODO: Exclusively use select when Python minversion is 3.10
    ep = entry_points()
    if hasattr(ep, "select"):
        populate_entry_points(ep.select(group="astropy.modeling"))
    else:
        populate_entry_points(ep.get("astropy.modeling", []))


_populate_ep()
