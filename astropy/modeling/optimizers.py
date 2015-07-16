# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Optimization algorithms used in `~astropy.modeling.fitting`.
"""

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)
import warnings
import abc
import numpy as np
from ..extern import six
from ..utils.exceptions import AstropyUserWarning

__all__ = ["Optimization", "SLSQP", "Simplex"]

# Maximum number of iterations
DEFAULT_MAXITER = 100

# Step for the forward difference approximation of the Jacobian
DEFAULT_EPS = np.sqrt(np.finfo(float).eps)

#Default requested accuracy
DEFAULT_ACC = 1e-07

DEFAULT_BOUNDS = (-10 ** 12, 10 ** 12)


@six.add_metaclass(abc.ABCMeta)
class Optimization(object):
    """
    Base class for optimizers.

    Parameters
    ----------
    opt_method : callable
        Implements optimization method

    Notes
    -----
    The base Optimizer does not support any constraints by default; individual
    optimizers should explicitly set this list to the specific constraints
    it supports.

    """

    supported_constraints = []

    def __init__(self, opt_method):
        self._opt_method = opt_method
        self._maxiter = DEFAULT_MAXITER
        self._eps = DEFAULT_EPS
        self._acc = DEFAULT_ACC

    @property
    def maxiter(self):
        """Maximum number of iterations"""
        return self._maxiter

    @maxiter.setter
    def maxiter(self, val):
        """Set maxiter"""
        self._maxiter = val

    @property
    def eps(self):
        """Step for the forward difference approximation of the Jacobian"""
        return self._eps

    @eps.setter
    def eps(self, val):
        """Set eps value"""
        self._eps = val

    @property
    def acc(self):
        """Requested accuracy"""
        return self._acc

    @acc.setter
    def acc(self, val):
        """Set accuracy"""
        self._acc = val

    def __repr__(self):
        fmt = "{0}()".format(self.__class__.__name__)
        return fmt

    @property
    def opt_method(self):
        return self._opt_method

    @abc.abstractmethod
    def __call__(self):
        raise NotImplementedError("Subclasses should implement this method")


class SLSQP(Optimization):
    """
    Sequential Least Squares Programming optimization algorithm.

    The algorithm is described in [1]_. It supports tied and fixed
    parameters, as well as bounded constraints. Uses
    `scipy.optimize.fmin_slsqp`.

    References
    ----------
    .. [1] http://www.netlib.org/toms/733
    """
    supported_constraints = ['bounds', 'eqcons', 'ineqcons', 'fixed', 'tied']

    def __init__(self):
        from scipy.optimize import fmin_slsqp
        super(SLSQP, self).__init__(fmin_slsqp)
        self.fit_info = {
            'final_func_val': None,
            'numiter': None,
            'exit_mode': None,
            'message': None
        }

    def __call__(self, objfunc, initval, fargs, **kwargs):
        """
        Run the solver.

        Parameters
        ----------
        objfunc : callable
            objection function
        initval : iterable
            initial guess for the parameter values
        fargs : tuple
            other arguments to be passed to the statistic function
        kwargs : dict
            other keyword arguments to be passed to the solver

        """
        kwargs['iter'] = kwargs.pop('maxiter', self._maxiter)

        if 'epsilon' not in kwargs:
            kwargs['epsilon'] = self._eps
        if 'acc' not in kwargs:
            kwargs['acc'] = self._acc

        # set the values of constraints to match the requirements of fmin_slsqp
        model = fargs[0]
        pars = [getattr(model, name) for name in model.param_names]
        bounds = [par.bounds for par in pars if par.fixed != True and
                  par.tied == False]
        bounds = np.asarray(bounds)
        for i in bounds:
            if i[0] is None:
                i[0] =  DEFAULT_BOUNDS[0]
            if i[1] is None:
                i[1] = DEFAULT_BOUNDS[1]
        # older versions of scipy require this array to be float
        bounds = np.asarray(bounds, dtype=np.float)
        eqcons = np.array(model.eqcons)
        ineqcons = np.array(model.ineqcons)
        fitparams, final_func_val, numiter, exit_mode, mess = self.opt_method(
            objfunc, initval, args=fargs, full_output=True,
            bounds=bounds, eqcons=eqcons, ieqcons=ineqcons,
            **kwargs)

        self.fit_info['final_func_val'] = final_func_val
        self.fit_info['numiter'] = numiter
        self.fit_info['exit_mode'] = exit_mode
        self.fit_info['message'] = mess

        if exit_mode != 0:
            warnings.warn("The fit may be unsuccessful; check "
                          "fit_info['message'] for more information.",
                          AstropyUserWarning)

        return fitparams, self.fit_info


class Simplex(Optimization):
    """
    Neald-Mead (downhill simplex) algorithm [1].

    This algorithm only uses function values, not derivatives.
    Uses `scipy.optimize.fmin`.

    .. [1] Nelder, J.A. and Mead, R. (1965), "A simplex method for function
           minimization", The Computer Journal, 7, pp. 308-313
    """

    supported_constraints = ['bounds', 'fixed', 'tied']

    def __init__(self):
        from scipy.optimize import fmin as simplex
        super(Simplex, self).__init__(simplex)
        self.fit_info = {
            'final_func_val': None,
            'numiter': None,
            'exit_mode': None,
            'num_function_calls': None
        }

    def __call__(self, objfunc, initval, fargs, **kwargs):
        """
        Run the solver.

        Parameters
        ----------
        objfunc : callable
            objection function
        initval : iterable
            initial guess for the parameter values
        fargs : tuple
            other arguments to be passed to the statistic function
        kwargs : dict
            other keyword arguments to be passed to the solver

        """
        if 'maxiter' not in kwargs:
            kwargs['maxiter'] = self._maxiter
        if 'acc' in kwargs:
            self._acc = kwargs['acc']
            kwargs.pop['acc']
        if 'xtol' in kwargs:
            self._acc = kwargs['xtol']
            kwargs.pop['xtol']
  
        fitparams, final_func_val, numiter, funcalls, exit_mode = self.opt_method(
            objfunc, initval, args=fargs, xtol=self._acc,
            full_output=True, **kwargs)
        self.fit_info['final_func_val'] = final_func_val
        self.fit_info['numiter'] = numiter
        self.fit_info['exit_mode'] = exit_mode
        self.fit_info['num_function_calls'] = funcalls
        if self.fit_info['exit_mode'] == 1:
            warnings.warn("The fit may be unsuccessful; "
                          "Maximum number of function evaluations reached.",
                          AstropyUserWarning)
        if self.fit_info['exit_mode'] == 2:
            warnings.warn("The fit may be unsuccessful; "
                          "Maximum number of iterations reached.",
                          AstropyUserWarning)
        return fitparams, self.fit_info
