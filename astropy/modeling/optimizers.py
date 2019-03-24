# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Optimization algorithms used in `~astropy.modeling.fitting`.
"""

import warnings
import abc
import numpy as np
from astropy.utils.exceptions import AstropyUserWarning

__all__ = ["Optimization", "SLSQP", "Simplex", "Minimize"]

# Maximum number of iterations
DEFAULT_MAXITER = 100

# Step for the forward difference approximation of the Jacobian
DEFAULT_EPS = np.sqrt(np.finfo(float).eps)

# Default requested accuracy
DEFAULT_ACC = 1e-07

DEFAULT_BOUNDS = (-10 ** 12, 10 ** 12)


class Optimization(metaclass=abc.ABCMeta):
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
        super().__init__(fmin_slsqp)
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
        # Get the verbosity level
        disp = kwargs.pop('verblevel', None)

        # set the values of constraints to match the requirements of fmin_slsqp
        model = fargs[0]
        pars = [getattr(model, name) for name in model.param_names]
        bounds = [par.bounds for par in pars if not (par.fixed or par.tied)]
        bounds = np.asarray(bounds)
        for i in bounds:
            if i[0] is None:
                i[0] = DEFAULT_BOUNDS[0]
            if i[1] is None:
                i[1] = DEFAULT_BOUNDS[1]
        # older versions of scipy require this array to be float
        bounds = np.asarray(bounds, dtype=float)
        eqcons = np.array(model.eqcons)
        ineqcons = np.array(model.ineqcons)
        fitparams, final_func_val, numiter, exit_mode, mess = self.opt_method(
            objfunc, initval, args=fargs, full_output=True, disp=disp,
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
    Neald-Mead (downhill simplex) algorithm.

    This algorithm [1]_ only uses function values, not derivatives.
    Uses `scipy.optimize.fmin`.

    References
    ----------
    .. [1] Nelder, J.A. and Mead, R. (1965), "A simplex method for function
       minimization", The Computer Journal, 7, pp. 308-313
    """

    supported_constraints = ['bounds', 'fixed', 'tied']

    def __init__(self):
        from scipy.optimize import fmin as simplex
        super().__init__(simplex)
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
            kwargs.pop('acc')
        if 'xtol' in kwargs:
            self._acc = kwargs['xtol']
            kwargs.pop('xtol')
        # Get the verbosity level
        disp = kwargs.pop('verblevel', None)

        fitparams, final_func_val, numiter, funcalls, exit_mode = self.opt_method(
            objfunc, initval, args=fargs, xtol=self._acc, disp=disp,
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


class Minimize(Optimization):
    """
    Interface for `scipy.optimize.minimize`.


    Parameters
    ----------
    method : str or callable(custom)
        Selects minimization method.

    supported_constraints : list
        The constraint types supported by the Optimization method. By default fixed and tied parameters are supported but
        these can be set for custom methods. For defined methods this will be changed automatically.


    """

    jac_required = False

    def __init__(self, method='slsqp', supported_constraints=None):
        from scipy.optimize import minimize
        super(Minimize, self).__init__(minimize)
        method = method.lower()
        if callable(method):
            if supported_constraints is None:
                warnings.warn("No supported_constraints set, by default the optimizer will only support fixed and tied parameters."
                              "If you wish this optimizer to support more then you change this", AstropyUserWarning)

                self.supported_constraints = ['fixed', 'tied']
            else:
                self.supported_constraints = supported_constraints
        else:
            if method in ['nelder-mead', 'powell', 'cg', 'bfgs', 'newton-cg', 'dogleg', 'trust-ncg']:
                self.supported_constraints = ['fixed', 'tied']
            elif method in ['l-bfgs-b', 'tnc']:
                self.supported_constraints = ['fixed', 'tied', 'bounds']
            elif method == 'cobyla':
                self.supported_constraints = ['fixed', 'tied', 'ineqcons']
            elif method in ['slsqp']:
                self.supported_constraints = ['fixed', 'tied', 'eqcons', 'ineqcons', 'bounds']

            if ['newton-cg', 'trust-ncg']:
                self.jac_required = True

        self._method = method
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
            kwargs['options'] = {'maxiter': self._maxiter}
        else:
            kwargs['options'] = {'maxiter': kwargs['maxiter']}
            kwargs.pop('maxiter')

        if 'acc' in kwargs:
            self._acc = kwargs['acc']
            kwargs.pop('acc')
        if 'xtol' in kwargs:
            self._acc = kwargs['xtol']
            kwargs.pop('xtol')

        model = fargs[0]
        pars = [getattr(model, name) for name in model.param_names]

        if 'bounds' in self.supported_constraints:
            bounds = [par.bounds for par in pars if (par.fixed is not True and par.tied is False)]
            bounds = np.asarray(bounds)
            for i in bounds:
                if i[0] is None:
                    i[0] = DEFAULT_BOUNDS[0]
                if i[1] is None:
                    i[1] = DEFAULT_BOUNDS[1]
                # older versions of scipy require this array to be float
                kwargs['bounds'] = np.asarray(bounds, dtype=np.float)

        if 'eqcons' in self.supported_constraints and np.array(model.eqcons)>0:
            if not 'constraints' in kwargs:
                kwargs['constraints']=[]
            for eq in model.eqcons:
                kwargs['constraints'].append({"type":"eq","fun":eq })

        if 'ineqcons' in self.supported_constraints and np.array(model.ineqcons)>0:
            if not 'constraints' in kwargs:
                kwargs['constraints']=[]
            for ineq in model.ineqcons:
                kwargs['constraints'].append({"type":"ineq","fun":ineq })

        res = self.opt_method(objfunc, initval, method=self._method, args=fargs, tol=self._acc, **kwargs)
        self.scipy_opt_result = res

        if res['status'] == 1:
            warnings.warn("The fit may be unsuccessful; "
                          "Maximum number of function evaluations reached.",
                          AstropyUserWarning)
        if res['status'] == 2:
            warnings.warn("The fit may be unsuccessful; "
                          "Maximum number of iterations reached.",
                          AstropyUserWarning)
        return res['x'], {info_name: res[res_name] for info_name, res_name in zip(['final_func_val', 'numiter', 'num_function_calls', 'exit_mode'],
                                                                                  ['fun', 'nit', 'nfev', 'status']) if res_name in res}
