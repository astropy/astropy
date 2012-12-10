# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module provides wrappers, called Fitters, around some Numpy and Scipy
fitting functions. All Fitters take an instance of `ParametricModel` as input
and define a __call__ method which fits the model to the data and changes the
model's parameters attribute. The idea is to make this extensible and allow
users to easily add other fitters.

Linear fitting is done using Numpy's linalg.lstsq function.
There are currently two non-linear fitters which use leastsq and slsqp functions
in scipy.optimize.

"""
from __future__ import division, print_function
import abc
import numpy as np
from numpy import linalg

from scipy import optimize
import operator
from .util import pmapdomain

__all__ = ['LinearLSQFitter', 'NonLinearLSQFitter', 'SLSQPFitter', 
           'JointFitter']

MAXITER = 100
EPS = np.sqrt(np.finfo(float).eps)

constraintsdef = {'NonLinearLSQFitter': ['fixed', 'tied', 'bounds'],
                  'SLSQPFitter': ['bounds', 'eqcons', 'ineqcons', 'fixed', 'tied'],
                  'LinearLSQFitter': ['fixed'],
                  }

class ModelLinearityException(Exception):
    """
    Called when a linear model is passed to a non-linear fitter and vice versa.
    """
    pass
        
class Fitter(object):
    """
    Base class for all fitters.
    
    The purpose of this class is to deal with constraints
    """
    def __init__(self, model):
        __metaclass__ = abc.ABCMeta
        
        self.model = model
        self._validate_constraints()
        if any(self.model.constraints._fixed.values()) or \
           any(self.model.constraints._tied.values()):
            self._fitpars = self.model.constraints.fitpars[:]
        else:
            self._fitpars = self.model._parameters[:]
        if self.model.deriv is None:
            self.dfunc = None
        else:
            self.dfunc = self._wrap_deriv
    
    @property
    def fitpars(self):
        if any(self.model.constraints._fixed.values()) or \
           any(self.model.constraints._tied.values()):
            return self.model.constraints.fitpars
        else:
            return self.model._parameters
        
    @fitpars.setter
    def fitpars(self, fps, set_bounds=None):
        """
        Returns a list of parameters to be passed to the fitting
        algorithm. This is either the model.parameters (if there are
        no constrints or a modified version of model.parameters which takes
        into account constraints.

        Different fitters deal with bounds in a different way. Some set bounds 
        internally as part of the algorithm and some don't. 
        If a function set_bounds is provided, bounds will be dealt with here.
        This function should be set in the individual fitters.
        Parameters
        ----------
        set_bounds: callable
            
        """
        if any(self.model.constraints._fixed.values()) or \
           any(self.model.constraints._tied.values()):
            self.model.constraints.fitpars = fps
            self._fitpars[:] = self.model.constraints.fitpars
        else:
            try:
                self._set_bounds(fps)
            except NotImplementedError:
                pass
            self._fitpars[:] = fps
            self.model._parameters._update(fps)

    def _set_bounds(self, pars):
        raise NotImplementedError
    
    def _wrap_deriv(self, p, x, y, z=None):
        fixed_and_tied = [name for name in self.model.constraints.fixed if  self.model.constraints.fixed[name]]
        fixed_and_tied.extend([name for name in self.model.constraints.tied if  self.model.constraints.tied[name]])
        if fixed_and_tied:
            pars =  self.model.constraints.modelpars
            if z is None:
                fullderiv = self.model.deriv(pars, x, y)
            else:
                fullderiv = self.model.deriv(pars, x, y, z)
            ind = range(len(self.model.parnames))
            for name in fixed_and_tied:
                index = self.model.parnames.index(name)
                ind.remove(index)
            res = np.empty((fullderiv.shape[0], fullderiv.shape[1]-len(ind)))
            res = fullderiv[:, ind]
            return res
        else:
            pars = self.model._parameters
            if z is None:
                return self.model.deriv(pars, x, y)
            else:
                return self.model.deriv(pars, x, y, z)

        if z is None:
            res = self.model.deriv(pars, x, y)
        else:
            res = self.model.deriv(pars, x, y, z)
        return res
     
    def _validate_constraints(self):
        fname = self.__class__.__name__
        try:
            c = constraintsdef[fname]
        except KeyError:
            print("%s does not support fitting with constraints" % fname)
            raise
        if any(self.model.constraints._fixed.values()) and 'fixed' not in c:
            raise ValueError("%s cannot handle fixed parameter constraints "\
                                            % fname)
        if any(self.model.constraints._tied.values()) and 'tied' not in c:
            raise ValueError("%s cannot handle tied parameter constraints "\
                                            % fname)
        if any(c != (-1E12, 1E12) for c in self.model.constraints._bounds.values()) \
                and 'bounds' not in c:            
            raise ValueError("%s cannot handle bound parameter constraints"
                             % fname)
        if self.model.constraints._eqcons and 'eqcons' not in c:
            raise ValueError("%s cannot handle equality constraints but "
                             "eqcons given" % fname)
        if self.model.constraints._ineqcons and 'ineqcons' not in c:
            raise ValueError("%s cannot handle inequality constraints but "
                             "ineqcons given" % fname)
    
    @property
    def covar(self):
        return None
    
    @abc.abstractmethod
    def __call__(self):
        raise NotImplementedError

class LinearLSQFitter(Fitter):
    """
    A class representing a linear least square fitting
    
    Uses numpy.linalg.lstsq to do the fitting.
    Given a model and data, fits the model to the data and changes the
    model's parameters. Keeps a dictionary of auxiliary fitting information.
    
    """
    def __init__(self, model):
        """
        Parameters
        ----------
        model: an instance of `fitting.models.ParametricModel`
        
        """
        self.model = model
        if not self.model.linear:
            raise ModelLinearityException('Model is not linear in parameters, '
                                    'linear fit methods should not be used.')
        self.fit_info = {'residuals': None,
                         'rank': None,
                         'singular_values': None,
                         'pars': None
                        }
        super(LinearLSQFitter, self).__init__(model)
        
    def _deriv_with_constraints(self, x, y=None):
        if y is None:
            d = self.model.deriv(x)
        else:
            d = self.model.deriv(x, y)
        fixed = [name for name in self.model.constraints.fixed if  self.model.constraints.fixed[name]]
        ind = range(len(self.model.parnames))
        for name in fixed:
            index = self.model.parnames.index(name)
            ind.remove(index)
        res = d[:,ind]
        return res
        
    def __call__(self, x, y, z=None, w=None, rcond=None):
        """
        Parameters
        ----------
        x: array
           input coordinates
        y: array
           input coordinates
        z: array (optional)
           input coordinates
        w: array (optional)
           weights
        rcond:  float, optional
                Cut-off ratio for small singular values of `a`.
                Singular values are set to zero if they are smaller than `rcond`
                times the largest singular value of `a`.
       
        """
        multiple = False
        x = np.asarray(x) + 0.0
        y = np.asarray(y) + 0.0
        
        if self.model.ndim == 2 and z is None:
            raise ValueError("Expected x, y and z for a 2 dimensional model.")
        
        if z is None:
            if x.shape[0] != y.shape[0]:
                raise ValueError("x and y should have the same length")
            if y.ndim == 2:
                assert y.shape[1] == self.model._parameters.paramdim, (
                    "Number of data sets (Y array is expected to equal "
                    "the number of parameter sets")
            # map domain into window
            if not self.model.domain:
                self.model.domain = [x.min(), x.max()]
            if not self.model.window:
                self.model.window = [-1, 1]
            xnew = pmapdomain(x, self.model.domain, self.model.window)
            if any(self.model.constraints.fixed.values()):
                lhs = self._deriv_with_constraints(xnew)
            else:
                lhs = self.model.deriv(xnew)
            if len(y.shape) == 2:
                rhs = y
                multiple = y.shape[1]
            else:
                rhs = y
        else:
            if x.shape != y.shape:
                raise ValueError("Expected x and y to have the same shape")
            if x.shape[-1] != z.shape[-1]:
                raise ValueError("x and z should have equal last dimensions")
            
            # map domain into window
            if self.model.xdomain is None:
                self.model.xdomain = [x.min(), x.max()]
            if self.model.ydomain is None:
                self.model.ydomain = [y.min(), y.max()]
            if self.model.xwindow is None:
                self.model.xwindow = [-1., 1.]       
            if self.model.ywindow is None:
                self.model.ywindow = [-1., 1.]
            
            xnew = pmapdomain(x, self.model.xdomain, self.model.xwindow)
            ynew = pmapdomain(y, self.model.ydomain, self.model.ywindow)
            
            if any(self.model.constraints.fixed.values()):
                lhs = self._deriv_with_constraints(xnew, ynew)
            else:
                lhs = self.model.deriv(xnew, ynew)
            if len(z.shape) == 3:
                rhs = np.array([i.flatten() for i in z]).T
                multiple = z.shape[0]
            else:
                rhs = z.flatten()

        if w is not None:
            w = np.asarray(w) + 0.0
            if len(x) != len(y):
                raise ValueError("x and w should have the same length")
            if rhs.ndim == 2:
                lhs *= w[:, np.newaxis]
                rhs *= w[:, np.newaxis]
            else:
                lhs *= w[:, np.newaxis]
                rhs *= w

        if not multiple and self.model._parameters.paramdim > 1:
            raise ValueError("Attempting to fit a 1D data set to a model "
                             "with multiple parameter sets")
        if rcond is None:
            rcond = len(x)*np.finfo(x.dtype).eps
        
        scl = (lhs*lhs).sum(0)
        lacoef, resids, rank, sval = linalg.lstsq(lhs/scl, rhs, rcond)
        
        self.fit_info['residuals'] = resids
        self.fit_info['rank'] = rank
        self.fit_info['singular_values'] = sval
        
        self.model._parameters._changed = True
        # If y.ndim > model.ndim we are doing a simultanious 1D fitting 
        # of several 1D arrays. Otherwise the model is 2D.
        #if y.ndim > self.model.ndim:
        if multiple:
            self.model._parameters.paramdim = multiple
        lacoef = (lacoef.T/scl).T
        self.fit_info['pars'] = lacoef
        if rank != self.model._order:
            print("The fit may be poorly conditioned\n")
        self.fitpars = lacoef.flatten()[:]
        
class NonLinearLSQFitter(Fitter):
    """
    A wrapper around scipy.optimize.leastsq
    Supports tied and frozen parameters.
    
    """
    def __init__(self, model):
        """
        Parameters
        ----------
        model: a fittable :class: `models.ParametricModel`
               model to fit to data
        fixed: iterable
               a tuple of parameter names to be held fixed during fitting
        tied:  dict
               keys are parameter names
               values are callable/function providing a relationship
               between parameters. Currently the callable takes a model 
               instance as an argument.
               In the example below xcen is tied to the value of xsigma
               
               def tie_center(model):
                   xcen = 50*model.xsigma
                   return xcen
        
               tied ={'xcen':tie_center}
        
        Raises
        ------
        ModelLinearityException 
            A linear model is passed to a nonlinear fitter
            
        """
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
            raise ModelLinearityException('Model is linear in parameters, '
                            'non-linear fitting methods should not be used.')
            
    def errorfunc(self, fps, *args):
        self.fitpars = fps
        meas = args[0]
        return np.ravel(self.model(*args[1:]) - meas)
    
    def _set_bounds(self, fitpars):
        if any(c != (-1E12, 1E12) for c in self.model.constraints.bounds.values()):
            bounds=[self.model.constraints.bounds[par] for par in self.model.parnames]
            [setattr(self.model, name, par if par>b[0] else b[0]) for name, par, b in \
                zip(self.model.parnames, fitpars, bounds)]
            [setattr(self.model, name, par if par<b[1] else b[1]) for name, par, b in \
                zip(self.model.parnames, fitpars, bounds)]
    
    @property
    def covar(self):
        """
        Calculate the covariance matrix 
        (doesn't take into account constraints)
        """
        n = len(self.model.parameters)
        # construct the permutation matrix
        P = np.take(np.eye(n), self.fit_info['ipvt']-1, 0)
        # construct the R matrix as in JP = QR
        r = np.triu(self.fit_info['fjac'].T[:n, :])
        r_pt = np.dot(r, P.T)
        p_rt = np.dot(P, r.T)
        try:
            return np.dual.inv(np.dot(p_rt, r_pt))
        except:
            print("Could not construct a covariance matrix")
            return None

    def __call__(self, x, y, z=None, w=None, maxiter=MAXITER, epsilon=EPS):
        """
        Parameters
        ----------
        x: array
           input coordinates
        y: array
           input coordinates
        z: array (optional)
           input coordinates
        w: array (optional
           weights
        maxiter: int
            maximum number of iterations
        epsilon: float
            A suitable step length for the forward-difference 
            approximation of the Jacobian (if model.fjac=None). If 
            epsfcn is less than the machine precision, it is 
            assumed that the relative errors in the functions are 
            of the order of the machine precision.

        """
        x = np.asarray(x) + 0.0
        
        if self.model._parameters.paramdim != 1:
            # for now only single data sets ca be fitted
            raise ValueError("NonLinearLSQFitter can only fit one "
                             "data set at a time")

        if z is None:
            if x.shape[0] != y.shape[0]:
                raise ValueError("x and y should have the same shape")
            meas = np.asarray(y) +0.0
            farg = (meas, x)
        else:
            if x.shape != z.shape:
                raise ValueError("x, y and z should have the same shape")
            y = np.asarray(y) + 0.0
            meas = np.asarray(z) + 0.0
            farg = (meas, x, y)
 
        self.fitpars, status, dinfo, mess, ierr = optimize.leastsq(
                    self.errorfunc, self.fitpars, args=farg, Dfun=self.dfunc,
                    maxfev=maxiter, epsfcn=epsilon, full_output=True)
        
        self.fit_info.update(dinfo)
        self.fit_info['status'] = status
        self.fit_info['message'] = mess
        self.fit_info['ierr'] = ierr

class SLSQPFitter(Fitter):
    """
    Sequential Least Squares Programming optimization algorithm [6]_
    
    Supports tied and frozen parameters, as well as bounds
    
    References
    ----------
    .. [6] http://www.netlib.org/toms/733
    
    """
    def __init__(self, model):
        """
        Parameters
        ----------
        model: a fittable :class: `models.ParametricModel`
               model to fit to data
       fixed: iterable
                  a tuple of parameter names to be held fixed during fitting
        tied: dict
                keys are parameter names
                values are callable/function providing a relationship
                between parameters. Currently the callable takes a model 
                instance as an argument.
               In the example below xcen is tied to the value of xsigma 
               
               def tied(model):
                   xcen = 50*model.xsigma
                   return xcen
        
               tied ={'xcen':tie_center}
        bounds: dict
                keys: parameter names
                values:  list of length 2 giving the desired range for hte parameter
        eqcons: list
                 A list of functions of length n such that
                 eqcons[j](x0,*args) == 0.0 in a successfully optimized
                 problem.
        ineqcons: list
                   A list of functions of length n such that
                   ieqcons[j](x0,*args) >= 0.0 in a successfully optimized
                   problem.
        
        
        Raises
        ------
        ModelLinearityException
            A linear model is passed to a nonlinear fitter
            
        """
        super(SLSQPFitter, self).__init__(model)
        if self.model.linear:
            raise ModelLinearityException('Model is linear in parameters, '
                         'non-linear fitting methods should not be used.')
        
        self.fit_info = {'final_func_val': None,
                       'numiter': None,
                       'exit_mode': None,
                       'message': None
                       }
        
    def errorfunc(self, fps, *args):
        """
        Compute the sum of the squared residuals
        """
        meas = args[0]
        self.fitpars = fps
        res = self.model(*args[1:]) - meas
        return np.sum(res**2)

            
    def __call__(self, x, y , z=None, w=None, verblevel=0, 
                 maxiter=MAXITER, epsilon=EPS):
        """
        Parameters
        ----------
        x: array
           input coordinates
        y: array
           input coordinates
        z: array (optional)
           input coordinates
        w: array (optional)
           weights
        verblevel: int
                   0-silent
                   1-print summary upon completion, 
                   2-print summary after each iteration
        maxiter - maximum number of iterations
        epsilon : float
                  The step size for finite-difference derivative 
                  estimates.
        """
        x = np.asarray(x) + 0.0
        
        if self.model._parameters.paramdim != 1:
            # for now only single data sets ca be fitted
            raise ValueError("NonLinearLSQFitter can only fit "
                             "one data set at a time")

        if not z:
            if x.shape[0] != y.shape[0]:
                raise ValueError("x and y should have the same shape")
            meas = np.asarray(y) + 0.0 
            fargs = (meas, x)
        else:
            if x.shape != z.shape:
                raise ValueError("x, y and z should have the same shape")
            y = np.asarray(y) + 0.0
            meas = np.asarray(z) + 0.0
            fargs = (meas, x, y)
        p0 = self.model._parameters[:]
        bounds = [self.model.constraints.bounds[par] for par in self.model.parnames]
        self.fitpars, final_func_val, numiter, exit_mode, mess = optimize.fmin_slsqp(
            self.errorfunc, p0, args=fargs, disp=verblevel, full_output=1,
            bounds=bounds, eqcons=self.model.constraints.eqcons, 
            ieqcons=self.model.constraints.ineqcons, iter=maxiter, acc=1.E-6, 
            epsilon=EPS)
        self.fit_info['final_func_val'] = final_func_val
        self.fit_info['numiter'] = numiter
        self.fit_info['exit_mode'] = exit_mode
        self.fit_info['message']  = mess
        
class JointFitter(object):
    """
    Fit models which share a parameter
    
    For example, fit two gaussians to two data sets but keep 
    the FWHM the same.
    """
    def __init__(self, models, jointparameters, initvals):
        """
        Parameters
        ----------
        models: list
        jointparameters: list
        initvals: list
        """
        self.models = list(models)
        self.initvals = list(initvals)
        self.jointpars = jointparameters
        self._verify_input()
        for m in self.jointpars.keys():
            m.set_joint_parameters(self.jointpars[m])
        self.fitpars = self._model_to_fit_pars()
        
        # a list of model.ndim
        self.modeldims = [m.ndim for m in self.models]
        # sum all model dimensions
        self.ndim = np.sum(self.modeldims)
    
    def _model_to_fit_pars(self):
        fpars = []
        fpars.extend(self.initvals)
        for model in self.models:
            pars = model._parameters[:]
            for pname in model.joint:
                sl = model._parameters.parinfo[pname][0]
                del pars[sl]
            fpars.extend(pars)
        return fpars
        
    def errorfunc(self, fps, *args):
        """
        fps - fitted parameters
        args - args is always passed as a tuple from optimize.leastsq
        """
        lstsqargs = list(args[:])
        fitted = []
        fitpars = list(fps[:])
        numjp = len(self.initvals)
        # make a separate list of the joint fitted parameters
        jointfitpars = fitpars[:numjp]
        del fitpars[:numjp]
        
        for model in self.models:
            margs = lstsqargs[:model.ndim+1]
            del lstsqargs[:model.ndim+1]
            #separate each model separately fitted parameters
            numfp = len(model._parameters) - len(model.joint)
            mfpars = fitpars[:numfp]
            
            del fitpars[:numfp]
            #recreate the model parameters
            mpars = []
            for pname in model.parnames:
                if pname in model.joint:
                    index = model.joint.index(pname)
                    # should do this with slices in case the 
                    # parameter is not a number
                    mpars.extend([jointfitpars[index]])
                else:
                    sl = model._parameters.parinfo[pname][0]
                    plen = sl.stop - sl.start
                    mpars.extend(mfpars[:plen])
                    del mfpars[:plen]
            modelfit = model.eval(margs[:-1], mpars)
            fitted.extend(modelfit-margs[-1])
        return np.ravel(fitted)
    
    def _verify_input(self):
        assert(len(self.models)>1)
        assert(len(self.jointpars.keys()) >=2)
        for j in self.jointpars.keys():
            assert(len(self.jointpars[j]) == len(self.initvals))
      
    def __call__(self, *args):
        assert(len(args) == reduce(lambda x, y: x+1 + y+1, self.modeldims))
        self.fitpars[:], s = optimize.leastsq(self.errorfunc, self.fitpars, 
                                              args=args) 
        
        fpars = self.fitpars[:]
        numjp = len(self.initvals)
        #make a separate list of the joint fitted parameters
        jointfitpars = fpars[:numjp]
        del fpars[:numjp]
        
        for model in self.models:
            # extract each model's fitted parameters
            numfp = len(model._parameters) - len(model.joint)
            mfpars = fpars[:numfp]
            
            del fpars[:numfp]
            # recreate the model parameters
            mpars = []
            for pname in model.parnames:
                if pname in model.joint:
                    index = model.joint.index(pname)
                    # should do this with slices in case the parameter is not a number
                    mpars.extend([jointfitpars[index]])
                else:
                    sl = model._parameters.parinfo[pname][0]
                    plen = sl.stop - sl.start
                    mpars.extend(mfpars[:plen])
                    del mfpars[:plen]
            model._parameters[:] = np.array(mpars)

