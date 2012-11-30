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
    def __init__(self, model, fixed=None, tied=None, bounds=None, 
                           eqcons=None, ineqcons=None):
        __metaclass__ = abc.ABCMeta
        
        self.model = model
        self.constraints = Constraints(self, fixed=fixed, tied=tied, bounds=bounds,
                                                        eqcons=eqcons, ineqcons=ineqcons)
        if self.constraints.pmask:
            self._fitpars = self.constraints.fitpars[:]
        else:
            self._fitpars = self.model._parameters[:]
        
   
    @property
    def fitpars(self):
        if self.constraints.pmask:
            return self.constraints.fitpars
        else:
            return self.model._parameters
        
    @fitpars.setter
    def fitpars(self, fps):
        if self.constraints.pmask:
            self.constraints.fitpars = fps
            self._fitpars[:] = self.constraints.fitpars
        else:
            self._fitpars[:] = fps
            self.model._parameters._update(fps)
            
    @property
    def modelpars(self):
        """
        modelpars is set through constraints.fitpars 
        or through model._parameters
        """
        if self.constraints.pmask:
            return self.constraints.modelpars
        else:
            return self.model._parameters
        
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
    def __init__(self, model, fixed=None):
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
        super(LinearLSQFitter, self).__init__(model, fixed=fixed)
        
    def _deriv_with_constraints(self, x, y=None):
        if y is None:
            d = self.model.deriv(x)
        else:
            d = self.model.deriv(x, y)
        ind = []
        for item in self.constraints.fixed:
                ind.append(self.model.parnames.index(item))
        ind.sort()
        for i in range(len(ind)):
            d = np.delete(d, ind[i]-i, 1)
        return d
        
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
            #xnew = self.model._set_domain_1d(x)
            if not self.model.domain:
                self.model.domain = [x.min(), x.max()]
            if not self.model.window:
                self.model.window = [-1, 1] #self.model.domain[:]
            xnew = pmapdomain(x, self.model.domain, self.model.window)
            if self.constraints and self.constraints.fixed:
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
            
            if self.constraints and self.constraints.fixed:
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
        #self.model._parameters = lacoef.flatten()[:]
        self.fitpars = lacoef.flatten()[:]
        #self.model._parameters._update(lacoef.flatten())
        
class NonLinearLSQFitter(Fitter):
    """
    A wrapper around scipy.optimize.leastsq
    
    Supports tied and frozen parameters.
    
    """
    def __init__(self, model, fixed=None, tied=None):
        """
        Parameters
        ----------
        model: a fittable :class: `models.ParametricModel`
               model to fit to data
        ixed: iterable
                  a tuple of parameter names to be held fixed during fitting
        tied: dict
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
                        
        super(NonLinearLSQFitter, self).__init__(model, fixed=fixed, tied=tied)
        if self.model.linear:
            raise ModelLinearityException('Model is linear in parameters, '
                            'non-linear fitting methods should not be used.')
               
    def errorfunc(self, fps, *args):
        self.fitpars = fps
        meas = args[0]
        return np.ravel(self.model(*args[1:]) - meas)
    
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
            self.errorfunc, self.fitpars, args=farg, Dfun=self.model.fjac,
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
    def __init__(self, model, fixed=None, tied=None, bounds=None,
                            eqcons=None, ineqcons=None):
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
        super(SLSQPFitter, self).__init__(model, fixed=fixed, tied=tied, bounds=bounds, eqcons=eqcons, 
                        ineqcons=ineqcons)
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
    
    def _validate_constraints(self):
        """
        default values for the constraints arguments
        """
        if self.constraints:
            if not self.constraints.eqcons:
                self.constraints.eqcons = []
            if not self.constraints.ineqcons:
                self.constraints.ineqcons = []
        else:
            self.constraints.bounds = []
            self.constraints.eqcons = []
            self.constraints.ineqcons = []
            
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
        self._validate_constraints()
        
        self.fitpars, final_func_val, numiter, exit_mode, mess = optimize.fmin_slsqp(
            self.errorfunc, p0, args=fargs, disp=verblevel, full_output=1,
            bounds=self.constraints._range, eqcons=self.constraints.eqcons, 
            ieqcons=self.constraints.ineqcons, iter=maxiter, acc=1.E-6, 
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

class Constraints(object):
    """
    Fitting constraints

    """
    fitters = {'NonLinearLSQFitter': ['fixed', 'tied'],
                'SLSQPFitter': ['bounds', 'eqcons', 'ineqcons', 'fixed', 'tied'],
                'LinearLSQFitter': ['fixed'],
                }
    def __init__(self, fitter, fixed=(), tied={}, bounds={},
                            eqcons=[], ineqcons=[]):
        """
        Parameters
        ----------
        fitter: object which supports fitting
        fixed: iterable
                  a tuple of parameter names to be held fixed during fitting
        tied: dict
                keys are parameter names
                values are callable/function providing a relationship
                between parameters. Currently the callable takes a model 
                instance as an argument.
               In the example below xcen is tied to the value of xsigma
               
               def tie_center(model):
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
        """
        self.fitter = fitter
        
        self._fixed = fixed
        self._tied = tied
        self._bounds = bounds
        self._eqcons = eqcons
        self._ineqcons = ineqcons
        self._validate_fitter()
        
        self.pmask = {}
        #it is unclear what it means to have more than one parameter 
        #constraint on the same parameter 
        # currently the behaviour is that it picks one of them
        if self._fixed:
            self.pmask.update({}.fromkeys(self._fixed, False))
            if self._tied:
                self.pmask.update(self._tied)
        elif self._tied:
            self.pmask = self._tied.copy()
        # bounds internally are converted to self._range and this is used as an
        # argument to the fiting routine
        self._range = []
        if self.pmask:
            self._fitpars = self._model_to_fit_pars()
        else:
            self._fitpars = self.fitter.model._parameters[:]
        
    def __repr__(self):
        fmt = "Constraints(%s, fixed=%s, tied=%s, bounds=%s, eqcons=%s, ineqcons=%s)" % \
                        (self.fitter.__class__.__name__, repr(self.fixed), repr(self.tied), 
                        repr(self.bounds), repr(self.eqcons), repr(self.ineqcons))
        return fmt
    
    def __str__(self):
        fmt = "%s \nConstraints:\n" % self.fitter.__class__.__name__
        if self.fixed:
            fmt += 'Fixed parameters: %s\n' % str(self.fixed)
        if self.tied:
            fmt += 'Tied parameters: %s\n' % str(self.tied)
        if self.bounds:
            fmt += 'Bound parameters: %s\n' % str(self.bould)
        if self.eqcons:
            fmt += "Equality constraints: %s\n" % str(self.eqcons)
        if self.ineqcons:
            fmt += "Inequality constraints: %s\n" % str(self.ineqcons)
        return fmt
    
    def _validate_fitter(self):
        fname = self.fitter.__class__.__name__
        try:
            c = self.fitters[fname]
        except KeyError:
            print("%s does not support fitting with constraints" % fname)
            raise
        if self._fixed and 'fixed' not in c:
            raise ValueError("%s cannot handle fixed parameter constraints "\
                                            % fname)
        if self._tied and 'tied' not in c:
            raise ValueError("%s cannot handle tied parameter constraints "\
                                            % fname)
        if self._bounds and 'bounds' not in c:            
            raise ValueError("%s cannot handle bound parameter constraints"
                             % fname)
        if self._eqcons and 'eqcons' not in c:
            raise ValueError("%s cannot handle equality constraints but "
                             "eqcons given" % fname)
        if self._ineqcons and 'ineqcons' not in c:
            raise ValueError("%s cannot handle inequality constraints but "
                             "ineqcons given" % fname)
    
    @property
    def fixed(self):
        return self._fixed
    
    @fixed.setter
    def fixed(self, fixedparlist):
        self._validate_fitter()
        if self._fixed:
            for key in self._fixed:
                self.pmask.pop(key)
        self._fixed = tuple(fixedparlist)
        self.pmask.update({}.fromkeys(fixedparlist, False))
        
    @property
    def tied(self):
        return self._tied
    
    @tied.setter
    def tied(self, funcdict):
        self._validate_fitter()
        self._tied = funcdict
        self.pmask.update(funcdict)
      
    @property
    def eqcons(self):
        return self._eqcons
    
    @eqcons.setter
    def eqcons(self, eqlist):
        self._validate_fitter()
        self._eqcons = eqlist
        
    @property
    def ineqcons(self):
        return self._ineqcons
    
    @ineqcons.setter
    def ineqcons(self, ineqlist):
        self._validate_fitter()
        self._ineqcons = ineqlist
        
    @property
    def bounds(self):
        return self._bounds
    
    @bounds.setter
    def bounds(self, b):
        self._validate_fitter()
        bl = []
        assert(isinstance(b, dict))
        defaultval = (-1.E12, 1.E12)
        for pname in self.fitter.model.parnames:
            if pname in b.keys():
                bl.append(tuple(b[pname]))
            else:
                bl.append(defaultval)
        self._range = bl[:]

    def _add_range(self, mask):
        for item in mask.items():
            if operator.isSequenceType(item[1]) and len(item[1]) == 2:
                self.boundary.update({item[0]: item[1]})
                self.fitter.pmask.pop(item[0])
                
    @property
    def fitpars(self):
        """
        Return the fitted parameters
        """
        return self._fitpars
    
    @fitpars.setter
    def fitpars(self, fp):
        """
        Used by the fitting routine to set fitpars and modelpars
        """
        self._fitpars[:] = fp
        self.modelpars = fp
        
    @property
    def modelpars(self):
        """
        Return the model parameters
        """
        return self.fitter.model._parameters
    
    @modelpars.setter
    def modelpars(self, fp):
        """
        Sets model parameters.
        
        Takes into account any constant or tied parameters
        and inserts them into the list of fitted parameters.
        """
        fitpars = list(fp[:])
        mpars = []
        for par in self.fitter.model.parnames:
            if par in self.pmask.keys():
                if self.pmask[par] == False:
                    mpars.extend(getattr(self.fitter.model, par))
                else:
                    mpars.extend([self.pmask[par](self.fitter.model)])
            else:
                sl = self.fitter.model._parameters.parinfo[par][0]
                plen = sl.stop - sl.start
                mpars.extend(fitpars[:plen])
                del fitpars[:plen]
        self.fitter.model._parameters._update(np.array(mpars))
            
    def _model_to_fit_pars(self):
        """
        Create a set of parameters to be fitted.
        These may be a subset of the model parameters, if some
        of them are help constant or tied.
        """
        pars = self.fitter.model._parameters[:]
        for item in self.pmask.items():
            if isinstance(item[1], bool) or operator.isCallable(item[1]):
                sl = self.fitter.model._parameters.parinfo[item[0]][0]
                del pars[sl]
        return pars
    
