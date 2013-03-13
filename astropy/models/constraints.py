from __future__ import division, print_function
import numpy as np
class Constraints(object):
    """
    Fitting constraints

    Should not be instanciated by users directly. 
    Instead constraints should be passed to an instance of Model.
    
    """
    
    def __init__(self, model, fixed={}, tied={}, bounds={},
                            eqcons=[], ineqcons=[]):
        """
        Parameters
        ----------
        fitter : object which supports fitting
        fixed : dict
            {parameter_name: True| False}
            Specify parameters which should be kept fixed during fitting.
            The keys of the dictionary are parameter names. If the value is
            True the parameter is not varied. The default is False.
        tied : dict
            Specify a relationship between parameters.
            The keys of the dictionary are parameter names.
            Values are a callable providing a relationship
            between parameters. Currently the callable takes a model 
            instance as an argument.
        bounds : dict
            Specify bounds on parameters.
            Keys  are parameter names.
            Values  are a list of length 2 giving the desired range for the parameter.
        eqcons : list
            A list of functions of length n such that
            eqcons[j](x0,*args) == 0.0 in a successfully optimized
            problem.
        ineqcons : list
            A list of functions of length n such that
            ieqcons[j](x0,*args) >= 0.0 is a successfully optimized
            problem.
            
        Examples
        ---------
        >>> def tie_center(model):
        ...         xcen = 50 * model.xsigma
        ...         return xcen
        >>> tied_parameters  ={'xcen': tie_center}
        
        Specify that 'xcen' is a tied parameter in one of two ways:
        
        >>> g1=models.Gauss1DModel(amplitude=10, xcen=5, xsigma=.3, tied=tied)

        or

        >>> g1=models.Gauss1DModel(amplitude=10, xcen=5, xsigma=.3)
        >>> g1.xcen.tied
        False
        >>> g1.xcen.tied = tie_center
        >>> g1.xcen.tied
        <function tie_center at 0x395ab0>

        Fixed parameters:
        
        >>> g1=models.Gauss1DModel(amplitude=10, xcen=5, xsigma=.3, fixed={'xsigma':True})
        >>> g1.xsigma.fixed
        True

        or
        
        >>> g1=models.Gauss1DModel(amplitude=10, xcen=5, xsigma=.3)
        >>> g1.xsigma.fixed
        False
        >>> g1.xsigma.fixed=True
        >>> g1.xsigma.fixed
        True

        """
        self.model = model
        
        self._fixed = fixed
        self._tied = tied
        self._bounds = bounds
        self.set_range(bounds)
        self._eqcons = eqcons
        self._ineqcons = ineqcons
        if any(self._fixed.values()) or any(self._tied.values()):
            self._fitpars = self._model_to_fit_pars()
        else:
            self._fitpars = self.model._parameters[:]
        
    def __str__(self):
        fmt = ""
        if any(self._fixed.values()):
            fixedstr = [par for par in self._fixed if self._fixed[par]]
            fmt += "fixed={0}".format(fixedstr)
        if any(self._tied.values()):
            tiedstr = [par+": "+self._tied[par].__name__+"()" for
                       par in self._tied if self._tied[par]]
            fmt += "tied={0}".format(tiedstr)
        if not all([(-np.inf, np.inf)==b for b in self._bounds.values()]):
            boundsstr = [par+":"+str(self._bounds[par]) for
                         par in self._bounds if
                         self._bounds[par]!= (-np.inf, np.inf)]
            fmt += "bounds={0}".format(boundsstr)
        if self._eqcons:
            fmt += "eqcons={0}".format(self._eqcons)
        if self._ineqcons:
            fmt += "ineqcons={0}".format(self._ineqcons)
        name = "Constraints({0}, ".format(self.model.__class__.__name__)
        if fmt:
            fmt = name + fmt +")"
        return fmt
    
    def __repr__(self):
        fmt = "<Constraints({0}, ".format(self.model.__class__.__name__)
        if self._fixed:
            fmt += 'fixed={0}, '.format(self._fixed)
        if self._tied:
            fmt += 'tied={0}, '.format(self._tied)
        if self._bounds:
            fmt += 'bounds={0}, '.format(self._bounds)
        if self._eqcons:
            fmt += "eqcons={0}, ".format(self._eqcons)
        if self._ineqcons:
            fmt += "ineqcons={0}".format(self._ineqcons)
        fmt += ")>"
        return fmt
    
    
    @property
    def fixed(self):
        return self._fixed
    
    @property
    def tied(self):
        return self._tied
      
    @property
    def eqcons(self):
        return self._eqcons
        
    @property
    def ineqcons(self):
        return self._ineqcons
        
    @property
    def bounds(self):
        return self._bounds

    def set_range(self, b):
        for key in b:
            b[key] = tuple(b[key])
        self._bounds.update(b)
        
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
        return self.model._parameters
    
    @modelpars.setter
    def modelpars(self, fp):
        """
        Sets model parameters.
        
        Takes into account any constant or tied parameters
        and inserts them into the list of fitted parameters.
        """
        fitpars = list(fp[:])
        mpars = []
        for par in self.model.parnames:
            if self.fixed[par]:
                mpars.extend(getattr(self.model, par))
            elif self.tied[par]:
                mpars.extend([self.tied[par](self.model)])
            else:
                sl = self.model._parameters.parinfo[par][0]
                plen = sl.stop - sl.start
                mpars.extend(fitpars[:plen])
                del fitpars[:plen]
        self.model.parameters = mpars
            
    def _model_to_fit_pars(self):
        """
        Create a set of parameters to be fitted.
        These may be a subset of the model parameters, if some
        of them are held constant or tied.
        """
        pars = self.model._parameters[:]
        for item in self.model.parnames:
            if self._fixed[item] or self.tied[item]:
                sl = self.model._parameters.parinfo[item][0]
                del pars[sl]
        return pars
    
    def _update(self):
        if any(self._fixed.values()) or any(self._tied.values()):
            self._fitpars = self._model_to_fit_pars()
        else:
            self._fitpars = self.model._parameters[:]
        
