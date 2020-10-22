# Licensed under a 3-clause BSD style license - see LICENSE.rst
# pylint: disable=invalid-name

"""
Optimization algorithms used in `~astropy.modeling.fitting`.
"""

import numpy as np

from astropy.modeling.optimizers import Optimization

try:
    import emcee
except ImportError:
    HAS_EMCEE = False
else:
    HAS_EMCEE = True


all = ["EmceeSampler"]


class EmceeSampler(Optimization):
    """
    Interface to emcee sampler.
    """

    supported_constraints = ["bounds", "fixed", "tied"]

    def __init__(self):
        super().__init__(emcee)
        self.fit_info = {"perparams": None, "samples": None, "sampler": None}

    @staticmethod
    def _get_best_fit_params(sampler):
        """
        Determine the best fit parameters given an emcee sampler object
        """
        # very likely a faster way
        max_lnp = -1e6
        nwalkers, nsteps = sampler.lnprobability.shape
        for k in range(nwalkers):
            tmax_lnp = np.nanmax(sampler.lnprobability[k])
            if tmax_lnp > max_lnp:
                max_lnp = tmax_lnp
                (indxs,) = np.where(sampler.lnprobability[k] == tmax_lnp)
                fit_params_best = sampler.chain[k, indxs[0], :]

        return fit_params_best

    def __call__(self, objfunc, initval, fargs, nsteps, save_samples=None, **kwargs):
        """
        Run the sampler.

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
        # optresult = self.opt_method(objfunc, initval, args=fargs)
        # fitparams = optresult['x']

        ndim = len(initval)
        nwalkers = 2 * ndim
        pos = initval + 1e-4 * np.random.randn(nwalkers, ndim)

        # ensure all the walkers start within the bounds
        model = fargs[0]
        for cp in pos:
            k = 0
            for cname in model.param_names:
                if not model.fixed[cname]:
                    if model.bounds[cname][0] is not None:
                        if cp[k] < model.bounds[cname][0]:
                            cp[k] = model.bounds[cname][0]
                    if model.bounds[cname][1] is not None:
                        if cp[k] > model.bounds[cname][1]:
                            cp[k] = model.bounds[cname][1]
                    # only non-fixed parameters are in initval
                    # so only increment k when non-fixed
                    k += 1

        # Set up the backend
        save_backend = None
        if save_samples:
            # Don't forget to clear it in case the file already exists
            save_backend = emcee.backends.HDFBackend(save_samples)
            save_backend.reset(nwalkers, ndim)

        sampler = self.opt_method.EnsembleSampler(
            nwalkers, ndim, objfunc, backend=save_backend, args=fargs
        )
        sampler.run_mcmc(pos, nsteps, progress=True)
        samples = sampler.get_chain()

        fitparams = self._get_best_fit_params(sampler)
        self.fit_info["sampler"] = sampler
        self.fit_info["samples"] = samples

        return fitparams, self.fit_info
