"""
Helpers for overriding numpy functions in
`~astropy.time.Time.__array_function__`.
"""
import numpy as np

from astropy.units.quantity_helper.function_helpers import FunctionAssigner

# TODO: Fill this in with functions that don't make sense for times
UNSUPPORTED_FUNCTIONS = {}
# Functions that return the final result of the numpy function
CUSTOM_FUNCTIONS = {}

custom_functions = FunctionAssigner(CUSTOM_FUNCTIONS)


@custom_functions(helps={np.linspace})
def linspace(tstart, tstop, *args, **kwargs):
    from astropy.time import Time
    if isinstance(tstart, Time):
        if not isinstance(tstop, Time):
            return NotImplemented

    if kwargs.get('retstep'):
        offsets, step = np.linspace(np.zeros(tstart.shape), np.ones(tstop.shape), *args, **kwargs)
        tdelta = tstop - tstart
        return tstart + tdelta * offsets, tdelta * step
    else:
        offsets = np.linspace(np.zeros(tstart.shape), np.ones(tstop.shape), *args, **kwargs)
        return tstart + (tstop - tstart) * offsets
