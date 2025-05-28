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

    if kwargs.get("retstep"):
        offsets, step = np.linspace(
            np.zeros(tstart.shape), np.ones(tstop.shape), *args, **kwargs
        )
        tdelta = tstop - tstart
        return tstart + tdelta * offsets, tdelta * step
    else:
        offsets = np.linspace(
            np.zeros(tstart.shape), np.ones(tstop.shape), *args, **kwargs
        )
        return tstart + (tstop - tstart) * offsets


@custom_functions
def zeros_like(a, dtype=None, order="K", subok=True, shape=None, *, device=None):
    """Create a new Time object set to J2000 with the properties of a.

    The other parameters allow to override what is used to create
    ``jd1`` and ``jd2`` (but one cannot override ``dtype``).
    """
    if dtype is not None:
        raise ValueError("Cannot set dtype for Time creation.")
    if shape is None:
        shape = a.shape
    jd2000 = 2451544.5  # Arbitrary JD value J2000.0 that will work with ERFA
    jd1 = np.full_like(a.jd1, jd2000, shape=shape, order=order, device=device)
    jd2 = np.zeros_like(a.jd2, shape=shape, order=order, device=device)
    tmp = a.__class__(jd1, jd2, format="jd", scale=a.scale, copy=False)
    # Change format to the correct one and add back other attributes.
    # This will break if location was an array that cannot be broadcast to
    # the new shape, but that seems reasonable; it is unclear what the user
    # wants in that case.
    return a.__class__(
        tmp,
        format=a.format,
        location=a.location,
        precision=a.precision,
        in_subfmt=a.in_subfmt,
        out_subfmt=a.out_subfmt,
        copy=False,
    )
