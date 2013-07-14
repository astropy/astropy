import numpy as np

def _is_unity(value):
    x = value.decompose()
    return (len(x.bases) == 0 and x.scale == 1.0)

UFUNC_HELPERS = {}

# In this file, we implement the logic that determines for a given ufunc and
# input what the required scaling of the input and what the output unit will be.

# The functions below take a single argument, which is the quantity upon which
# the ufunc is being used. The output of the function should be two values: the
# scale by which the input needs to be multiplied before being passed to the
# ufunc, and the unit the output will be in.

# ufuncs that return a value with the same unit as the input

helper_invariant = lambda f, unit: (1., unit)

UFUNC_HELPERS[np.absolute] = helper_invariant
UFUNC_HELPERS[np.fabs] = helper_invariant
UFUNC_HELPERS[np.conj] = helper_invariant
UFUNC_HELPERS[np.conjugate] = helper_invariant
UFUNC_HELPERS[np.negative] = helper_invariant
UFUNC_HELPERS[np.spacing] = helper_invariant
UFUNC_HELPERS[np.rint] = helper_invariant
UFUNC_HELPERS[np.floor] = helper_invariant
UFUNC_HELPERS[np.ceil] = helper_invariant
UFUNC_HELPERS[np.trunc] = helper_invariant

# ufuncs handled as special cases

UFUNC_HELPERS[np.sqrt] = lambda f, unit: (1., unit ** 0.5)
UFUNC_HELPERS[np.square] = lambda f, unit: (1., unit ** 2)
UFUNC_HELPERS[np.reciprocal] = lambda f, unit: (1., unit ** -1)

# ufuncs that require dimensionless input and and give dimensionless output

def helper_dimensionless_to_dimensionless(f, unit):
    from . import dimensionless_unscaled
    try:
        scale = unit.to(dimensionless_unscaled)
    except:
        raise TypeError("Can only apply '{0}' function to "
                        "dimensionless quantities"
                        .format(f.__name__))
    return scale, dimensionless_unscaled
    
UFUNC_HELPERS[np.exp] = helper_dimensionless_to_dimensionless
UFUNC_HELPERS[np.expm1] = helper_dimensionless_to_dimensionless
UFUNC_HELPERS[np.exp2] = helper_dimensionless_to_dimensionless
UFUNC_HELPERS[np.log] = helper_dimensionless_to_dimensionless
UFUNC_HELPERS[np.log10] = helper_dimensionless_to_dimensionless
UFUNC_HELPERS[np.log2] = helper_dimensionless_to_dimensionless
UFUNC_HELPERS[np.log1p] = helper_dimensionless_to_dimensionless

# ufuncs that require dimensionless input and give output in radians

def helper_dimensionless_to_radian(f, unit):
    from .si import radian
    from . import dimensionless_unscaled
    try:
        scale = unit.to(dimensionless_unscaled)
    except:
        raise TypeError("Can only apply '{0}' function to "
                        "dimensionless quantities"
                        .format(f.__name__))
    return scale, radian
    
UFUNC_HELPERS[np.arccos] = helper_dimensionless_to_radian
UFUNC_HELPERS[np.arcsin] = helper_dimensionless_to_radian
UFUNC_HELPERS[np.arctan] = helper_dimensionless_to_radian
UFUNC_HELPERS[np.arccosh] = helper_dimensionless_to_radian
UFUNC_HELPERS[np.arcsinh] = helper_dimensionless_to_radian
UFUNC_HELPERS[np.arctanh] = helper_dimensionless_to_radian

# ufuncs that require input in degrees and give output in radians

def helper_degree_to_radian(f, unit):
    from .si import degree, radian
    try:
        scale = unit.to(degree)
    except:
        raise TypeError("Can only apply '{0}' function to "
                        "quantities with angle units"
                        .format(f.__name__))
    return scale, radian
        
UFUNC_HELPERS[np.radians] = helper_degree_to_radian
UFUNC_HELPERS[np.deg2rad] = helper_degree_to_radian

# ufuncs that require input in radians and give output in degrees

def helper_radian_to_degree(f, unit):
    from .si import degree, radian
    try:
        scale = unit.to(radian)
    except:
        raise TypeError("Can only apply '{0}' function to "
                        "quantities with angle units"
                        .format(f.__name__))
    return scale, degree
        
UFUNC_HELPERS[np.degrees] = helper_degree_to_radian
UFUNC_HELPERS[np.rad2deg] = helper_degree_to_radian

# ufuncs that require input in radians and give dimensionless output

def helper_radian_to_dimensionless(f, unit):
    from .si import radian
    from . import dimensionless_unscaled
    try:
        scale = unit.to(radian)
    except:
        raise TypeError("Can only apply '{0}' function to "
                        "quantities with angle units"
                        .format(f.__name__))
    return scale, dimensionless_unscaled
        
UFUNC_HELPERS[np.cos] = helper_radian_to_dimensionless
UFUNC_HELPERS[np.sin] = helper_radian_to_dimensionless
UFUNC_HELPERS[np.tan] = helper_radian_to_dimensionless
UFUNC_HELPERS[np.cosh] = helper_radian_to_dimensionless
UFUNC_HELPERS[np.sinh] = helper_radian_to_dimensionless
UFUNC_HELPERS[np.tanh] = helper_radian_to_dimensionless

# ufuncs that require dimensionless_unscaled input and return non-quantities

def helper_dimensionless_to_none(f, unit):
    from . import dimensionless_unscaled
    if not _is_unity(unit):
        raise TypeError("Can only apply '{0}' function to "
                        "unscaled dimensionless quantities"
                        .format(f.__name__))
    return 1., None
        
UFUNC_HELPERS[np.modf] = helper_dimensionless_to_none
UFUNC_HELPERS[np.frexp] = helper_dimensionless_to_none
