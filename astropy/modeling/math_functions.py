# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Define Numpy Ufuncs as Models.
"""
import numpy as np

from astropy.modeling.core import Model

trig_ufuncs = [
    "sin",
    "cos",
    "tan",
    "arcsin",
    "arccos",
    "arctan",
    "arctan2",
    "hypot",
    "sinh",
    "cosh",
    "tanh",
    "arcsinh",
    "arccosh",
    "arctanh",
    "deg2rad",
    "rad2deg",
]


math_ops = [
    "add",
    "subtract",
    "multiply",
    "logaddexp",
    "logaddexp2",
    "true_divide",
    "floor_divide",
    "negative",
    "positive",
    "power",
    "remainder",
    "fmod",
    "divmod",
    "absolute",
    "fabs",
    "rint",
    "exp",
    "exp2",
    "log",
    "log2",
    "log10",
    "expm1",
    "log1p",
    "sqrt",
    "square",
    "cbrt",
    "reciprocal",
    "divide",
    "mod",
]


supported_ufuncs = trig_ufuncs + math_ops


# These names are just aliases for other ufunc objects
# in the numpy API.  The alias name must occur later
# in the lists above.
alias_ufuncs = {
    "divide": "true_divide",
    "mod": "remainder",
}


class _NPUfuncModel(Model):
    _is_dynamic = True

    def __init__(self, **kwargs):
        super().__init__(**kwargs)


def _make_class_name(name):
    """Make a ufunc model class name from the name of the ufunc."""
    return name[0].upper() + name[1:] + "Ufunc"


def ufunc_model(name):
    """Define a Model from a Numpy ufunc name."""
    ufunc = getattr(np, name)
    nin = ufunc.nin
    nout = ufunc.nout
    if nin == 1:
        separable = True

        def evaluate(self, x):
            return self.func(x)

    else:
        separable = False

        def evaluate(self, x, y):
            return self.func(x, y)

    klass_name = _make_class_name(name)

    members = {
        "n_inputs": nin,
        "n_outputs": nout,
        "func": ufunc,
        "linear": False,
        "fittable": False,
        "_separable": separable,
        "_is_dynamic": True,
        "evaluate": evaluate,
    }

    klass = type(str(klass_name), (_NPUfuncModel,), members)
    klass.__module__ = "astropy.modeling.math_functions"
    return klass


__all__ = []

for name in supported_ufuncs:
    if name in alias_ufuncs:
        klass_name = _make_class_name(name)
        alias_klass_name = _make_class_name(alias_ufuncs[name])
        globals()[klass_name] = globals()[alias_klass_name]
        __all__.append(klass_name)
    else:
        m = ufunc_model(name)
        klass_name = m.__name__
        globals()[klass_name] = m
        __all__.append(klass_name)
