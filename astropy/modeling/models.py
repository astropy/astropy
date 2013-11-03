# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Creates a common namespace for all pre-defined models.
"""

from .projections import *
from .rotations import *
from .polynomial import *
from .functional_models import *
from .powerlaws import *

"""
Attach a docstring explaining constraints to all models which support them.

Note: add new models to this list
"""

CONSTRAINTS_DOC = """
    fixed: a dict
        a dictionary ``{parameter_name: boolean}`` of parameters to not be
        varied during fitting. True means the parameter is held fixed.
        Alternatively the `~astropy.modeling.parameters.Parameter.fixed`
        property of a parameter may be used.
    tied: dict
        a dictionary ``{parameter_name: callable}`` of parameters which are
        linked to some other parameter. The dictionary values are callables
        providing the linking relationship.  Alternatively the
        `~astropy.modeling.parameters.Parameter.tied` property of a parameter
        may be used.
    bounds: dict
        a dictionary ``{parameter_name: boolean}`` of lower and upper bounds of
        parameters. Keys  are parameter names. Values  are a list of length 2
        giving the desired range for the parameter.  Alternatively the
        `~astropy.modeling.parameters.Parameter.min` and
        `~astropy.modeling.parameters.Parameter.max` properties of a parameter
        may be used.
    eqcons: list
        A list of functions of length ``n`` such that ``eqcons[j](x0,*args) ==
        0.0`` in a successfully optimized problem.
    ineqcons : list
        A list of functions of length ``n`` such that ``ieqcons[j](x0,*args) >=
        0.0`` is a successfully optimized problem.
"""


MODELS_WITH_CONSTRAINTS = [
    AiryDisk2D, Beta1D, Beta2D, Box1D, Box2D,
    Const1D, Const2D, Disk2D, Gaussian1D, Gaussian2D,
    Linear1D, Lorentz1D, MexicanHat1D, MexicanHat2D,
    PowerLaw1D, Sine1D, Trapezoid1D, TrapezoidDisk2D,
    Chebyshev1D, Chebyshev2D, Legendre2D, Legendre1D,
    Polynomial1D, Polynomial2D
]


for item in MODELS_WITH_CONSTRAINTS:
    item.__doc__ += CONSTRAINTS_DOC
