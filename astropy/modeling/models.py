# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Creates a common namespace for all pre-defined models.
"""

# pylint: disable=unused-wildcard-import, unused-import, wildcard-import

from . import math_functions as math  # noqa: F401, F403
from .core import custom_model, fix_inputs, hide_inverse  # pylint: disable=W0611 # noqa: F401, F403
from .functional_models import *  # noqa: F401, F403
from .mappings import *  # noqa: F401, F403
from .physical_models import *  # noqa: F401, F403
from .polynomial import *  # noqa: F401, F403
from .powerlaws import *  # noqa: F401, F403
from .projections import *  # noqa: F401, F403
from .rotations import *  # noqa: F401, F403
from .spline import *  # noqa: F401, F403
from .tabular import *  # noqa: F401, F403

# Attach a docstring explaining constraints to all models which support them.
# Note: add new models to this list

CONSTRAINTS_DOC = """
    Other Parameters
    ----------------
    fixed : a dict, optional
        A dictionary ``{parameter_name: boolean}`` of parameters to not be
        varied during fitting. True means the parameter is held fixed.
        Alternatively the `~astropy.modeling.Parameter.fixed`
        property of a parameter may be used.
    tied : dict, optional
        A dictionary ``{parameter_name: callable}`` of parameters which are
        linked to some other parameter. The dictionary values are callables
        providing the linking relationship.  Alternatively the
        `~astropy.modeling.Parameter.tied` property of a parameter
        may be used.
    bounds : dict, optional
        A dictionary ``{parameter_name: value}`` of lower and upper bounds of
        parameters. Keys are parameter names. Values are a list or a tuple
        of length 2 giving the desired range for the parameter.
        Alternatively, the
        `~astropy.modeling.Parameter.min` and
        `~astropy.modeling.Parameter.max` properties of a parameter
        may be used.
    eqcons : list, optional
        A list of functions of length ``n`` such that ``eqcons[j](x0,*args) ==
        0.0`` in a successfully optimized problem.
    ineqcons : list, optional
        A list of functions of length ``n`` such that ``ieqcons[j](x0,*args) >=
        0.0`` is a successfully optimized problem.
"""


MODELS_WITH_CONSTRAINTS = [
    AiryDisk2D,  # noqa: F405
    Moffat1D,  # noqa: F405
    Moffat2D,  # noqa: F405
    Box1D,  # noqa: F405
    Box2D,  # noqa: F405
    Const1D,  # noqa: F405
    Const2D,  # noqa: F405
    Ellipse2D,  # noqa: F405
    Disk2D,  # noqa: F405
    Gaussian1D,  # noqa: F405
    Gaussian2D,  # noqa: F405
    Linear1D,  # noqa: F405
    Lorentz1D,  # noqa: F405
    RickerWavelet1D,  # noqa: F405
    RickerWavelet2D,  # noqa: F405
    PowerLaw1D,  # noqa: F405
    Sersic1D,  # noqa: F405
    Sersic2D,  # noqa: F405
    Sine1D,  # noqa: F405
    Cosine1D,  # noqa: F405
    Tangent1D,  # noqa: F405
    ArcSine1D,  # noqa: F405
    ArcCosine1D,  # noqa: F405
    ArcTangent1D,  # noqa: F405
    Trapezoid1D,  # noqa: F405
    TrapezoidDisk2D,  # noqa: F405
    Chebyshev1D,  # noqa: F405
    Chebyshev2D,  # noqa: F405
    Hermite1D,  # noqa: F405
    Hermite2D,  # noqa: F405
    Legendre2D,  # noqa: F405
    Legendre1D,  # noqa: F405
    Polynomial1D,  # noqa: F405
    Polynomial2D,  # noqa: F405
    Voigt1D,  # noqa: F405
    KingProjectedAnalytic1D,  # noqa: F405
    NFW,  # noqa: F405
]


for item in MODELS_WITH_CONSTRAINTS:
    if isinstance(item.__doc__, str):
        item.__doc__ += CONSTRAINTS_DOC
