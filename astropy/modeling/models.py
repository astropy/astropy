# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Creates a common namespace for all pre-defined models.
"""

from . import math_functions as math
from .core import custom_model, fix_inputs, hide_inverse
from .functional_models import (
    Scale,
    Shift,
    Multiply,
    Planar2D,
    Exponential1D,
    Logarithmic1D,
    RedshiftScaleFactor,
    Ring2D,
    AiryDisk2D,
    Moffat1D,
    Moffat2D,
    Box1D,
    Box2D,
    Const1D,
    Const2D,
    Ellipse2D,
    Disk2D,
    Gaussian1D,
    Gaussian2D,
    Linear1D,
    Lorentz1D,
    RickerWavelet1D,
    RickerWavelet2D,
    Sersic1D,
    Sersic2D,
    Sine1D,
    Cosine1D,
    Tangent1D,
    ArcSine1D,
    ArcCosine1D,
    ArcTangent1D,
    Trapezoid1D,
    TrapezoidDisk2D,
    Voigt1D,
    KingProjectedAnalytic1D,
)
from .mappings import (
    Identity,
    Mapping,
    UnitsMapping
)
from .physical_models import (
    BlackBody,
    NFW,
    Drude1D,
    Plummer1D
)
from .polynomial import (
    SIP,
    Chebyshev1D,
    Chebyshev2D,
    Hermite1D,
    Hermite2D,
    Legendre2D,
    Legendre1D,
    Polynomial1D,
    Polynomial2D
)
from .powerlaws import (
    PowerLaw1D,
    BrokenPowerLaw1D,
    SmoothlyBrokenPowerLaw1D, 
    ExponentialCutoffPowerLaw1D, 
    LogParabola1D, 
    Schechter1D, 
)
from .projections import (
    Pix2Sky_TAN
)
from .rotations import (
    RotateCelestial2Native, 
    RotateNative2Celestial, 
    Rotation2D, 
    EulerAngleRotation, 
    RotationSequence3D, 
    SphericalRotationSequence
)
from .spline import (
    Spline1D,
    SplineInterpolateFitter, 
    SplineSmoothingFitter, 
    SplineExactKnotsFitter, 
    SplineSplrepFitter, 
)
from .tabular import (
    Tabular1D,
    Tabular2D,
    tabular_model
)

__all__ = [
    'math',
    'Plummer1D',
    'custom_model',
    'fix_inputs',
    'hide_inverse',
    'Scale',
    'Shift',
    'Multiply',
    'Planar2D',
    'Exponential1D',
    'Logarithmic1D',
    'RedshiftScaleFactor',
    'Ring2D',
    'AiryDisk2D',
    'Moffat1D',
    'Moffat2D',
    'Box1D',
    'Box2D',
    'Const1D',
    'Const2D',
    'Ellipse2D',
    'Disk2D',
    'Gaussian1D',
    'Gaussian2D',
    'Linear1D',
    'Lorentz1D',
    'RickerWavelet1D',
    'RickerWavelet2D',
    'Sersic1D',
    'Sersic2D',
    'Sine1D',
    'Cosine1D',
    'Tangent1D',
    'ArcSine1D',
    'ArcCosine1D',
    'ArcTangent1D',
    'Trapezoid1D',
    'TrapezoidDisk2D',
    'Voigt1D',
    'KingProjectedAnalytic1D',
    'Identity',
    'Mapping',
    'UnitsMapping',
    'BlackBody',
    'NFW',
    'Drude1D',
    'SIP',
    'OrthoPolynomialBase',
    'Chebyshev1D',
    'Chebyshev2D',
    'Hermite1D',
    'Hermite2D',
    'Legendre2D',
    'Legendre1D',
    'Polynomial1D',
    'Polynomial2D',
    'PowerLaw1D',
    'BrokenPowerLaw1D',
    'SmoothlyBrokenPowerLaw1D', 
    'ExponentialCutoffPowerLaw1D', 
    'LogParabola1D', 
    'Schechter1D', 
    'Pix2Sky_TAN',
    'Rotation2D',
    'RotateNative2Celestial',
    'RotateCelestial2Native',
    'EulerAngleRotation',
    'RotationSequence3D', 
    'SphericalRotationSequence',
    'Spline1D',
    'SplineInterpolateFitter', 
    'SplineSmoothingFitter', 
    'SplineExactKnotsFitter', 
    'SplineSplrepFitter', 
    'Tabular1D',
    'Tabular2D',
    'tabular_model',
]

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
    AiryDisk2D,
    Moffat1D,
    Moffat2D,
    Box1D,
    Box2D,
    Const1D,
    Const2D,
    Ellipse2D,
    Disk2D,
    Gaussian1D,
    Gaussian2D,
    Linear1D,
    Lorentz1D,
    RickerWavelet1D,
    RickerWavelet2D,
    PowerLaw1D,
    Sersic1D,
    Sersic2D,
    Sine1D,
    Cosine1D,
    Tangent1D,
    ArcSine1D,
    ArcCosine1D,
    ArcTangent1D,
    Trapezoid1D,
    TrapezoidDisk2D,
    Chebyshev1D,
    Chebyshev2D,
    Hermite1D,
    Hermite2D,
    Legendre2D,
    Legendre1D,
    Polynomial1D,
    Polynomial2D,
    Voigt1D,
    KingProjectedAnalytic1D,
    NFW,
]


for item in MODELS_WITH_CONSTRAINTS:
    if isinstance(item.__doc__, str):
        item.__doc__ += CONSTRAINTS_DOC