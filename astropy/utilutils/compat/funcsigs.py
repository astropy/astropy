from inspect import signature, Parameter, Signature, BoundArguments

__all__ = ['BoundArguments', 'Parameter', 'Signature', 'signature']

import warnings
from ..exceptions import AstropyDeprecationWarning

warnings.warn("astropy.utils.compat.funcsigs is now deprecated - "
              "use inspect instead", AstropyDeprecationWarning)
