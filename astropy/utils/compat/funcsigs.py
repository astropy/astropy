import warnings
from inspect import Parameter, Signature, BoundArguments, signature

from ..exceptions import AstropyDeprecationWarning

__all__ = ['BoundArguments', 'Parameter', 'Signature', 'signature']


warnings.warn("astropy.utils.compat.funcsigs is now deprecated - "
              "use inspect instead", AstropyDeprecationWarning)
