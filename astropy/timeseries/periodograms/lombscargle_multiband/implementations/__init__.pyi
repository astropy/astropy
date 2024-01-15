# Licensed under a 3-clause BSD style license - see LICENSE.rst
from .main import (
    available_methods as available_methods,
    lombscargle_multiband as lombscargle_multiband,
)
from .mbfast_impl import lombscargle_mbfast as lombscargle_mbfast
from .mbflex_impl import lombscargle_mbflex as lombscargle_mbflex
from . import (
    main as main,
    mbflex_impl as mbflex_impl,
    mle as mle,
    mbfast_impl as mbfast_impl,
)
