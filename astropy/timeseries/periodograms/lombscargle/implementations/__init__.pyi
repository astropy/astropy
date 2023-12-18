# Licensed under a 3-clause BSD style license - see LICENSE.rst
from .chi2_impl import lombscargle_chi2 as lombscargle_chi2
from .fast_impl import lombscargle_fast as lombscargle_fast
from .fastchi2_impl import lombscargle_fastchi2 as lombscargle_fastchi2
from .main import (
    available_methods as available_methods,
    lombscargle as lombscargle,
)
from .scipy_impl import lombscargle_scipy as lombscargle_scipy
from .slow_impl import lombscargle_slow as lombscargle_slow
from . import (
    chi2_impl as chi2_impl,
    fast_impl as fast_impl,
    fastchi2_impl as fastchi2_impl,
    main as main,
    mle as mle,
    scipy_impl as scipy_impl,
    slow_impl as slow_impl,
    utils as utils,
    cython_impl as cython_impl,
    tests as tests,
)
