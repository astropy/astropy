# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This subpackage contains statistical tools provided for or used by Astropy.

While the `scipy.stats` package contains a wide range of statistical
tools, it is a general-purpose package, and is missing some that are
particularly useful to astronomy or are used in an atypical way in
astronomy. This package is intended to provide such functionality, but
*not* to replace `scipy.stats` if its implementation satisfies
astronomers' needs.

"""

from . import bayesian_blocks as _bb
from . import (
    biweight,
    circstats,
    funcs,
    info_theory,
    jackknife,
    sigma_clipping,
    spatial,
)
from . import histogram as _hist
from .bayesian_blocks import *
from .biweight import *
from .circstats import *
from .funcs import *
from .histogram import *
from .info_theory import *
from .jackknife import *
from .sigma_clipping import *
from .spatial import *

__all__ = []
__all__ += funcs.__all__
__all__ += biweight.__all__
__all__ += sigma_clipping.__all__
__all__ += jackknife.__all__
__all__ += circstats.__all__
__all__ += _bb.__all__
__all__ += _hist.__all__
__all__ += info_theory.__all__
__all__ += spatial.__all__
