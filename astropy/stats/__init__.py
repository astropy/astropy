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
from . import biweight, circstats, funcs
from . import histogram as _hist
from . import info_theory, jackknife, sigma_clipping, spatial
from .bayesian_blocks import *
from .biweight import *
from .bls import *
from .circstats import *
from .funcs import *
from .histogram import *
from .info_theory import *
from .jackknife import *
from .lombscargle import *
from .sigma_clipping import *
from .spatial import *

# This is to avoid importing deprecated modules in subpackage star import
__all__ = []
__all__.extend(funcs.__all__)
__all__.extend(biweight.__all__)
__all__.extend(sigma_clipping.__all__)
__all__.extend(jackknife.__all__)
__all__.extend(circstats.__all__)
__all__.extend(_bb.__all__)
__all__.extend(_hist.__all__)
__all__.extend(info_theory.__all__)
__all__.extend(spatial.__all__)
