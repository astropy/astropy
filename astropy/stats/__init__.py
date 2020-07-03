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

from . import funcs
from .funcs import *  # noqa
from . import biweight
from .biweight import *  # noqa
from . import sigma_clipping
from .sigma_clipping import *  # noqa
from . import jackknife
from .jackknife import *  # noqa
from . import circstats
from .circstats import *  # noqa
from . import bayesian_blocks as _bb
from .bayesian_blocks import *  # noqa
from . import histogram as _hist
from .histogram import *  # noqa
from . import info_theory
from .info_theory import *  # noqa
from . import spatial
from .spatial import *  # noqa
from .lombscargle import *  # noqa
from .bls import *  # noqa

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
