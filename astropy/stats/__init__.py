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

from .funcs import *  # noqa
from .biweight import *  # noqa
from .sigma_clipping import *  # noqa
from .jackknife import *  # noqa
from .circstats import *  # noqa
from .bayesian_blocks import *  # noqa
from .histogram import *  # noqa
from .info_theory import *  # noqa
from .lombscargle import *  # noqa
from .spatial import *  # noqa
from .bls import *  # noqa
