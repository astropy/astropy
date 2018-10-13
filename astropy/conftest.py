# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This file contains pytest configuration settings that are astropy-specific
(i.e.  those that would not necessarily be shared by affiliated packages
making use of astropy's test runner).
"""
from .extern.six.moves import builtins

from .tests.pytest_plugins import *

try:
    import matplotlib
except ImportError:
    HAS_MATPLOTLIB = False
else:
    HAS_MATPLOTLIB = True

enable_deprecations_as_exceptions(include_astropy_deprecations=False)


if HAS_MATPLOTLIB:
    matplotlib.use('Agg')


matplotlibrc_cache = {}


def pytest_configure(config):
    builtins._pytest_running = True
    # do not assign to matplotlibrc_cache in function scope
    if HAS_MATPLOTLIB:
        matplotlibrc_cache.update(matplotlib.rcParams)
        matplotlib.rcdefaults()


def pytest_unconfigure(config):
    builtins._pytest_running = False
    # do not assign to matplotlibrc_cache in function scope
    if HAS_MATPLOTLIB:
        matplotlib.rcParams.update(matplotlibrc_cache)
        matplotlibrc_cache.clear()


PYTEST_HEADER_MODULES['Cython'] = 'cython'
