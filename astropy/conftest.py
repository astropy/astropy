# this contains imports plugins that configure py.test for astropy tests.
# by importing them here in conftest.py they are discoverable by py.test
# no matter how it is invoked within the astropy tree.

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
    # do not assign to matplotlibrc_cache in function scope
    if HAS_MATPLOTLIB:
        matplotlibrc_cache.update(matplotlib.rcParams)
        matplotlib.rcdefaults()


def pytest_unconfigure(config):
    # do not assign to matplotlibrc_cache in function scope
    if HAS_MATPLOTLIB:
        matplotlib.rcParams.update(matplotlibrc_cache)
        matplotlibrc_cache.clear()


PYTEST_HEADER_MODULES['Cython'] = 'cython'
