# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This file contains pytest configuration settings that are astropy-specific
(i.e.  those that would not necessarily be shared by affiliated packages
making use of astropy's test runner).
"""
import os
import builtins
import tempfile

from astropy.tests.plugins.display import PYTEST_HEADER_MODULES
from astropy.tests.helper import enable_deprecations_as_exceptions

try:
    import matplotlib
except ImportError:
    HAS_MATPLOTLIB = False
else:
    HAS_MATPLOTLIB = True

enable_deprecations_as_exceptions(
    include_astropy_deprecations=False,
    # This is a workaround for the OpenSSL deprecation warning that comes from
    # the `requests` module. It only appears when both asdf and sphinx are
    # installed. This can be removed once pyopenssl 1.7.20+ is released.
    modules_to_ignore_on_import=['requests'])

if HAS_MATPLOTLIB:
    matplotlib.use('Agg')

matplotlibrc_cache = {}


def pytest_configure(config):
    builtins._pytest_running = True
    # do not assign to matplotlibrc_cache in function scope
    if HAS_MATPLOTLIB:
        matplotlibrc_cache.update(matplotlib.rcParams)
        matplotlib.rcdefaults()

    # Make sure we use temporary directories for the config and cache
    # so that the tests are insensitive to local configuration. Note that this
    # is also set in the test runner, but we need to also set it here for
    # things to work properly in parallel mode

    builtins._xdg_config_home_orig = os.environ.get('XDG_CONFIG_HOME')
    builtins._xdg_cache_home_orig = os.environ.get('XDG_CACHE_HOME')

    os.environ['XDG_CONFIG_HOME'] = tempfile.mkdtemp('astropy_config')
    os.environ['XDG_CACHE_HOME'] = tempfile.mkdtemp('astropy_cache')

    os.mkdir(os.path.join(os.environ['XDG_CONFIG_HOME'], 'astropy'))
    os.mkdir(os.path.join(os.environ['XDG_CACHE_HOME'], 'astropy'))


def pytest_unconfigure(config):
    builtins._pytest_running = False
    # do not assign to matplotlibrc_cache in function scope
    if HAS_MATPLOTLIB:
        matplotlib.rcParams.update(matplotlibrc_cache)
        matplotlibrc_cache.clear()

    if builtins._xdg_config_home_orig is None:
        os.environ.pop('XDG_CONFIG_HOME')
    else:
        os.environ['XDG_CONFIG_HOME'] = builtins._xdg_config_home_orig

    if builtins._xdg_cache_home_orig is None:
        os.environ.pop('XDG_CACHE_HOME')
    else:
        os.environ['XDG_CACHE_HOME'] = builtins._xdg_cache_home_orig


PYTEST_HEADER_MODULES['Cython'] = 'cython'
PYTEST_HEADER_MODULES['Scikit-image'] = 'skimage'
