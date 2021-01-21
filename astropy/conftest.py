# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This file contains pytest configuration settings that are astropy-specific
(i.e.  those that would not necessarily be shared by affiliated packages
making use of astropy's test runner).
"""
import os
import builtins
import sys
import tempfile
import warnings

try:
    from pytest_astropy_header.display import PYTEST_HEADER_MODULES, TESTED_VERSIONS
except ImportError:
    PYTEST_HEADER_MODULES = {}
    TESTED_VERSIONS = {}

import pytest

from astropy import __version__
from astropy.tests.helper import enable_deprecations_as_exceptions

try:
    # This is needed to silence a warning from matplotlib caused by
    # PyInstaller's matplotlib runtime hook.  This can be removed once the
    # issue is fixed upstream in PyInstaller, and only impacts us when running
    # the tests from a PyInstaller bundle.
    # See https://github.com/astropy/astropy/issues/10785
    if getattr(sys, 'frozen', False) and hasattr(sys, '_MEIPASS'):
        # The above checks whether we are running in a PyInstaller bundle.
        warnings.filterwarnings("ignore", "(?s).*MATPLOTLIBDATA.*",
                                category=UserWarning)
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

matplotlibrc_cache = {}


@pytest.fixture
def ignore_matplotlibrc():
    # This is a fixture for tests that use matplotlib but not pytest-mpl
    # (which already handles rcParams)
    from matplotlib import pyplot as plt
    with plt.style.context({}, after_reset=True):
        yield


def pytest_configure(config):
    from astropy.utils.iers import conf as iers_conf

    # Disable IERS auto download for testing
    iers_conf.auto_download = False

    builtins._pytest_running = True
    # do not assign to matplotlibrc_cache in function scope
    if HAS_MATPLOTLIB:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            matplotlibrc_cache.update(matplotlib.rcParams)
            matplotlib.rcdefaults()
            matplotlib.use('Agg')

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

    config.option.astropy_header = True
    PYTEST_HEADER_MODULES['PyERFA'] = 'erfa'
    PYTEST_HEADER_MODULES['Cython'] = 'cython'
    PYTEST_HEADER_MODULES['Scikit-image'] = 'skimage'
    PYTEST_HEADER_MODULES['asdf'] = 'asdf'
    TESTED_VERSIONS['Astropy'] = __version__


def pytest_unconfigure(config):
    from astropy.utils.iers import conf as iers_conf

    # Undo IERS auto download setting for testing
    iers_conf.reset('auto_download')

    builtins._pytest_running = False
    # do not assign to matplotlibrc_cache in function scope
    if HAS_MATPLOTLIB:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
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


def pytest_terminal_summary(terminalreporter):
    """Output a warning to IPython users in case any tests failed."""

    try:
        get_ipython()
    except NameError:
        return

    if not terminalreporter.stats.get('failed'):
        # Only issue the warning when there are actually failures
        return

    terminalreporter.ensure_newline()
    terminalreporter.write_line(
        'Some tests are known to fail when run from the IPython prompt; '
        'especially, but not limited to tests involving logging and warning '
        'handling.  Unless you are certain as to the cause of the failure, '
        'please check that the failure occurs outside IPython as well.  See '
        'https://docs.astropy.org/en/stable/known_issues.html#failing-logging-'
        'tests-when-running-the-tests-in-ipython for more information.',
        yellow=True, bold=True)
