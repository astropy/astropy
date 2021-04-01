# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""This file contains pytest configuration settings that are astropy-specific
(i.e.  those that would not necessarily be shared by affiliated packages
making use of astropy's test runner).

This is the file usually picked up by pytest and tox.
"""
import os
import builtins
import sys
import tempfile
import warnings

import hypothesis

from astropy import __version__
from astropy.tests.helper import enable_deprecations_as_exceptions

try:
    from pytest_astropy_header.display import PYTEST_HEADER_MODULES, TESTED_VERSIONS
except ImportError:
    PYTEST_HEADER_MODULES = {}
    TESTED_VERSIONS = {}

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


# This has to be in the root dir or it will not display in CI.
def pytest_report_header(config):
    # This gets added after the pytest-astropy-header output.
    return (f'ARCH_ON_CI: {os.environ.get("ARCH_ON_CI", "undefined")}\n'
            f'IS_CRON: {os.environ.get("IS_CRON", "undefined")}\n')


# Tell Hypothesis that we might be running slow tests, to print the seed blob
# so we can easily reproduce failures from CI, and derive a fuzzing profile
# to try many more inputs when we detect a scheduled build or when specifically
# requested using the HYPOTHESIS_PROFILE=fuzz environment variable or
# `pytest --hypothesis-profile=fuzz ...` argument.

hypothesis.settings.register_profile(
    'ci', deadline=None, print_blob=True, derandomize=True
)
hypothesis.settings.register_profile(
    'fuzzing', deadline=None, print_blob=True, max_examples=1000
)
default = 'fuzzing' if (os.environ.get('IS_CRON') == 'true' and os.environ.get('ARCH_ON_CI') not in ('aarch64', 'ppc64le')) else 'ci'  # noqa: E501
hypothesis.settings.load_profile(os.environ.get('HYPOTHESIS_PROFILE', default))

# Make sure we use temporary directories for the config and cache
# so that the tests are insensitive to local configuration.

os.environ['XDG_CONFIG_HOME'] = tempfile.mkdtemp('astropy_config')
os.environ['XDG_CACHE_HOME'] = tempfile.mkdtemp('astropy_cache')

os.mkdir(os.path.join(os.environ['XDG_CONFIG_HOME'], 'astropy'))
os.mkdir(os.path.join(os.environ['XDG_CACHE_HOME'], 'astropy'))

# Note that we don't need to change the environment variables back or remove
# them after testing, because they are only changed for the duration of the
# Python process, and this configuration only matters if running pytest
# directly, not from e.g. an IPython session.
