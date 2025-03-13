# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This file contains pytest configuration settings that are astropy-specific
(i.e.  those that would not necessarily be shared by affiliated packages
making use of astropy's test runner).
"""

import builtins
import os
import sys
import tempfile
import warnings
from pathlib import Path

try:
    from pytest_astropy_header.display import PYTEST_HEADER_MODULES, TESTED_VERSIONS
except ImportError:
    PYTEST_HEADER_MODULES = {}
    TESTED_VERSIONS = {}

import pytest

from astropy import __version__
from astropy.utils.compat.optional_deps import HAS_MATPLOTLIB

if HAS_MATPLOTLIB:
    import matplotlib as mpl

matplotlibrc_cache = {}


@pytest.fixture
def ignore_matplotlibrc():
    # This is a fixture for tests that use matplotlib but not pytest-mpl
    # (which already handles rcParams)
    from matplotlib import pyplot as plt

    with plt.style.context({}, after_reset=True):
        yield


@pytest.fixture
def fast_thread_switching():
    """Fixture that reduces thread switching interval.

    This makes it easier to provoke race conditions.
    """
    old = sys.getswitchinterval()
    sys.setswitchinterval(1e-6)
    yield
    sys.setswitchinterval(old)


def pytest_configure(config):
    # Ensure number of columns and lines is deterministic for testing
    from astropy import conf

    conf.max_width = 80
    conf.max_lines = 24

    # Disable IERS auto download for testing
    from astropy.utils.iers import conf as iers_conf

    iers_conf.auto_download = False

    builtins._pytest_running = True
    # do not assign to matplotlibrc_cache in function scope
    if HAS_MATPLOTLIB:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            matplotlibrc_cache.update(mpl.rcParams)
            mpl.rcdefaults()
            mpl.use("Agg")

    # Make sure we use temporary directories for the config and cache
    # so that the tests are insensitive to local configuration. Note that this
    # is also set in the test runner, but we need to also set it here for
    # things to work properly in parallel mode

    builtins._xdg_config_home_orig = os.environ.get("XDG_CONFIG_HOME")
    builtins._xdg_cache_home_orig = os.environ.get("XDG_CACHE_HOME")

    os.environ["XDG_CONFIG_HOME"] = tempfile.mkdtemp("astropy_config")
    os.environ["XDG_CACHE_HOME"] = tempfile.mkdtemp("astropy_cache")

    Path(os.environ["XDG_CONFIG_HOME"]).joinpath("astropy").mkdir()
    Path(os.environ["XDG_CACHE_HOME"]).joinpath("astropy").mkdir()

    config.option.astropy_header = True
    PYTEST_HEADER_MODULES["PyERFA"] = "erfa"
    PYTEST_HEADER_MODULES["Cython"] = "cython"
    PYTEST_HEADER_MODULES["Scikit-image"] = "skimage"
    PYTEST_HEADER_MODULES["asdf-astropy"] = "asdf_astropy"
    TESTED_VERSIONS["Astropy"] = __version__

    # Limit the number of threads used by each worker when pytest-xdist is in
    # use.  Lifted from https://github.com/scipy/scipy/pull/14441
    # and https://github.com/scikit-learn/scikit-learn/pull/25918
    try:
        from threadpoolctl import threadpool_limits
    except ImportError:
        pass
    else:
        xdist_worker_count = os.environ.get("PYTEST_XDIST_WORKER_COUNT")
        if xdist_worker_count is not None:
            # use number of physical cores, assume hyperthreading
            max_threads = os.cpu_count() // 2
            threads_per_worker = max(max_threads // int(xdist_worker_count), 1)
            threadpool_limits(threads_per_worker)


def pytest_unconfigure(config):
    # Undo settings related to number of lines/columns to show
    from astropy import conf

    conf.reset("max_width")
    conf.reset("max_lines")

    # Undo IERS auto download setting for testing
    from astropy.utils.iers import conf as iers_conf

    iers_conf.reset("auto_download")

    builtins._pytest_running = False
    # do not assign to matplotlibrc_cache in function scope
    if HAS_MATPLOTLIB:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            mpl.rcParams.update(matplotlibrc_cache)
            matplotlibrc_cache.clear()

    if builtins._xdg_config_home_orig is None:
        os.environ.pop("XDG_CONFIG_HOME")
    else:
        os.environ["XDG_CONFIG_HOME"] = builtins._xdg_config_home_orig

    if builtins._xdg_cache_home_orig is None:
        os.environ.pop("XDG_CACHE_HOME")
    else:
        os.environ["XDG_CACHE_HOME"] = builtins._xdg_cache_home_orig


def pytest_terminal_summary(terminalreporter):
    """Output a warning to IPython users in case any tests failed."""
    try:
        # this name is only defined if running within ipython/jupyter
        __IPYTHON__  # noqa: B018
    except NameError:
        return

    if not terminalreporter.stats.get("failed"):
        # Only issue the warning when there are actually failures
        return

    terminalreporter.ensure_newline()
    terminalreporter.write_line(
        "Some tests may fail when run from the IPython prompt; "
        "especially, but not limited to tests involving logging and warning "
        "handling.  Unless you are certain as to the cause of the failure, "
        "please check that the failure occurs outside IPython as well.  See "
        "https://docs.astropy.org/en/stable/known_issues.html#failing-logging-"
        "tests-when-running-the-tests-in-ipython for more information.",
        yellow=True,
        bold=True,
    )
