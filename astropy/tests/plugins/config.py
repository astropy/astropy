# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This plugin provides customization of configuration and cache directories used
by pytest.
"""
import datetime
import locale
import os
import sys
from collections import OrderedDict

import pytest

from ...config.paths import set_temp_config, set_temp_cache
from ...utils.argparse import writeable_directory
from ..helper import treat_deprecations_as_exceptions

import importlib.machinery as importlib_machinery


# these pytest hooks allow us to mark tests and run the marked tests with
# specific command line options.

def pytest_addoption(parser):

    parser.addoption("--astropy-config-dir", nargs='?', type=writeable_directory,
                     help="specify directory for storing and retrieving the "
                          "Astropy configuration during tests (default is "
                          "to use a temporary directory created by the test "
                          "runner); be aware that using an Astropy config "
                          "file other than the default can cause some tests "
                          "to fail unexpectedly")

    parser.addoption("--astropy-cache-dir", nargs='?', type=writeable_directory,
                     help="specify directory for storing and retrieving the "
                          "Astropy cache during tests (default is "
                          "to use a temporary directory created by the test "
                          "runner)")
    parser.addini("astropy_config_dir",
                  "specify directory for storing and retrieving the "
                  "Astropy configuration during tests (default is "
                  "to use a temporary directory created by the test "
                  "runner); be aware that using an Astropy config "
                  "file other than the default can cause some tests "
                  "to fail unexpectedly", default=None)

    parser.addini("astropy_cache_dir",
                  "specify directory for storing and retrieving the "
                  "Astropy cache during tests (default is "
                  "to use a temporary directory created by the test "
                  "runner)", default=None)

def pytest_configure(config):
    treat_deprecations_as_exceptions()

def pytest_runtest_setup(item):
    config_dir = item.config.getini('astropy_config_dir')
    cache_dir = item.config.getini('astropy_cache_dir')

    # Command-line options can override, however
    config_dir = item.config.getoption('astropy_config_dir') or config_dir
    cache_dir = item.config.getoption('astropy_cache_dir') or cache_dir

    # We can't really use context managers directly in py.test (although
    # py.test 2.7 adds the capability), so this may look a bit hacky
    if config_dir:
        item.set_temp_config = set_temp_config(config_dir)
        item.set_temp_config.__enter__()
    if cache_dir:
        item.set_temp_cache = set_temp_cache(cache_dir)
        item.set_temp_cache.__enter__()

def pytest_runtest_teardown(item, nextitem):
    if hasattr(item, 'set_temp_cache'):
        item.set_temp_cache.__exit__()
    if hasattr(item, 'set_temp_config'):
        item.set_temp_config.__exit__()
