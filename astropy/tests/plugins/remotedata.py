# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
These plugins modify the behavior of py.test and are meant to be imported
into conftest.py in the root directory.
"""

import pytest
from .disable_internet import turn_off_internet, turn_on_internet


def pytest_addoption(parser):

    # The following means that if --remote-data is not specified, the default
    # is 'none', but if it is specified without arguments (--remote-data), it
    # defaults to '--remote-data=any'.
    parser.addoption("--remote-data", nargs="?", const='any', default='none',
                     help="run tests with online data")


def pytest_configure(config):
    config.getini('markers').append(
        'remote_data: Run tests that require data from remote servers')

    # Monkeypatch to deny access to remote resources unless explicitly told
    # otherwise

    if config.getoption('remote_data') != 'any':
        turn_off_internet(verbose=config.option.verbose,
                          allow_astropy_data=config.getoption('remote_data') == 'astropy')


def pytest_unconfigure():
    """
    Cleanup post-testing
    """
    # restore internet connectivity (only lost if remote_data=False and
    # turn_off_internet previously called)
    # this is harmless / does nothing if socket connections were never disabled
    turn_on_internet()


def pytest_runtest_setup(item):

    remote_data = item.keywords.get('remote_data')

    remote_data_config = item.config.getvalue("remote_data")

    if remote_data is not None:

        source = remote_data.kwargs.get('source', 'any')

        if source not in ('astropy', 'any'):
            raise ValueError("source should be 'astropy' or 'any'")

        if remote_data_config == 'none':
            pytest.skip("need --remote-data option to run")
        elif remote_data_config == 'astropy':
            if source == 'any':
                pytest.skip("need --remote-data option to run")
