# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This plugin provides command-line options for controlling whether and how tests
make use of online data.
"""
import pytest
from .disable_internet import turn_off_internet, turn_on_internet
from distutils.version import StrictVersion


def pytest_addoption(parser):

    # The following means that if --remote-data is not specified, the default
    # is 'none', but if it is specified without arguments (--remote-data), it
    # defaults to '--remote-data=any'.
    parser.addoption(
        "--remote-data", nargs="?", const='any', default='none',
        help="run tests with online data")

    parser.addini(
        'remote_data_strict',
        "If 'True', tests will fail if they attempt to access the internet "
        "but are not explicitly marked with 'remote_data'",
        type="bool", default=False)


def pytest_configure(config):
    config.getini('markers').append(
        'remote_data: Apply to tests that require data from remote servers')
    config.getini('markers').append(
        'internet_off: Apply to tests that should only run when '
        'network access is deactivated')

    strict_check = bool(config.getini('remote_data_strict'))

    remote_data = config.getoption('remote_data')
    if remote_data not in ['astropy', 'any', 'github', 'none']:
        raise pytest.UsageError(
            "'{}' is not a valid source for remote data".format(remote_data))

    # Monkeypatch to deny access to remote resources unless explicitly told
    # otherwise
    if strict_check:
        if remote_data == 'github':  # GitHub only
            turn_off_internet(
                verbose=config.option.verbose, allow_github_data=True)
        elif remote_data != 'any':  # Astropy + GitHub or nothing
            turn_off_internet(
                verbose=config.option.verbose,
                allow_astropy_data=(remote_data == 'astropy'))


def pytest_unconfigure():
    """
    Cleanup post-testing
    """
    # restore internet connectivity (only lost if remote_data=False and
    # turn_off_internet previously called)
    # this is harmless / does nothing if socket connections were never disabled
    turn_on_internet()


def pytest_runtest_setup(item):

    if StrictVersion(pytest.__version__) < StrictVersion("3.6"):
        remote_data = item.get_marker('remote_data')
        internet_off = item.get_marker('internet_off')
    else:
        remote_data = item.get_closest_marker('remote_data')
        internet_off = item.get_closest_marker('internet_off')

    remote_data_config = item.config.getvalue("remote_data")

    if remote_data is not None and internet_off is not None:
        raise ValueError("remote_data and internet_off are not compatible")

    if remote_data is not None:
        if len(remote_data.args) > 0:
            source = remote_data.args[0]
        else:
            source = remote_data.kwargs.get('source', 'any')

        if source not in ('astropy', 'any', 'github'):
            raise ValueError("source should be 'astropy', 'any', or 'github'")

        if remote_data_config == 'none':
            pytest.skip("need --remote-data option to run")
        elif remote_data_config == 'github':
            if source in ('any', 'astropy'):
                pytest.skip("need --remote-data or "
                            "--remote-data=astropy option to run")
        elif remote_data_config == 'astropy':
            # When --remote-data=astropy is given, skip tests simply
            # marked as --remote-data
            if source == 'any':
                pytest.skip("need --remote-data option to run")

    if internet_off is not None:
        if remote_data_config != 'none':
            pytest.skip("run this test only when network access is disabled")
