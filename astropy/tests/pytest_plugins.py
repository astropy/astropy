# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
These plugins modify the behavior of py.test and are meant to be imported
into conftest.py in the root directory.
"""

from .helper import pytest

# these pytest hooks allow us to mark tests and run the marked tests with
# specific command line options.
def pytest_addoption(parser):
    parser.addoption("--remote-data", action="store_true",
        help="run tests with online data")


def pytest_runtest_setup(item):
    if ('remote_data' in item.keywords and
        not item.config.getvalue("remote_data")):
        pytest.skip("need --remote-data option to run")


def pytest_report_header(config):
    from .. import __version__
    s = "\nTesting Astropy version {0}.\n".format(__version__)
    s += "Running tests in {0}.\n".format(" ".join(config.args))

    special_opts = ["remote_data", "pep8"]
    opts = []
    for op in special_opts:
        if getattr(config.option, op, None):
            opts.append(op)
    if opts:
        s += "Using Astropy options: {0}.\n".format(" ".join(opts))

    return s
