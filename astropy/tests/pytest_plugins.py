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


@pytest.fixture(autouse=True)
def modarg(request):
    """Sets up environment variables to fake the config and cache
    directories, then removes the temporary directories.

    Does nothing if we are inside the sphinx testing command, as it
    should have already done this for us.
    """
    import os
    import shutil
    import tempfile

    # check if we're inside the distutils test command, which sets the
    # _ASTROPY_TEST_ builtin
    try:
        _ASTROPY_TEST_
        insidetestcmd = True
    except NameError:
        insidetestcmd = False

    if not insidetestcmd:
        oldconfigdir = os.environ.get('XDG_CONFIG_HOME')
        oldcachedir = os.environ.get('XDG_CACHE_HOME')
        os.environ['XDG_CONFIG_HOME'] = tempfile.mkdtemp('astropy_config')
        os.environ['XDG_CACHE_HOME'] = tempfile.mkdtemp('astropy_cache')
        os.mkdir(os.path.join(os.environ['XDG_CONFIG_HOME'], 'astropy'))
        os.mkdir(os.path.join(os.environ['XDG_CACHE_HOME'], 'astropy'))

        def teardown():
            #wipe the config/cache tmpdirs and restore the envars
            shutil.rmtree(os.environ['XDG_CONFIG_HOME'])
            shutil.rmtree(os.environ['XDG_CACHE_HOME'])
            if oldconfigdir is None:
                del os.environ['XDG_CONFIG_HOME']
            else:
                os.environ['XDG_CONFIG_HOME'] = oldconfigdir
            if oldcachedir is None:
                del os.environ['XDG_CACHE_HOME']
            else:
                os.environ['XDG_CACHE_HOME'] = oldcachedir

        request.addfinalizer(teardown)
