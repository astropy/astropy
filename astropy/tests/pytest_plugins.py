# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
These plugins modify the behavior of py.test and are meant to be imported
into conftest.py in the root directory.
"""

import io
import os
import sys
import traceback

from .helper import pytest

PY3K = sys.version_info[0] >= 3

# these pytest hooks allow us to mark tests and run the marked tests with
# specific command line options.
def pytest_addoption(parser):
    parser.addoption("--remote-data", action="store_true",
        help="run tests with online data")

# Open file detection.
#
# This works by replacing the built-in "open" with one that keeps
# track of all of the currently open files.  When each test is started,
# this list is copied.  At the end of the test, that list is compared
# with the set of currently open files.  If the latter is larger, the
# newly opened files that were not closed are displayed in an error
# message and the test fails.
#
# This is not thread-safe.  We're not currently running our tests
# multi-threaded, but that is worth noting.
#
# Note that this doesn't work with Python 3.  Since Python 3 has a
# hierarchy of classes for I/O written in C, it doesn't seem possible
# to monkey-patch the base file class with our own in order to keep
# track of opening and closing like we do here.  Further investigation
# and perhaps a wholly different approach is needed.  In the meantime,
# this functionality is just turned off for Python 3.
if not PY3K:
    import __builtin__
    _open_files = dict()
    _old_file = __builtin__.file

    class _new_file(_old_file):
        def __init__(self, *args, **kwargs):
            self.x = args[0]
            _old_file.__init__(self, *args, **kwargs)
            __builtin__.open = _old_open
            _open_files[self] = traceback.format_stack()
            __builtin__.open = _new_open

        def close(self):
            _old_file.close(self)
            try:
                del _open_files[self]
            except KeyError:
                pass

    _old_open = __builtin__.open
    def _new_open(*args, **kwargs):
        return _new_file(*args, **kwargs)
    __builtin__.open = _new_open
    __builtin__.file = _new_file


def pytest_runtest_setup(item):
    # Store a list of the currently opened files so we can compare
    # against them when the test is done.
    if not PY3K:
        item.open_files = dict(_open_files)

    if ('remote_data' in item.keywords and
        not item.config.getvalue("remote_data")):
        pytest.skip("need --remote-data option to run")


if not PY3K:
    def pytest_runtest_teardown(item, nextitem):
        # a "skipped" test will not have been called with
        # pytest_runtest_setup, so therefore won't have an
        # "open_files" member
        if not hasattr(item, 'open_files'):
            return

        start_open_files = item.open_files
        del item.open_files

        # This works in tandem with the test_open_file_detection test to
        # ensure that it creates one extra open file.
        if item.name == 'test_open_file_detection':
            assert len(start_open_files) + 1 == len(_open_files)
            return

        if len(start_open_files) != len(_open_files):
            not_closed = {}
            for fd, stack in _open_files.items():
                # astropy.log files are allowed to continue to exist
                # between test runs
                if os.path.basename(fd.name) == 'astropy.log':
                    continue

                if fd not in start_open_files:
                    not_closed[fd.name] = stack

            if len(not_closed):
                msg = ['File(s) not closed:']
                for name, stack in not_closed.items():
                    msg.append('  {0} opened from:'.format(name))
                    msg.extend('    {0}'.format(x) for x in stack)
                raise AssertionError('\n'.join(msg))


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
