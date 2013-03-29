# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
These plugins modify the behavior of py.test and are meant to be imported
into conftest.py in the root directory.
"""

import io
import os
import sys

from .helper import pytest

PY3K = sys.version_info[0] >= 3

# these pytest hooks allow us to mark tests and run the marked tests with
# specific command line options.
def pytest_addoption(parser):
    parser.addoption("--remote-data", action="store_true",
        help="run tests with online data")
    parser.addoption("--open-files", action="store_true",
        help="fail if any test leaves files open")


# Open file detection.
#
# This works by calling out to lsof to get the list of open files held
# by the process both before and after the test.  If something is
# still open after the test that wasn't open before the test, an
# AssertionError is raised.
#
# This is not thread-safe.  We're not currently running our tests
# multi-threaded, but that is worth noting.

SUPPORTS_OPEN_FILE_DETECTION = (
    sys.platform in ('linux', 'linux2', 'darwin'))


def _get_open_file_list():
    import imp
    import subprocess
    fsencoding = sys.getfilesystemencoding()

    sproc = subprocess.Popen(
        ['lsof -F0 -n -p {0}'.format(os.getpid())],
        shell=True, stdout=subprocess.PIPE)
    output = sproc.communicate()[0].strip()
    files = []
    for line in output.split(b'\n'):
        columns = line.split(b'\0')
        mapping = {}
        for column in columns:
            if len(column) >= 2:
                mapping[column[0:1]] = column[1:]

        if (mapping.get(b'f') and
            mapping.get(b'a', b' ') != b' ' and
            mapping.get(b't') == b'REG'):
            # Ignore extension modules -- they may be imported by a
            # test but are never again closed by the runtime.  That's
            # ok.
            for suffix, mode, filetype in imp.get_suffixes():
                if mapping[b'n'].decode(fsencoding).endswith(suffix):
                    break
            else:
                files.append(mapping[b'n'])

    return set(files)


def pytest_runtest_setup(item):
    # Store a list of the currently opened files so we can compare
    # against them when the test is done.
    if SUPPORTS_OPEN_FILE_DETECTION and item.config.getvalue('open_files'):
        item.open_files = _get_open_file_list()

    if ('remote_data' in item.keywords and
        not item.config.getvalue("remote_data")):
        pytest.skip("need --remote-data option to run")


if SUPPORTS_OPEN_FILE_DETECTION:
    def pytest_runtest_teardown(item, nextitem):
        # a "skipped" test will not have been called with
        # pytest_runtest_setup, so therefore won't have an
        # "open_files" member
        if (not item.config.getvalue('open_files') or
            not hasattr(item, 'open_files')):
            return

        start_open_files = item.open_files
        del item.open_files

        open_files = _get_open_file_list()

        # This works in tandem with the test_open_file_detection test to
        # ensure that it creates one extra open file.
        if item.name == 'test_open_file_detection':
            assert len(start_open_files) + 1 == len(open_files)
            return

        not_closed = set()
        for filename in open_files:
            # astropy.log files are allowed to continue to exist
            # between test runs
            if os.path.basename(filename) == 'astropy.log':
                continue

            if filename not in start_open_files:
                not_closed.add(filename)

        if len(not_closed):
            msg = [u'File(s) not closed:']
            for name in not_closed:
                msg.append(u'  {0}'.format(
                    name.decode(sys.getfilesystemencoding())))
            raise AssertionError(u'\n'.join(msg))


def pytest_report_header(config):

    from .. import __version__

    s = "\nRunning tests with Astropy version {0}.\n".format(__version__)
    s += "Running tests in {0}.\n\n".format(" ".join(config.args))

    from platform import platform
    s += "Platform: {0}\n\n".format(platform())
    s += "Executable: {0}\n\n".format(sys.executable)
    s += "Full Python Version: \n{0}\n\n".format(sys.version)

    import numpy
    s += "Numpy: {0}\n".format(numpy.__version__)

    try:
        import scipy
        s += "Scipy: {0}\n".format(scipy.__version__)
    except:
        s += "Scipy: not available\n"

    try:
        import matplotlib
        s += "Matplotlib: {0}\n".format(matplotlib.__version__)
    except:
        s += "Matplotlib: not available\n"

    try:
        import h5py
        s += "h5py: {0}\n".format(h5py.__version__)
    except:
        s += "h5py: not available\n"

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
