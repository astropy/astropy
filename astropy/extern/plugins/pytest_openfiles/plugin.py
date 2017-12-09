# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This plugin provides support for testing whether file-like objects are properly
closed.
"""
import imp
import os

try:
    import importlib.machinery as importlib_machinery
except ImportError:
    importlib_machinery = None


def pytest_addoption(parser):

    parser.addoption("--open-files", action="store_true",
                     help="fail if any test leaves files open")

    parser.addini("open_files_ignore",
                  "when used with the --open-files option, allows "
                  "specifying names of files that may be ignored when "
                  "left open between tests--files in this list are matched "
                  "may be specified by their base name (ignoring their full "
                  "path) or by absolute path", type="args", default=())

# Open file detection.
#
# This works by calling out to psutil to get the list of open files
# held by the process both before and after the test.  If something is
# still open after the test that wasn't open before the test, an
# AssertionError is raised.
#
# This is not thread-safe.  We're not currently running our tests
# multi-threaded, but that is worth noting.


def _get_open_file_list():
    import psutil
    files = []
    p = psutil.Process()

    if importlib_machinery is not None:
        suffixes = tuple(importlib_machinery.all_suffixes())
    else:
        suffixes = tuple(info[0] for info in imp.get_suffixes())

    files = [x.path for x in p.open_files() if not x.path.endswith(suffixes)]

    return set(files)


def pytest_runtest_setup(item):
    # Store a list of the currently opened files so we can compare
    # against them when the test is done.
    if item.config.getvalue('open_files'):
        item.open_files = _get_open_file_list()


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
    open_files_ignore = item.config.getini('open_files_ignore')
    for filename in open_files:
        ignore = False

        for ignored in open_files_ignore:
            if not os.path.isabs(ignored):
                if os.path.basename(filename) == ignored:
                    ignore = True
                    break
            else:
                if filename == ignored:
                    ignore = True
                    break

        if ignore:
            continue

        if filename not in start_open_files:
            not_closed.add(filename)

    if len(not_closed):
        msg = ['File(s) not closed:']
        for name in not_closed:
            msg.append('  {0}'.format(name))
        raise AssertionError('\n'.join(msg))
