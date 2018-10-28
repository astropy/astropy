# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This plugin provides customization of the header displayed by pytest for
reporting purposes.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import os
import sys
import datetime
import locale
import math
from collections import OrderedDict

from ..helper import ignore_warnings
from ...utils.introspection import resolve_name


PYTEST_HEADER_MODULES = OrderedDict([('Numpy', 'numpy'),
                                     ('Scipy', 'scipy'),
                                     ('Matplotlib', 'matplotlib'),
                                     ('h5py', 'h5py'),
                                     ('Pandas', 'pandas')])

# This always returns with Astropy's version
from ... import __version__
TESTED_VERSIONS = OrderedDict([('Astropy', __version__)])


def pytest_report_header(config):

    try:
        stdoutencoding = sys.stdout.encoding or 'ascii'
    except AttributeError:
        stdoutencoding = 'ascii'

    args = config.args

    # TESTED_VERSIONS can contain the affiliated package version, too
    if len(TESTED_VERSIONS) > 1:
        for pkg, version in TESTED_VERSIONS.items():
            if pkg != 'Astropy':
                s = "\nRunning tests with {0} version {1}.\n".format(
                    pkg, version)
    else:
        s = "\nRunning tests with Astropy version {0}.\n".format(
            TESTED_VERSIONS['Astropy'])

    # Per https://github.com/astropy/astropy/pull/4204, strip the rootdir from
    # each directory argument
    if hasattr(config, 'rootdir'):
        rootdir = str(config.rootdir)
        if not rootdir.endswith(os.sep):
            rootdir += os.sep

        dirs = [arg[len(rootdir):] if arg.startswith(rootdir) else arg
                for arg in args]
    else:
        dirs = args

    s += "Running tests in {0}.\n\n".format(" ".join(dirs))

    s += "Date: {0}\n\n".format(datetime.datetime.now().isoformat()[:19])

    from platform import platform
    plat = platform()
    if isinstance(plat, bytes):
        plat = plat.decode(stdoutencoding, 'replace')
    s += "Platform: {0}\n\n".format(plat)
    s += "Executable: {0}\n\n".format(sys.executable)
    s += "Full Python Version: \n{0}\n\n".format(sys.version)

    s += "encodings: sys: {0}, locale: {1}, filesystem: {2}".format(
        sys.getdefaultencoding(),
        locale.getpreferredencoding(),
        sys.getfilesystemencoding())
    s += '\n'

    s += "byteorder: {0}\n".format(sys.byteorder)
    s += "float info: dig: {0.dig}, mant_dig: {0.dig}\n\n".format(
        sys.float_info)

    for module_display, module_name in PYTEST_HEADER_MODULES.items():
        try:
            with ignore_warnings(DeprecationWarning):
                module = resolve_name(module_name)
        except ImportError:
            s += "{0}: not available\n".format(module_display)
        else:
            try:
                version = module.__version__
            except AttributeError:
                version = 'unknown (no __version__ attribute)'
            s += "{0}: {1}\n".format(module_display, version)

    # Helpers version
    try:
        from ...version import astropy_helpers_version
    except ImportError:
        pass
    else:
        s += "astropy_helpers: {0}\n".format(astropy_helpers_version)

    special_opts = ["remote_data", "pep8"]
    opts = []
    for op in special_opts:
        op_value = getattr(config.option, op, None)
        if op_value:
            if isinstance(op_value, str):
                op = ': '.join((op, op_value))
            opts.append(op)
    if opts:
        s += "Using Astropy options: {0}.\n".format(", ".join(opts))

    return s


def pytest_terminal_summary(terminalreporter):
    """Output a warning to IPython users in case any tests failed."""

    try:
        get_ipython()
    except NameError:
        return

    if not terminalreporter.stats.get('failed'):
        # Only issue the warning when there are actually failures
        return

    terminalreporter.ensure_newline()
    terminalreporter.write_line(
        'Some tests are known to fail when run from the IPython prompt; '
        'especially, but not limited to tests involving logging and warning '
        'handling.  Unless you are certain as to the cause of the failure, '
        'please check that the failure occurs outside IPython as well.  See '
        'http://docs.astropy.org/en/stable/known_issues.html#failing-logging-'
        'tests-when-running-the-tests-in-ipython for more information.',
        yellow=True, bold=True)
