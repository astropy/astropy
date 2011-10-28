# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module contains a number of utilities for use during
setup/build/packaging that are useful to astropy as a whole.
"""

import os
import sys
from warnings import warn

from distutils.dist import Distribution
from distutils.errors import DistutilsError

try:
    import Cython
    HAVE_CYTHON = True
except ImportError:
    HAVE_CYTHON = False


def get_debug_option():
    # Pre-parse the Distutils command-line options and config files to
    # determine whether or not the --debug option was set.
    # Normally the only options that allow the debug option are build,
    # build_ext, and build_clib.  So for the sake of generality it's best to
    # just use the build command.
    dist = Distribution()
    try:
        dist.parse_config_files()
        dist.parse_command_line()
    except DistutilsError:
        # Let distutils handle this itself
        return None

    debug_build_cmds = ['build', 'build_ext', 'build_clib']

    # First ensure that one of the commands that accepts --debug is being run
    for cmd in debug_build_cmds:
        if cmd in dist.commands:
            break
    else:
        return None

    for cmd in debug_build_cmds:
        cmd_opts = dist.get_option_dict(cmd)
        if 'debug' in cmd_opts:
            debug = bool(cmd_opts['debug'][1])
            break
    else:
        debug = False

    try:
        from astropy.version import debug as current_debug
    except ImportError:
        current_debug = None

    if current_debug is not None and current_debug != debug:
        # Force rebuild of extension modules
        sys.argv.extend(['build', '--force'])

    return debug


def iter_setup_packages():
    # Find all of the setup_package.py modules, import them, and add them
    # to the setup_packages list.
    for root, dirs, files in os.walk('astropy'):
        if 'setup_package.py' in files:
            name = root.replace(os.path.sep, '.') + '.setup_package'
            module = import_module(name)
            yield module


def iter_pyx_files():
    for dirpath, dirnames, filenames in os.walk('astropy'):
        modbase = dirpath.replace(os.sep, '.')
        for fn in filenames:
            if fn.endswith('.pyx'):
                fullfn = os.path.join(dirpath, fn)
                # Package must match file nam
                extmod = modbase + '.' + fn[:-4]
                yield (extmod, fullfn)


def get_cython_extensions():
    # Look for Cython files - compile with Cython if it is not a release
    # and Cython is installed. Otherwise, use the .c files that live next
    # to the Cython files.
    from astropy.version import release

    ext_modules = []
    if not release and HAVE_CYTHON:
        # Add .pyx files
        for extmod, pyxfn in iter_pyx_files():
            ext_modules.append(Extension(extmod, [pyxfn]))
    else:
        # Add .c files
        for extmod, pyxfn in iter_pyx_files():
            cfn = pyxfn[:-4] + '.c'
            if os.path.exists(cfn):
                ext_modules.append(Extension(extmod, [cfn]))
            else:
                warn('Could not find Cython-generated C extension {0} - '
                     'The {1} module will be skipped.'.format(cfn, extmod))
    return ext_modules


def write_if_different(filename, data):
    """
    Write *data* to *filename*, if the content of the file is
    different.
    """
    assert isinstance(data, bytes)

    if os.path.exists(filename):
        with open(filename, 'rb') as fd:
            original_data = fd.read()
    else:
        original_data = None

    if original_data != data:
        with open(filename, 'wb') as fd:
            fd.write(data)


def check_numpy():
    """
    Check that Numpy is installed and it is of the minimum version we
    require.
    """
    import numpy

    major, minor, rest = numpy.__version__.split(".", 2)
    if (int(major), int(minor)) < (1, 3):
        raise ImportError("numpy version 1.3 or later must be installed to build astropy")


def get_numpy_include_path():
    """
    Gets the path to the numpy headers
    """
    import numpy

    try:
        numpy_include = numpy.get_include()
    except AttributeError:
        numpy_include = numpy.get_numpy_include()
    return numpy_include



################################################################################
# Backport of importlib.import_module from 3.x.  This backport was provided by
# Brett Cannon and downloaded from here:
#    http://pypi.python.org/pypi/importlib/1.0.1
# While not critical (and in no way guaranteed!), it would be nice to keep this
# code compatible with Python 2.3.

try:
    from importlib import import_module
except ImportError:

    def _resolve_name(name, package, level):
        """Return the absolute name of the module to be imported."""
        if not hasattr(package, 'rindex'):
            raise ValueError("'package' not set to a string")
        dot = len(package)
        for x in xrange(level, 1, -1):
            try:
                dot = package.rindex('.', 0, dot)
            except ValueError:
                raise ValueError("attempted relative import beyond top-level "
                                  "package")
        return "%s.%s" % (package[:dot], name)


    def import_module(name, package=None):
        """Import a module.

        The 'package' argument is required when performing a relative
        import. It specifies the package to use as the anchor point
        from which to resolve the relative import to an absolute
        import.
        """
        if name.startswith('.'):
            if not package:
                raise TypeError(
                    "relative imports require the 'package' argument")
            level = 0
            for character in name:
                if character != '.':
                    break
                level += 1
            name = _resolve_name(name[level:], package, level)
        __import__(name)
        return sys.modules[name]
