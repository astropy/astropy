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
from distutils.core import Extension

try:
    import Cython
    HAVE_CYTHON = True
except ImportError:
    HAVE_CYTHON = False


def get_distutils_option(option, commands):
    """
    Returns the value of the given distutils option.  If the option is not
    set, returns None.

    Parameters
    ----------
    option : str
        The name of the option

    commands : list of str
        The list of commands on which this option is available
    """
    # Pre-parse the Distutils command-line options and config files to
    # if the option is set.
    dist = Distribution()
    try:
        dist.parse_config_files()
        dist.parse_command_line()
    except DistutilsError:
        # Let distutils handle this itself
        return None

    for cmd in commands:
        if cmd in dist.commands:
            break
    else:
        return None

    for cmd in commands:
        cmd_opts = dist.get_option_dict(cmd)
        if option in cmd_opts:
            return cmd_opts[option][1]
    else:
        return None


def get_debug_option():
    debug = bool(get_distutils_option(
        'debug', ['build', 'build_ext', 'build_clib']))

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


def iter_pyx_files(srcdir):
    for dirpath, dirnames, filenames in os.walk(srcdir):
        modbase = dirpath.replace(os.sep, '.')
        for fn in filenames:
            if fn.endswith('.pyx'):
                fullfn = os.path.join(dirpath, fn)
                # Package must match file name
                extmod = modbase + '.' + fn[:-4]
                yield (extmod, fullfn)


def get_cython_extensions(srcdir,prevextensions=tuple(),extincludedirs=None):
    """ Looks for Cython files and generates Extensions if needed.
    
    Parameters
    ----------
    srcdir : str
        Path to the root of the source directory to search.
    prevextensions: list of `~distutils.core.Extension` objects
        The extensions that are already defined.  Any .pyx files already here 
        will be ignored.
    extincludedirs : list of str or None
        Directories to include as the `include_dirs` argument to the generated
        `~distutils.core.Extension` objects.
    
    Returns
    -------
    exts : list of `~distutils.core.Extension` objects
        The new extensions that are needed to compile all .pyx files (does not
        include any already in `prevextensions`).
    """
    
    prevpyxpaths = []
    for ext in prevextensions:
        for s in ext.sources:
            if s.endswith('.pyx'):
                prevpyxpaths.append(os.path.realpath(s))
    
    ext_modules = []
    for extmod, pyxfn in iter_pyx_files(srcdir):
        if os.path.realpath(pyxfn) not in prevpyxpaths:
            ext_modules.append(Extension(extmod, [pyxfn],
                                         include_dirs=extincludedirs))
                
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


def adjust_compiler():
    """
    This function detects broken compilers and switches to another.  If
    the environment variable CC is explicitly set, or a compiler is
    specified on the commandline, no override is performed -- the purpose
    here is to only override a default compiler.

    The specific compilers with problems are:
        - The default compiler in XCode-4.2, llvm-gcc-4.2,
          segfaults when compiling wcslib.

    The set of broken compilers can be updated by changing the
    compiler_mapping variable.
    """
    if 'CC' in os.environ:
        return

    if get_distutils_option(
        'compiler', ['build', 'build_ext', 'build_clib']) is not None:
        return

    from distutils import ccompiler
    import subprocess

    compiler_mapping = {
        'i686-apple-darwin11-llvm-gcc-4.2': 'clang'
        }

    c = ccompiler.new_compiler()
    output = subprocess.check_output(c.compiler + ['--version'])
    version = output.split()[0]
    if version in compiler_mapping:
        os.environ['CC'] = compiler_mapping[version]


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
