# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module contains a number of utilities for use during
setup/build/packaging that are useful to astropy as a whole.
"""

from __future__ import absolute_import

import imp
import os
import shutil
import sys

from distutils import log
from distutils.dist import Distribution
from distutils.errors import DistutilsError
from distutils.core import Extension
from distutils.log import warn

from .tests.helper import astropy_test


try:
    import Cython
    HAVE_CYTHON = True
except ImportError:
    HAVE_CYTHON = False

try:
    from numpy import get_include as get_numpy_include
    numpy_includes = get_numpy_include()
except ImportError:
    numpy_includes = []

try:
    from sphinx.setup_command import BuildDoc

    class AstropyBuildSphinx(BuildDoc):
        """ A version of the ``build_sphinx`` command that uses the
        version of Astropy that is built by the setup ``build`` command,
        rather than whatever is installed on the system - to build docs
        against the installed version, run ``make html`` in the
        ``astropy/docs`` directory.

        This also automatically creates the docs/_static directories -
        this is needed because github won't create the _static dir
        because it has no tracked files.
        """
        def run(self):
            from os.path import split, join
            from distutils.cmd import DistutilsOptionError

            # If possible, create the _static dir
            if self.build_dir is not None:
                # the _static dir should be in the same place as the _build dir
                # for Astropy
                basedir, subdir = split(self.build_dir)
                if subdir == '':  # the path has a trailing /...
                    basedir, subdir = split(basedir)
                staticdir = join(basedir, '_static')
                if os.path.isfile(staticdir):
                    raise DistutilsOptionError(
                        'Attempted to build_sphinx such in a location where' +
                        staticdir + 'is a file.  Must be a directory.')
                self.mkpath(staticdir)

            #Now make sure Astropy is built and inject it into the sphinx path
            build_cmd = self.reinitialize_command('build')
            build_cmd.inplace = 0
            self.run_command('build')
            build_cmd = self.get_finalized_command('build')
            new_path = os.path.abspath(build_cmd.build_lib)
            sys.path.insert(0, os.path.abspath(new_path))

            return BuildDoc.run(self)

except ImportError as e:
    if 'sphinx' in e.args[0]:  # Sphinx not present
        AstropyBuildSphinx = None
    else:
        raise


def get_distutils_option(option, commands):
    """ Returns the value of the given distutils option.

    Parameters
    ----------
    option : str
        The name of the option

    commands : list of str
        The list of commands on which this option is available

    Returns
    -------
    val : str or None
        the value of the given distutils option. If the option is not set,
        returns None.
    """
    # Pre-parse the Distutils command-line options and config files to
    # if the option is set.
    dist = Distribution({'script_name': os.path.basename(sys.argv[0]),
                         'script_args': sys.argv[1:]})
    try:
        dist.parse_config_files()
        dist.parse_command_line()
    except DistutilsError:
        # Let distutils handle this itself
        return None
    except AttributeError:
        # This seems to get thrown for ./setup.py --help
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


def get_compiler_option():
    """ Determines the compiler that will be used to build extension modules.

    Returns
    -------
    compiler : str
        The compiler option specificied for the build, build_ext, or build_clib
        command; or the default compiler for the platform if none was
        specified.

    """

    compiler = get_distutils_option('compiler',
                                    ['build', 'build_ext', 'build_clib'])
    if compiler is None:
        import distutils.ccompiler
        return distutils.ccompiler.get_default_compiler()

    return compiler


def get_debug_option():
    """ Determines if the build is in debug mode.

    Returns
    -------
    debug : bool
        True if the current build was started with the debug option, False
        otherwise.

    """

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


def update_package_files(srcdir, extensions, package_data, packagenames,
                         package_dirs):
    """ Extends existing extensions, package_data, packagenames and
    package_dirs collections by iterating through all packages in
    ``srcdir`` and locating a ``setup_package.py`` module.  This
    module can contain any of three functions: ``get_extensions()``,
    ``get_package_data()``, and ``get_legacy_alias()``.

    Each of those functions take no arguments.  ``get_extensions``
    returns a list of `distutils.extension.Extension` objects.
    ``get_package_data()`` returns a dict formatted as required by the
    ``package_data`` argument to ``setup()``.  ``get_legacy_alias()``
    should call `add_legacy_alias` and return its result.

    The purpose of this function is to allow subpackages to update the
    arguments to the package's ``setup()`` function in its setup.py
    script, rather than having to specify all extensions/package data
    directly in the setup.py.  It updates existing lists in the
    setup.py rather than returning new ones.  See Astropy's own
    ``setup.py`` for example usage and the Astropy development docs
    for more details.

    """

    from astropy.version import release

    # For each of the setup_package.py modules, extract any information that is
    # needed to install them.
    for setuppkg in iter_setup_packages(srcdir):
        # get_extensions must include any Cython extensions by their .pyx
        # filename.
        if hasattr(setuppkg, 'get_extensions'):
            extensions.extend(setuppkg.get_extensions())

        if hasattr(setuppkg, 'get_package_data'):
            package_data.update(setuppkg.get_package_data())
        if hasattr(setuppkg, 'get_legacy_alias'):
            pkg, dir = setuppkg.get_legacy_alias()
            if pkg is not None:
                packagenames.append(pkg)
                package_dirs[pkg] = dir

    # Locate any .pyx files not already specified, and add their extensions in.
    # The default include dirs include numpy to facilitate numerical work.
    extensions.extend(get_cython_extensions(srcdir, extensions,
                                            [numpy_includes]))

    # Now remove extensions that have the special name 'skip_cython', as they
    # exist Only to indicate that the cython extensions shouldn't be built
    for i, ext in reversed(list(enumerate(extensions))):
        if ext.name == 'skip_cython':
            del extensions[i]

    if not HAVE_CYTHON:
        # Replace .pyx with C-equivalents, unless c files are missing
        for idx, ext in reversed(list(enumerate(extensions))):
            for jdx, src in enumerate(ext.sources):
                if src.endswith('.pyx'):
                    pyxfn = src
                    cfn = src[:-4] + '.c'
                elif src.endswith('.c'):
                    pyxfn = src[:-2] + '.pyx'
                    cfn = src
                if os.path.isfile(pyxfn):
                    if os.path.isfile(cfn):
                        ext.sources[jdx] = cfn
                    else:
                        log.warn('Could not find c file {0} for {1}; skipping '
                                 'extension {2}.'.format(cfn, pyxfn, ext.name))
                        del extensions[idx]
                        break

    # On Microsoft compilers, we need to pass the '/MANIFEST'
    # commandline argument.  This was the default on MSVC 9.0, but is
    # now required on MSVC 10.0, but it doesn't seeem to hurt to add
    # it unconditionally.
    if get_compiler_option() == 'msvc':
        for ext in extensions:
            ext.extra_link_args.append('/MANIFEST')


def iter_setup_packages(srcdir):
    """ A generator that finds and imports all of the ``setup_package.py``
    modules in the source packages.

    Returns
    -------
    modgen : generator
        A generator that yields (modname, mod), where `mod` is the module and
        `modname` is the module name for the ``setup_package.py`` modules.

    """

    for root, dirs, files in os.walk(srcdir):
        if 'setup_package.py' in files:
            filename = os.path.join(root, 'setup_package.py')
            module = import_file(filename)
            yield module


def iter_pyx_files(srcdir):
    """ A generator that yields Cython source files (ending in '.pyx') in the
    source packages.

    Returns
    -------
    pyxgen : generator
        A generator that yields (extmod, fullfn) where `extmod` is the
        full name of the module that the .pyx file would live in based
        on the source directory structure, and `fullfn` is the path to
        the .pyx file.

    """
    for dirpath, dirnames, filenames in os.walk(srcdir):
        modbase = dirpath.replace(os.sep, '.')
        for fn in filenames:
            if fn.endswith('.pyx'):
                fullfn = os.path.join(dirpath, fn)
                # Package must match file name
                extmod = modbase + '.' + fn[:-4]
                yield (extmod, fullfn)


def get_cython_extensions(srcdir, prevextensions=tuple(), extincludedirs=None):
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
    """ Write `data` to `filename`, if the content of the file is different.

    Parameters
    ----------
    filename : str
        The file name to be written to.
    data : bytes
        The data to be written to `filename`.
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
    if (int(major), int(minor)) < (1, 4):
        msg = "numpy version 1.4 or later must be installed to build astropy"
        raise ImportError(msg)


def get_numpy_include_path():
    """
    Gets the path to the numpy headers.
    """
    import numpy

    try:
        numpy_include = numpy.get_include()
    except AttributeError:
        numpy_include = numpy.get_numpy_include()
    return numpy_include


_adjusted_compiler = False


def adjust_compiler():
    """
    This function detects broken compilers and switches to another.  If
    the environment variable CC is explicitly set, or a compiler is
    specified on the commandline, no override is performed -- the purpose
    here is to only override a default compiler.

    The specific compilers with problems are:

        * The default compiler in XCode-4.2, llvm-gcc-4.2,
          segfaults when compiling wcslib.

    The set of broken compilers can be updated by changing the
    compiler_mapping variable.  It is a list of 2-tuples where the
    first in the pair is a regular expression matching the version
    of the broken compiler, and the second is the compiler to change
    to.
    """

    global _adjusted_compiler
    if _adjusted_compiler:
        return

    # Whatever the result of this function is, it only needs to be run once
    _adjusted_compiler = True

    if 'CC' in os.environ:
        return

    if get_distutils_option(
        'compiler', ['build', 'build_ext', 'build_clib']) is not None:
        return

    from distutils import ccompiler
    import subprocess
    import re

    compiler_mapping = [
        (b'i686-apple-darwin[0-9]*-llvm-gcc-4.2', 'clang')
        ]

    c = ccompiler.new_compiler()
    # The MSVC ccompiler class doesn't have a `compiler` member.
    if not hasattr(c, 'compiler'):
        return
    process = subprocess.Popen(
        c.compiler + ['--version'], stdout=subprocess.PIPE)
    output = process.communicate()[0].strip()
    version = output.split()[0]
    for broken, fixed in compiler_mapping:
        if re.match(broken, version):
            os.environ['CC'] = fixed
            break


def is_in_build_mode():
    return __builtins__.get('_build_mode')


def set_build_mode():
    __builtins__['_build_mode'] = True


def setup_test_command(package_name):
    return type(package_name + '_test_command', (astropy_test,),
                {'package_name': package_name})


def import_file(filename):
    """
    Imports a module from a single file as if it doesn't belong to a
    particular package.
    """
    with open(filename, 'U') as fd:
        # If we specify a fully qualified name here, we'll get a
        # number of "Parent module 'astropy' not found while handling
        # absolute import" warnings.  Since this function is just used
        # to import setup_package.py files without importing their
        # parent packages, we can just give the name of the module.
        # In general, however, this is not something one would want to
        # do.
        name = os.path.splitext(os.path.basename(filename))[0]
        return imp.load_module(name, fd, filename, ('.py', 'U', 1))


def get_legacy_alias_dir():
    return os.path.join('build', 'legacy-aliases')


legacy_shim_template = """
# This is generated code.  DO NOT EDIT!

import warnings
warnings.warn(
    "{old_package} is deprecated.  Use {new_package} instead.",
    DeprecationWarning)

import pkgutil
__path__ = pkgutil.extend_path(__path__, "{new_package}")

from {new_package} import *
from astropy import __version__
__version__ = {equiv_version!r} + '-' + __version__
{extras}

_is_astropy_legacy_alias = True
"""


def add_legacy_alias(old_package, new_package, equiv_version, extras={}):
    """
    Adds a legacy alias that makes *pkgfrom* also importable as
    *pkgto*.

    For example::

       add_legacy_alias('astropy.io.vo', 'vo')

    If the legacy package is importable and it is not merely the
    compatibility shim, a warning is printed to the user, and the
    shim is not installed.

    Parameters
    ----------
    old_package : str
        The old namespace.  Must be a single name (i.e. not have `.`).

    new_package : str
        The new namespace, specified using `.` as a delimiter

    equiv_version : str
        The equivalent version of the old package.  Code using the
        legacy shim may do a version check, and this version should be
        based on the version of the legacy package, not the version of
        astropy.

    extras : dict
        A dictionary of extra values to include in the legacy shim template;
        the keys should be the variable names, while the values will be written
        to the template in their repr() form, so they should generally be
        simple objects such as strings.

    Returns
    -------
    old_package, shim_dir : (str, str)
        The name of the alias package and its source directory in the
        file system (useful for adding to distutils' `package_dir` kwarg.
    """
    import imp

    found_legacy_module = False
    try:
        location = imp.find_module(old_package)
    except ImportError:
        pass
    else:
        # We want ImportError to raise here, because that means it was
        # found, but something else went wrong.

        # We could import the module here to determine if its "real"
        # or just a legacy alias.  However, importing the legacy alias
        # may cause importing of code within the astropy source tree,
        # which may require 2to3 to have been run.  It's safer to just
        # open the file and search for a string.
        filename = os.path.join(location[1], '__init__.py')
        if os.path.exists(filename):
            with open(filename, 'U') as fd:
                if '_is_astropy_legacy_alias' in fd.read():
                    found_legacy_module = True

    shim_dir = os.path.join(get_legacy_alias_dir(), old_package)

    if found_legacy_module:
        warn('-' * 60)
        warn("The legacy package '{0}' was found.".format(old_package))
        warn("To install astropy's compatibility layer instead, uninstall")
        warn("'{0}' and then reinstall astropy.".format(old_package))
        warn('-' * 60)

        if os.path.isdir(shim_dir):
            shutil.rmtree(shim_dir)
        return (None, None)

    if extras:
        extras = '\n'.join('{0} = {1!r}'.format(*v) for v in extras.items())
    else:
        extras = ''

    if not os.path.isdir(shim_dir):
        os.makedirs(shim_dir)
    content = legacy_shim_template.format(**locals()).encode('utf-8')
    write_if_different(
        os.path.join(shim_dir, '__init__.py'), content)

    return (old_package, shim_dir)
