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
import re
import shlex

from distutils import log
from distutils.dist import Distribution
from distutils.errors import DistutilsError
from distutils.core import Extension
from distutils.log import warn
from distutils.command.build import build as DistutilsBuild

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


class AstropyBuild(DistutilsBuild):
    """
    A custom 'build' command that adds the following options:

        `--disable-legacy` : Disable installing legacy shims
    """
    user_options = DistutilsBuild.user_options + [
        ('disable-legacy', None,
         "Don't install legacy shims")]

    boolean_options = DistutilsBuild.boolean_options + ['disable-legacy']

    def initialize_options(self):
        DistutilsBuild.initialize_options(self)
        # Set a default value for disable_legacy
        self.disable_legacy = False


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

        description = 'Build Sphinx documentation for Astropy environment'
        user_options = BuildDoc.user_options[:]
        user_options.append(('clean-docs', 'l', 'Completely clean previously '
                                                'builds, including auotmodapi-'
                                                'generated files before '
                                                'building new ones'))

        user_options.append(('no-intersphinx', 'n', 'Skip intersphinx, even if '
                                                    'conf.py says to use it'))
        boolean_options = BuildDoc.boolean_options[:]
        boolean_options.append('clean-docs')
        boolean_options.append('no-intersphinx')

        _self_iden_rex = re.compile(r"self\.([^\d\W][\w]+)", re.UNICODE)

        def initialize_options(self):
            BuildDoc.initialize_options(self)
            self.clean_docs = False
            self.no_intersphinx = False

        def finalize_options(self):
            from os.path import isdir
            from shutil import rmtree

            #Clear out previous sphinx builds, if requested
            if self.clean_docs:
                dirstorm = ['docs/_generated']
                if self.build_dir is None:
                    dirstorm.append('docs/_build')
                else:
                    dirstorm.append(self.build_dir)

                for d in dirstorm:
                    if isdir(d):
                        log.info('Cleaning directory ' + d)
                        rmtree(d)
                    else:
                        log.info('Not cleaning directory ' + d + ' because '
                                 'not present or not a directory')

            BuildDoc.finalize_options(self)

        def run(self):
            from os.path import split, join
            from distutils.cmd import DistutilsOptionError
            from subprocess import Popen, PIPE
            from textwrap import dedent
            from inspect import getsourcelines

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
                        'Attempted to build_sphinx in a location where' +
                        staticdir + 'is a file.  Must be a directory.')
                self.mkpath(staticdir)

            #Now make sure Astropy is built and determine where it was built
            build_cmd = self.reinitialize_command('build')
            build_cmd.inplace = 0
            self.run_command('build')
            build_cmd = self.get_finalized_command('build')
            build_cmd_path = os.path.abspath(build_cmd.build_lib)

            #Now generate the source for and spawn a new process that runs the
            #command.  This is needed to get the correct imports for the built
            #version

            runlines, runlineno = getsourcelines(BuildDoc.run)
            subproccode = dedent("""
            from sphinx.setup_command import *

            os.chdir('{srcdir}')
            sys.path.insert(0,'{build_cmd_path}')

            """).format(build_cmd_path=build_cmd_path, srcdir=self.source_dir)
            #runlines[1:] removes 'def run(self)' on the first line
            subproccode += dedent(''.join(runlines[1:]))

            # All "self.foo" in the subprocess code needs to be replaced by the
            # values taken from the current self in *this* process
            subproccode = AstropyBuildSphinx._self_iden_rex.split(subproccode)
            for i in range(1, len(subproccode), 2):
                iden = subproccode[i]
                val = getattr(self, iden)
                if iden.endswith('_dir'):
                    #Directories should be absolute, because the `chdir` call
                    #in the new process moves to a different directory
                    subproccode[i] = repr(os.path.abspath(val))
                else:
                    subproccode[i] = repr(val)
            subproccode = ''.join(subproccode)

            if self.no_intersphinx:
                #the confoverrides variable in sphinx.setup_command.BuildDoc can
                #be used to override the conf.py ... but this could well break
                #if future versions of sphinx change the internals of BuildDoc,
                #so remain vigilant!
                subproccode = subproccode.replace('confoverrides = {}',
                    'confoverrides = {\'intersphinx_mapping\':{}}')

            log.debug('Starting subprocess of {0} with python code:\n{1}\n'
                      '[CODE END])'.format(sys.executable, subproccode))

            proc = Popen([sys.executable], stdin=PIPE)
            proc.communicate(subproccode)
            if proc.returncode != 0:
                log.warn('Sphinx Documentation subprocess failed with return '
                         'code ' + str(proc.returncode))

except ImportError as e:
    if 'sphinx' in e.args[0]:  # Sphinx not present
        AstropyBuildSphinx = None
    else:
        raise


def get_distutils_display_options():
    """ Returns a set of all the distutils display options in their long and
    short forms.  These are the setup.py arguments such as --name or --version
    which print the project's metadata and then exit.

    Returns
    -------
    opts : set
        The long and short form display option arguments, including the - or --
    """

    short_display_opts = set('-' + o[1] for o in Distribution.display_options
                             if o[1])
    long_display_opts = set('--' + o[0] for o in Distribution.display_options)

    # Include -h and --help which are not explicitly listed in
    # Distribution.display_options (as they are handled by optparse)
    short_display_opts.add('-h')
    long_display_opts.add('--help')

    return short_display_opts.union(long_display_opts)


def is_distutils_display_option():
    """ Returns True if sys.argv contains any of the distutils display options
    such as --version or --name.
    """

    display_options = get_distutils_display_options()
    return bool(set(sys.argv[1:]).intersection(display_options))


cmdclassd = {}
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

    display_opts = get_distutils_display_options()
    args = [arg for arg in sys.argv[1:] if arg not in display_opts]

    # Pre-parse the Distutils command-line options and config files to
    # if the option is set.
    dist = Distribution({'script_name': os.path.basename(sys.argv[0]),
                         'script_args': args})
    dist.cmdclass.update(cmdclassd)

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

    if release or not HAVE_CYTHON:
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
                        raise IOError(
                            'Could not find C file {0} for Cython file {1} '
                            'when building extension {2}. '
                            'Cython must be installed to build from a git '
                            'checkout'.format(cfn, pyxfn, ext.name))

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

    from distutils import ccompiler, sysconfig
    import re

    compiler_mapping = [
        (b'i686-apple-darwin[0-9]*-llvm-gcc-4.2', 'clang')
        ]

    global _adjusted_compiler
    if _adjusted_compiler:
        return

    # Whatever the result of this function is, it only needs to be run once
    _adjusted_compiler = True

    if 'CC' in os.environ:

        # Check that CC is not set to llvm-gcc-4.2
        c_compiler = os.environ['CC']

        version = get_compiler_version(c_compiler)

        for broken, fixed in compiler_mapping:
            if re.match(broken, version):
                print("Compiler specified by CC environment variable ({0:s}: {1:s}) will fail to compile Astropy. Please set CC={2:s} and try again. You can do this for example by doing:\n\n    CC={2:s} python setup.py <command>\n\nwhere <command> is the command you ran.".format(c_compiler, version, fixed))
                sys.exit(1)

    if get_distutils_option(
        'compiler', ['build', 'build_ext', 'build_clib']) is not None:
        return

    compiler_type = ccompiler.get_default_compiler()

    if compiler_type == 'unix':

        # We have to get the compiler this way, as this is the one that is
        # used if os.environ['CC'] is not set. It is actually read in from
        # the Python Makefile. Note that this is not necessarily the same
        # compiler as returned by ccompiler.new_compiler()
        c_compiler = sysconfig.get_config_var('CC')

        version = get_compiler_version(c_compiler)

        for broken, fixed in compiler_mapping:
            if re.match(broken, version):
                os.environ['CC'] = fixed
                break

def get_compiler_version(compiler):

    import subprocess

    process = subprocess.Popen(
    shlex.split(compiler) + ['--version'], stdout=subprocess.PIPE)

    output = process.communicate()[0].strip()
    version = output.split()[0]

    return version


def is_in_build_mode():
    """
    Determines if the current package is being built.

    Returns
    -------
    buildmode : bool
        True if the current package is in the process of being built.

    See Also
    --------
    `set_build_mode`
    """
    #_ASTROPY_SETUP_ is added to the builtins in setup.py or astropy/__init__.py
    return _ASTROPY_SETUP_


def set_build_mode(val=True):
    """
    Sets whether or not the current package is being built.

    Parameters
    ----------
    val : bool
        Whether or not build mode should be activated.

    See Also
    --------
    `is_in_build_mode`
    """
    from sys import version_info

    if version_info[0] >= 3:
        import builtins
    else:
        import __builtin__ as builtins
    builtins._ASTROPY_SETUP_ = val


def setup_test_command(package_name):
    return type(package_name + '_test_command', (astropy_test,),
                {'package_name': package_name})


def import_file(filename):
    """
    Imports a module from a single file as if it doesn't belong to a
    particular package.
    """
    # Specifying a traditional dot-separated fully qualified name here
    # results in a number of "Parent module 'astropy' not found while
    # handling absolute import" warnings.  Using the same name, the
    # namespaces of the modules get merged together.  So, this
    # generates an underscore-separated name which is more likely to
    # be unique, and it doesn't really matter because the name isn't
    # used directly here anyway.
    with open(filename, 'U') as fd:
        name = '_'.join(
            os.path.relpath(os.path.splitext(filename)[0]).split(os.sep)[1:])
        return imp.load_module(name, fd, filename, ('.py', 'U', 1))


def get_legacy_alias_dir():
    return os.path.join('build', 'legacy-aliases')


legacy_shim_template = """
# This is generated code.  DO NOT EDIT!

from __future__ import absolute_import

# This implements a PEP 302 finder/loader pair that translates
# {old_package}.foo import {new_package}.foo.  This approach allows
# relative imports in astropy that go above the level of the
# {new_package} subpackage to work.
class Finder(object):
    def find_module(self, fullname, path=None):
        if fullname.startswith("{old_package}."):
            return self.Loader()

    class Loader(object):
        def load_module(self, fullname):
            import importlib
            fullname = fullname[len("{old_package}"):]
            return importlib.import_module(fullname, package="{new_package}")

import sys
sys.meta_path.append(Finder())
# Carefully clean up the namespace, since we can't use __all__ here
del sys
del Finder

import warnings
warnings.warn(
    "{old_package} is deprecated.  Use {new_package} instead.",
    DeprecationWarning)
del warnings

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

    # If legacy shims have been disabled at the commandline, simply do
    # nothing.
    if get_distutils_option('disable_legacy', ['build']):
        return (None, None)

    found_legacy_module = True
    try:
        location = imp.find_module(old_package)
    except ImportError:
        found_legacy_module = False
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
                    found_legacy_module = False

    shim_dir = os.path.join(get_legacy_alias_dir(), old_package)

    if found_legacy_module and not is_distutils_display_option():
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
