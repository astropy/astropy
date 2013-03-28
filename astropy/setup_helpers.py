# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module contains a number of utilities for use during
setup/build/packaging that are useful to astropy as a whole.
"""

from __future__ import absolute_import, print_function

import collections
import errno
import imp
import os
import re
import shlex
import shutil
import subprocess
import sys
import textwrap

from distutils import log
from distutils.dist import Distribution
from distutils.errors import DistutilsError, DistutilsFileError
from distutils.core import Extension
from distutils.core import Command
from distutils.command.sdist import sdist as DistutilsSdist
from setuptools.command.build_ext import build_ext as SetuptoolsBuildExt
from setuptools.command.build_py import build_py as SetuptoolsBuildPy

from setuptools.command.register import register as SetuptoolsRegister

from .tests.helper import astropy_test
from .utils import silence
from .utils.misc import walk_skip_hidden


try:
    import Cython
    HAVE_CYTHON = True
except ImportError:
    HAVE_CYTHON = False


try:
    import sphinx
    from sphinx.setup_command import BuildDoc as SphinxBuildDoc
    HAVE_SPHINX = True
except ImportError:
    HAVE_SPHINX = False


PY3 = sys.version_info[0] >= 3


_adjusted_compiler = False
def adjust_compiler(package):
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

        try:
            version = get_compiler_version(c_compiler)
        except OSError:
            msg = textwrap.dedent(
                    """
                    The C compiler set by the CC environment variable:

                        {compiler:s}

                    cannot be found or executed.
                    """.format(compiler=c_compiler))
            log.warn(msg)
            sys.exit(1)

        for broken, fixed in compiler_mapping:
            if re.match(broken, version):
                msg = textwrap.dedent(
                    """Compiler specified by CC environment variable
                    ({compiler:s}:{version:s}) will fail to compile {pkg:s}.
                    Please set CC={fixed:s} and try again.
                    You can do this, for example, by running:

                        CC={fixed:s} python setup.py <command>

                    where <command> is the command you ran.
                    """.format(compiler=c_compiler, version=version,
                               pkg=package, fixed=fixed))
                log.warn(msg)
                sys.exit(1)

        # If C compiler is set via CC, and isn't broken, we are good to go. We
        # should definitely not try accessing the compiler specified by
        # ``sysconfig.get_config_var('CC')`` lower down, because this may fail
        # if the compiler used to compile Python is missing (and maybe this is
        # why the user is setting CC). For example, the official Python 2.7.3
        # MacOS X binary was compled with gcc-4.2, which is no longer available
        # in XCode 4.
        return

    if get_distutils_build_option('compiler'):
        return

    compiler_type = ccompiler.get_default_compiler()

    if compiler_type == 'unix':

        # We have to get the compiler this way, as this is the one that is
        # used if os.environ['CC'] is not set. It is actually read in from
        # the Python Makefile. Note that this is not necessarily the same
        # compiler as returned by ccompiler.new_compiler()
        c_compiler = sysconfig.get_config_var('CC')

        try:
            version = get_compiler_version(c_compiler)
        except OSError:
            msg = textwrap.dedent(
                    """
                    The C compiler used to compile Python {compiler:s}, and
                    which is normally used to compile C extensions, is not
                    available. You can explicitly specifiy which compiler to
                    use by setting the CC environment variable, for example:

                        CC=gcc python setup.py <command>

                    or if you are using MacOS X, you can try:

                        CC=clang python setup.py <command>
                    """.format(compiler=c_compiler))
            log.warn(msg)
            sys.exit(1)


        for broken, fixed in compiler_mapping:
            if re.match(broken, version):
                os.environ['CC'] = fixed
                break


def get_compiler_version(compiler):

    process = subprocess.Popen(
    shlex.split(compiler) + ['--version'], stdout=subprocess.PIPE)

    output = process.communicate()[0].strip()
    version = output.split()[0]

    return version


def get_dummy_distribution():
    """Returns a distutils Distribution object used to instrument the setup
    environment before calling the actual setup() function.
    """

    global _registered_commands

    if _registered_commands is None:
        raise RuntimeError('astropy.setup_helpers.register_commands() must be '
                           'called before using '
                           'astropy.setup_helpers.get_dummy_distribution()')

    # Pre-parse the Distutils command-line options and config files to if
    # the option is set.
    dist = Distribution({'script_name': os.path.basename(sys.argv[0]),
                         'script_args': sys.argv[1:]})
    dist.cmdclass.update(_registered_commands)

    with silence():
        try:
            dist.parse_config_files()
            dist.parse_command_line()
        except (DistutilsError, AttributeError, SystemExit):
            # Let distutils handle DistutilsErrors itself AttributeErrors can
            # get raise for ./setup.py --help SystemExit can be raised if a
            # display option was used, for example
            pass

    return dist


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

    dist = get_dummy_distribution()

    for cmd in commands:
        cmd_opts = dist.command_options.get(cmd)
        if cmd_opts is not None and option in cmd_opts:
            return cmd_opts[option][1]
    else:
        return None


def get_distutils_build_option(option):
    """ Returns the value of the given distutils build option.

    Parameters
    ----------
    option : str
        The name of the option

    Returns
    -------
    val : str or None
        The value of the given distutils build option. If the option
        is not set, returns None.
    """
    return get_distutils_option(option, ['build', 'build_ext', 'build_clib'])


def get_distutils_install_option(option):
    """ Returns the value of the given distutils install option.

    Parameters
    ----------
    option : str
        The name of the option

    Returns
    -------
    val : str or None
        The value of the given distutils build option. If the option
        is not set, returns None.
    """
    return get_distutils_option(option, ['install'])


def get_distutils_build_or_install_option(option):
    """ Returns the value of the given distutils build or install option.

    Parameters
    ----------
    option : str
        The name of the option

    Returns
    -------
    val : str or None
        The value of the given distutils build or install option. If the
        option is not set, returns None.
    """
    return get_distutils_option(option, ['build', 'build_ext', 'build_clib',
                                         'install'])


def get_compiler_option():
    """ Determines the compiler that will be used to build extension modules.

    Returns
    -------
    compiler : str
        The compiler option specificied for the build, build_ext, or build_clib
        command; or the default compiler for the platform if none was
        specified.

    """

    compiler = get_distutils_build_option('compiler')
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

    try:
        from .version import debug as current_debug
    except ImportError:
        current_debug = None

    # Only modify the debug flag if one of the build commands was explicitly
    # run (i.e. not as a sub-command of something else)
    dist = get_dummy_distribution()
    if any(cmd in dist.commands for cmd in ['build', 'build_ext']):
        debug = bool(get_distutils_build_option('debug'))
    else:
        debug = bool(current_debug)

    if current_debug is not None and current_debug != debug:
        build_ext_cmd = dist.get_command_class('build_ext')
        build_ext_cmd.force_rebuild = True

    return debug


_registered_commands = None
def register_commands(package, version, release):
    global _registered_commands

    if _registered_commands is not None:
        return _registered_commands

    _registered_commands = {
        'test': generate_test_command(package),

         # Use distutils' sdist because it respects package_data.
         # setuptools/distributes sdist requires duplication of information in
         # MANIFEST.in
         'sdist': DistutilsSdist,

         # The exact form of the build_ext command depends on whether or not
         # we're building a release version
         'build_ext': generate_build_ext_command(release),

         # We have a custom build_py to generate the default configuration file
         'build_py': AstropyBuildPy,

         'register': AstropyRegister
    }

    try:
        import bdist_mpkg
    except ImportError:
        pass
    else:
        # Use a custom command to build a dmg (on MacOS X)
        _registered_commands['bdist_dmg'] = bdist_dmg


    if HAVE_SPHINX:
        _registered_commands['build_sphinx'] = AstropyBuildSphinx

    # Need to override the __name__ here so that the commandline options are
    # presented as being related to the "build" command, for example; normally
    # this wouldn't be necessary since commands also have a command_name
    # attribute, but there is a bug in distutils' help display code that it
    # uses __name__ instead of command_name. Yay distutils!
    for name, cls in _registered_commands.items():
        cls.__name__ = name

    # Add a few custom options; more of these can be added by specific packages
    # later
    for option in [
            ('enable-legacy',
             "Install legacy shims", True),
            ('use-system-libraries',
             "Use system libraries whenever possible", True)]:
        add_command_option('build', *option)
        add_command_option('install', *option)

    return _registered_commands


def generate_test_command(package_name):
    return type(package_name + '_test_command', (astropy_test,),
                {'package_name': package_name})


def generate_build_ext_command(release):
    """
    Creates a custom 'build_ext' command that allows for manipulating some of
    the C extension options at build time.  We use a function to build the
    class since the base class for build_ext may be different depending on
    certain build-time parameters (for example, we may use Cython's build_ext
    instead of the default version in distutils).

    Uses the default distutils.command.build_ext by default.
    """

    uses_cython = should_build_with_cython(release)

    if uses_cython:
        from Cython.Distutils import build_ext as basecls
    else:
        basecls = SetuptoolsBuildExt

    attrs = dict(basecls.__dict__)
    orig_run = getattr(basecls, 'run', None)
    orig_finalize = getattr(basecls, 'finalize_options', None)

    def finalize_options(self):
        if orig_finalize is not None:
            orig_finalize(self)

        # Generate
        if self.uses_cython:
            try:
                from Cython import __version__ as cython_version
            except ImportError:
                # This shouldn't happen if we made it this far
                cython_version = None

            if (cython_version is not None and
                    cython_version != self.uses_cython):
                self.force_rebuild = True
                # Update the used cython version
                self.uses_cython = cython_version

        # Regardless of the value of the '--force' option, force a rebuild if
        # the debug flag changed from the last build
        if self.force_rebuild:
            self.force = True

    def run(self):
        # For extensions that require 'numpy' in their include dirs, replace
        # 'numpy' with the actual paths
        np_include = get_numpy_include_path()
        for extension in self.extensions:
            if 'numpy' in extension.include_dirs:
                idx = extension.include_dirs.index('numpy')
                extension.include_dirs.insert(idx, np_include)
                extension.include_dirs.remove('numpy')

            # Replace .pyx with C-equivalents, unless c files are missing
            for jdx, src in enumerate(extension.sources):
                if src.endswith('.pyx'):
                    pyxfn = src
                    cfn = src[:-4] + '.c'
                elif src.endswith('.c'):
                    pyxfn = src[:-2] + '.pyx'
                    cfn = src

                if os.path.isfile(pyxfn):
                    if self.uses_cython:
                        extension.sources[jdx] = pyxfn
                    else:
                        if os.path.isfile(cfn):
                            extension.sources[jdx] = cfn
                        else:
                            msg = (
                                'Could not find C file {0} for Cython file '
                                '{1} when building extension {2}. '
                                'Cython must be installed to build from a '
                                'git checkout'.format(cfn, pyxfn,
                                                      extension.name))
                            raise IOError(errno.ENOENT, msg, cfn)

        if orig_run is not None:
            # This should always be the case for a correctly implemented
            # distutils command.
            orig_run(self)

        # Update cython_version.py if building with Cython
        try:
            from .version import cython_version
        except ImportError:
            cython_version = 'unknown'
        if self.uses_cython and self.uses_cython != cython_version:
            astropy_dir = os.path.relpath(os.path.dirname(__file__))
            cython_py = os.path.join(astropy_dir, 'cython_version.py')
            with open(cython_py, 'w') as f:
                f.write('# Generated file; do not modify\n')
                f.write('cython_version = {0!r}\n'.format(self.uses_cython))

            if os.path.isdir(self.build_lib):
                # The build/lib directory may not exist if the build_py command
                # was not previously run, which may sometimes be the case
                self.copy_file(cython_py,
                               os.path.join(self.build_lib, cython_py),
                               preserve_mode=False)


        if not self.distribution.is_pure() and os.path.isdir(self.build_lib):
            # Finally, generate the default astropy.cfg; this can only be done
            # after extension modules are built as some extension modules
            # include config items.  We only do this if it's not pure python,
            # though, because if it is, we already did it in build_py
            default_cfg = generate_default_config(
                    os.path.abspath(self.build_lib),
                    self.distribution.packages[0])
            if default_cfg:
                default_cfg = os.path.relpath(default_cfg)
                self.copy_file(default_cfg,
                               os.path.join(self.build_lib, default_cfg),
                               preserve_mode=False)

    attrs['run'] = run
    attrs['finalize_options'] = finalize_options
    attrs['force_rebuild'] = False
    attrs['uses_cython'] = uses_cython

    return type('build_ext', (basecls, object), attrs)


class AstropyBuildPy(SetuptoolsBuildPy):

    def run(self):
        # first run the normal build_py
        SetuptoolsBuildPy.run(self)

        if self.distribution.is_pure():
            # Generate the default astropy.cfg - we only do this here if it's
            # pure python.  Otherwise, it'll happen at the end of build_exp
            default_cfg = generate_default_config(
                    os.path.abspath(self.build_lib),
                    self.distribution.packages[0])
            if default_cfg:
                default_cfg = os.path.relpath(default_cfg)
                self.copy_file(default_cfg,
                               os.path.join(self.build_lib, default_cfg),
                               preserve_mode=False)


def generate_default_config(build_lib, package):
    config_path = os.path.relpath(package)
    filename = os.path.join(config_path, package + '.cfg')

    if os.path.exists(filename):
        log.info('regenerating default {0}.cfg file'.format(package))
    else:
        log.info('generating default {0}.cfg file'.format(package))

    if PY3:
        builtins = 'builtins'
    else:
        builtins = '__builtin__'

    subproccode = (
        'import {builtins};{builtins}._ASTROPY_SETUP_ = True;'
        'from astropy.config.configuration import generate_all_config_items;'
        'generate_all_config_items({pkgnm!r}, True, filename={filenm!r})')
    subproccode = subproccode.format(builtins=builtins,
                                     pkgnm=package,
                                     filenm=os.path.abspath(filename))

    # Note that cwd=build_lib--we're importing astropy from the build/ dir
    # but using the astropy/ source dir as the config directory
    proc = subprocess.Popen([sys.executable, '-c', subproccode],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                             cwd=build_lib)
    stdout, stderr = proc.communicate()

    if proc.returncode == 0 and os.path.exists(filename):
        return filename
    else:
        msg = ('Generation of default configuration item failed! Stdout '
               'and stderr are shown below.\n'
               'Stdout:\n{stdout}\nStderr:\n{stderr}').decode('UTF-8')
        log.error(msg.format(stdout=stdout.decode('UTF-8'),
                             stderr=stderr.decode('UTF-8')))


def add_command_option(command, name, doc, is_bool=False):
    """
    Add a custom option to a setup command.

    Issues a warning if the option already exists on that command.

    Parameters
    ----------
    command : str
        The name of the command as given on the command line

    name : str
        The name of the build option

    doc : str
        A short description of the option, for the `--help` message

    is_bool : bool, optional
        When `True`, the option is a boolean option and doesn't
        require an associated value.
    """

    dist = get_dummy_distribution()
    cmdcls = dist.get_command_class(command)

    attr = name.replace('-', '_')

    if hasattr(cmdcls, attr):
        raise RuntimeError(
            '{0!r} already has a {1!r} class attribute, barring {2!r} from '
            'being usable as a custom option name.'.format(cmdcls, attr, name))

    for idx, cmd in enumerate(cmdcls.user_options):
        if cmd[0] == name:
            log.warning('Overriding existing {0!r} option '
                        '{1!r}'.format(command, name))
            del cmdcls.user_options[idx]
            if name in cmdcls.boolean_options:
                cmdcls.boolean_options.remove(name)
            break

    cmdcls.user_options.append((name, None, doc))

    if is_bool:
        cmdcls.boolean_options.append(name)

    # Distutils' command parsing requires that a command object have an
    # attribute with the same name as the option (with '-' replaced with '_')
    # in order for that option to be recognized as valid
    setattr(cmdcls, attr, None)


class AstropyRegister(SetuptoolsRegister):
    """Extends the built in 'register' command to support a ``--hidden`` option
    to make the registered version hidden on PyPI by default.

    The result of this is that when a version is registered as "hidden" it can
    still be downloaded from PyPI, but it does not show up in the list of
    actively supported versions under http://pypi.python.org/pypi/astropy, and
    is not set as the most recent version.

    Although this can always be set through the web interface it may be more
    convenient to be able to specify via the 'register' command.  Hidden may
    also be considered a safer default when running the 'register' command,
    though this command uses distutils' normal behavior if the ``--hidden``
    option is omitted.
    """

    user_options = SetuptoolsRegister.user_options + [
        ('hidden', None, 'mark this release as hidden on PyPI by default')
    ]
    boolean_options = SetuptoolsRegister.boolean_options + ['hidden']

    def initialize_options(self):
        SetuptoolsRegister.initialize_options(self)
        self.hidden = False

    def build_post_data(self, action):
        data = SetuptoolsRegister.build_post_data(self, action)
        if action == 'submit' and self.hidden:
            data['_pypi_hidden'] = '1'
        return data

    def _set_config(self):
        # The original register command is buggy--if you use .pypirc with a
        # server-login section *at all* the repository you specify with the -r
        # option will be overwritten with either the repository in .pypirc or
        # with the default,
        # If you do not have a .pypirc using the -r option will just crash.
        # Way to go distutils

        # If we don't set self.repository back to a default value _set_config
        # can crash if there was a user-supplied value for this option; don't
        # worry, we'll get the real value back afterwards
        self.repository = 'pypi'
        SetuptoolsRegister._set_config(self)
        options = self.distribution.get_option_dict('register')
        if 'repository' in options:
            source, value = options['repository']
            # Really anything that came from setup.cfg or the command line
            # should override whatever was in .pypirc
            self.repository = value


if HAVE_SPHINX:
    class AstropyBuildSphinx(SphinxBuildDoc):
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
        user_options = SphinxBuildDoc.user_options[:]
        user_options.append(('warnings-returncode', 'w',
                             'Parses the sphinx output and sets the return '
                             'code to 1 if there are any warnings. Note that '
                             'this will cause the sphinx log to only update '
                             'when it completes, rather than continuously as '
                             'is normally the case.'))
        user_options.append(('clean-docs', 'l',
                             'Completely clean previous builds, including '
                             'automodapi-generated files before building new '
                             'ones'))
        user_options.append(('no-intersphinx', 'n',
                             'Skip intersphinx, even if conf.py says to use '
                             'it'))
        user_options.append(('open-docs-in-browser', 'o',
                             'Open the docs in a browser (using the '
                             'webbrowser module) if the build finishes '
                             'successfully.'))

        boolean_options = SphinxBuildDoc.boolean_options[:]
        boolean_options.append('warnings-returncode')
        boolean_options.append('clean-docs')
        boolean_options.append('no-intersphinx')
        boolean_options.append('open-docs-in-browser')

        _self_iden_rex = re.compile(r"self\.([^\d\W][\w]+)", re.UNICODE)

        def initialize_options(self):
            SphinxBuildDoc.initialize_options(self)
            self.clean_docs = False
            self.no_intersphinx = False
            self.open_docs_in_browser = False
            self.warnings_returncode = False

        def finalize_options(self):
            #Clear out previous sphinx builds, if requested
            if self.clean_docs:
                dirstorm = ['docs/_generated']
                if self.build_dir is None:
                    dirstorm.append('docs/_build')
                else:
                    dirstorm.append(self.build_dir)

                for d in dirstorm:
                    if os.path.isdir(d):
                        log.info('Cleaning directory ' + d)
                        shutil.rmtree(d)
                    else:
                        log.info('Not cleaning directory ' + d + ' because '
                                 'not present or not a directory')

            SphinxBuildDoc.finalize_options(self)

        def run(self):
            import atexit
            import webbrowser

            from os.path import split, join, abspath
            from distutils.cmd import DistutilsOptionError
            from subprocess import Popen, PIPE, STDOUT
            from inspect import getsourcelines
            from urllib import pathname2url

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

            runlines, runlineno = getsourcelines(SphinxBuildDoc.run)
            subproccode = textwrap.dedent("""
            from sphinx.setup_command import *

            os.chdir('{srcdir}')
            sys.path.insert(0,'{build_cmd_path}')

            """).format(build_cmd_path=build_cmd_path, srcdir=self.source_dir)
            #runlines[1:] removes 'def run(self)' on the first line
            subproccode += textwrap.dedent(''.join(runlines[1:]))

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

            # To return the number of warnings, we need to capture stdout. This
            # prevents a continuous updating at the terminal, but there's no
            # apparent way around this.
            if self.warnings_returncode:
                proc = Popen([sys.executable], stdin=PIPE, stdout=PIPE, stderr=STDOUT)
                stdo, stde = proc.communicate(subproccode)

                print(stdo)

                stdolines = stdo.split('\n')

                if stdolines[-2] == 'build succeeded.':
                    retcode = 0
                else:
                    retcode = 1

                if retcode != 0:
                    def overrideexitcode():
                        raise SystemExit(retcode)
                    atexit.register(overrideexitcode)

            else:
                proc = Popen([sys.executable], stdin=PIPE)
                proc.communicate(subproccode)

            if proc.returncode == 0:
                if self.open_docs_in_browser:
                    if self.builder == 'html':
                        absdir = abspath(self.builder_target_dir)
                        fileurl = 'file://' + pathname2url(join(absdir, 'index.html'))
                        webbrowser.open(fileurl)
                    else:
                        log.warn('open-docs-in-browser option was given, but '
                                 'the builder is not html! Ignogring.')
            else:
                log.warn('Sphinx Documentation subprocess failed with return '
                         'code ' + str(proc.returncode))


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


def update_package_files(srcdir, extensions, package_data, packagenames,
                         package_dirs):
    """ Extends existing extensions, package_data, packagenames and
    package_dirs collections by iterating through all packages in
    ``srcdir`` and locating a ``setup_package.py`` module.  This
    module can contain the following functions: ``get_extensions()``,
    ``get_package_data()``, ``get_legacy_alias()``,
    ``get_build_options()``, and ``get_external_libraries()``.

    Each of those functions take no arguments.  ``get_extensions``
    returns a list of `distutils.extension.Extension` objects.
    ``get_package_data()`` returns a dict formatted as required by the
    ``package_data`` argument to ``setup()``.  ``get_legacy_alias()``
    should call `add_legacy_alias` and return its result.
    ``get_build_options()`` returns a list of tuples describing the
    extra build options to add.  ``get_external_libraries()`` returns
    a list of libraries that can optionally be built using external
    dependencies.

    The purpose of this function is to allow subpackages to update the
    arguments to the package's ``setup()`` function in its setup.py
    script, rather than having to specify all extensions/package data
    directly in the setup.py.  It updates existing lists in the
    setup.py rather than returning new ones.  See Astropy's own
    ``setup.py`` for example usage and the Astropy development docs
    for more details.
    """

    # For each of the setup_package.py modules, extract any
    # information that is needed to install them.  The build options
    # are extracted first, so that their values will be available in
    # subsequent calls to `get_extensions`, etc.
    for setuppkg in iter_setup_packages(srcdir):
        if hasattr(setuppkg, 'get_build_options'):
            options = setuppkg.get_build_options()
            for option in options:
                add_command_option('build', *option)
        if hasattr(setuppkg, 'get_external_libraries'):
            libraries = setuppkg.get_external_libraries()
            for library in libraries:
                add_external_library(library)

    # Check if all the legacy packages are needed
    if get_distutils_build_or_install_option('enable_legacy'):
        installed = []
        for setuppkg in iter_setup_packages(srcdir):
            if hasattr(setuppkg, 'get_legacy_alias'):
                pkg, dir = setuppkg.get_legacy_alias()
                if dir is None:
                    installed.append(pkg)
                else:
                    packagenames.append(pkg)
                    package_dirs[pkg] = dir
        if len(installed) > 0:
            lines = [
                '-' * 60,
                'The compatibility packages cannot be installed because the',
                'following legacy packages are already installed:']
            for pkg in installed:
                lines.append("    * {0:s}".format(pkg))
            lines.extend([
                '',
                'The compatibility packages can only be installed if none of',
                'the corresponding legacy packages are present.',
                '-' * 60
                ])
            log.warn('\n'.join(lines))
            sys.exit(1)

    for setuppkg in iter_setup_packages(srcdir):
        # get_extensions must include any Cython extensions by their .pyx
        # filename.
        if hasattr(setuppkg, 'get_extensions'):
            extensions.extend(setuppkg.get_extensions())
        if hasattr(setuppkg, 'get_package_data'):
            package_data.update(setuppkg.get_package_data())

    # Locate any .pyx files not already specified, and add their extensions in.
    # The default include dirs include numpy to facilitate numerical work.
    extensions.extend(get_cython_extensions(srcdir, extensions, ['numpy']))

    # Now remove extensions that have the special name 'skip_cython', as they
    # exist Only to indicate that the cython extensions shouldn't be built
    for i, ext in reversed(list(enumerate(extensions))):
        if ext.name == 'skip_cython':
            del extensions[i]

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
    for root, dirs, files in walk_skip_hidden(srcdir):
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
    for dirpath, dirnames, filenames in walk_skip_hidden(srcdir):
        modbase = dirpath.replace(os.sep, '.')
        for fn in filenames:
            if fn.endswith('.pyx'):
                fullfn = os.path.join(dirpath, fn)
                # Package must match file name
                extmod = modbase + '.' + fn[:-4]
                yield (extmod, fullfn)


def should_build_with_cython(release=None):
    """Returns the previously used Cython version (or 'unknown' if not
    previously built) if Cython should be used to build extension modules from
    pyx files.  If the ``release`` parameter is not specified an attempt is
    made to determine the release flag from `astropy.version`.
    """

    if release is None:
        try:
            from .version import release
        except ImportError:
            pass

    try:
        from .version import cython_version
    except ImportError:
        cython_version = 'unknown'

    # Only build with Cython if, of course, Cython is installed, we're in a
    # development version (i.e. not release) or the Cython-generated source
    # files haven't been created yet (cython_version == 'unknown'). The latter
    # case can happen even when release is True if checking out a release tag
    # from the repository
    if HAVE_CYTHON and (not release or cython_version == 'unknown'):
        return cython_version
    else:
        return False


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

    # Vanilla setuptools and old versions of distribute include Cython files
    # as .c files in the sources, not .pyx, so we cannot simply look for
    # existing .pyx sources in the previous sources, but we should also check
    # for .c files with the same remaining filename. So we look for .pyx and
    # .c files, and we strip the extension.

    prevsourcepaths = []
    for ext in prevextensions:
        for s in ext.sources:
            if s.endswith(('.pyx', '.c')):
                prevsourcepaths.append(os.path.realpath(os.path.splitext(s)[0]))

    ext_modules = []
    for extmod, pyxfn in iter_pyx_files(srcdir):
        if os.path.realpath(os.path.splitext(pyxfn)[0]) not in prevsourcepaths:
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

    requirement_met = False

    try:
        import numpy
    except ImportError:
        pass
    else:
        major, minor, rest = numpy.__version__.split(".", 2)
        requirement_met = (int(major), int(minor)) >= (1, 4)

    if not requirement_met:
        msg = "numpy version 1.4 or later must be installed to build astropy"
        raise ImportError(msg)

    return numpy


def get_numpy_include_path():
    """
    Gets the path to the numpy headers.
    """

    numpy = check_numpy()

    try:
        numpy_include = numpy.get_include()
    except AttributeError:
        numpy_include = numpy.get_numpy_include()
    return numpy_include


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

       add_legacy_alias('astropy.io.votable', 'vo')

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

    # If legacy shims have not been enabled at the commandline, simply do
    # nothing.
    if not get_distutils_build_or_install_option('enable_legacy'):
        return (old_package, None)

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
        if os.path.isdir(shim_dir):
            shutil.rmtree(shim_dir)
        return (old_package, None)

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


class DistutilsExtensionArgs(collections.defaultdict):
    """
    A special dictionary whose default values are the empty list.

    This is useful for building up a set of arguments for
    `distutils.Extension` without worrying whether the entry is
    already present.
    """
    def __init__(self, *args, **kwargs):
        def default_factory():
            return []

        super(DistutilsExtensionArgs, self).__init__(
            default_factory, *args, **kwargs)

    def update(self, other):
        for key, val in other.items():
            self[key].extend(val)


def pkg_config(packages, default_libraries):
    """
    Uses pkg-config to update a set of distutils Extension arguments
    to include the flags necessary to link against the given packages.

    If the pkg-config lookup fails, default_libraries is applied to
    libraries.

    Parameters
    ----------
    packages : list of str
        A list of pkg-config packages to look up.

    default_libraries : list of str
        A list of library names to use if the pkg-config lookup fails.

    Returns
    -------
    config : dict
        A dictionary containing keyword arguments to
        `distutils.Extension`.  These entries include:

        - ``include_dirs``: A list of include directories
        - ``library_dirs``: A list of library directories
        - ``libraries``: A list of libraries
        - ``define_macros``: A list of macro defines
        - ``undef_macros``: A list of macros to undefine
        - ``extra_compile_args``: A list of extra arguments to pass to
          the compiler
    """

    flag_map = {'-I': 'include_dirs', '-L': 'library_dirs', '-l': 'libraries',
                '-D': 'define_macros', '-U': 'undef_macros'}
    command = "pkg-config --libs --cflags {0}".format(' '.join(packages)),

    result = DistutilsExtensionArgs()

    try:
        pipe = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
        output = pipe.communicate()[0].strip()
    except subprocess.CalledProcessError as e:
        lines = [
            "pkg-config failed.  This may cause the build to fail below.",
            "  command: {0}".format(e.cmd),
            "  returncode: {0}".format(e.returncode),
            "  output: {0}".format(e.output)
            ]
        log.warn('\n'.join(lines))
        result['libraries'].extend(default_libraries)
    else:
        if pipe.returncode != 0:
            lines = [
                "pkg-config could not lookup up package(s) {0}.".format(
                    ", ".join(packages)),
                "This may cause the build to fail below."
                ]
            log.warn('\n'.join(lines))
            result['libraries'].extend(default_libraries)
        else:
            for token in output.split():
                # It's not clear what encoding the output of
                # pkg-config will come to us in.  It will probably be
                # some combination of pure ASCII (for the compiler
                # flags) and the filesystem encoding (for any argument
                # that includes directories or filenames), but this is
                # just conjecture, as the pkg-config documentation
                # doesn't seem to address it.
                arg = token[:2].decode('ascii')
                value = token[2:].decode(sys.getfilesystemencoding())
                if arg in flag_map:
                    if arg == '-D':
                        value = tuple(value.split('=', 1))
                    result[flag_map[arg]].append(value)
                else:
                    result['extra_compile_args'].append(value)

    return result


def add_external_library(library):
    """
    Add a build option for selecting the internal or system copy of a library.

    Parameters
    ----------
    library : str
        The name of the library.  If the library is `foo`, the build
        option will be called `--use-system-foo`.
    """

    for command in ['build', 'build_ext', 'install']:
        add_command_option(command, 'use-system-' + library,
                           'Use the system {0} library'.format(library),
                           is_bool=True)


def use_system_library(library):
    """
    Returns `True` if the build configuration indicates that the given
    library should use the system copy of the library rather than the
    internal one.

    For the given library `foo`, this will be `True` if
    `--use-system-foo` or `--use-system-libraries` was provided at the
    commandline or in `setup.cfg`.

    Parameters
    ----------
    library : str
        The name of the library

    Returns
    -------
    use_system : bool
        `True` if the build should use the system copy of the library.
    """
    return (
        get_distutils_build_or_install_option('use_system_{0}'.format(library))
        or get_distutils_build_or_install_option('use_system_libraries'))


def filter_packages(packagenames):
    """
    Removes some packages from the package list that shouldn't be
    installed on the current version of Python.
    """

    if PY3:
        exclude = '_py2'
    else:
        exclude = '_py3'

    return [x for x in packagenames if not x.endswith(exclude)]


class bdist_dmg(Command):
    """
    The bdist_dmg command is used to produce the disk image containing the
    installer, and with a custom background and icon placement.
    """

    user_options = [
                    ('background=', 'b', "background image to use (should be 500x500px)"),
                    ('dist-dir=', 'd', "directory to put final built distributions in")
                   ]
    description = "Create a Mac OS X disk image with the package installer"

    def initialize_options(self):
        self.dist_dir = None
        self.background = None
        self.finalized = False

    def finalize_options(self):
        self.set_undefined_options('bdist', ('dist_dir', 'dist_dir'))
        self.finalized = True

    def run(self):

        pkg_dir = os.path.join(self.dist_dir, 'pkg')

        # Remove directory if it already exists
        if os.path.exists(pkg_dir):
            shutil.rmtree(pkg_dir)

        # First create the package installer with bdist_mpkg
        mpkg = self.reinitialize_command('bdist_mpkg', reinit_subcommands=1)
        mpkg.dist_dir = pkg_dir
        mpkg.ensure_finalized()
        mpkg.run()

        # Find the name of the pkg file. Since we removed the dist directory
        # at the start of the script, our pkg should be the only file there.
        files = os.listdir(pkg_dir)
        if len(files) != 1:
            raise DistutilsFileError(
                "Expected a single file in the {pkg_dir} "
                "directory".format(pkg_dir=pkg_dir))
        pkg_file = os.path.basename(files[0])
        pkg_name = os.path.splitext(pkg_file)[0]

        # Build the docs
        docs = self.reinitialize_command('build_sphinx', reinit_subcommands=1)
        docs.ensure_finalized()
        docs.run()

        # Copy over the docs to the dist directory
        shutil.copytree(os.path.join(docs.build_dir, 'html'),
                        os.path.join(pkg_dir, 'Documentation'))

        # Copy over the background to the disk image
        if self.background is not None:
            background_dir = os.path.join(pkg_dir, '.background')
            os.mkdir(background_dir)
            shutil.copy2(self.background,
                         os.path.join(background_dir, 'background.png'))

        # Start creating the volume
        dmg_path = os.path.join(self.dist_dir, pkg_name + '.dmg')
        dmg_path_tmp = os.path.join(self.dist_dir, pkg_name + '_tmp.dmg')
        volume_name = pkg_name

        # Remove existing dmg files
        if os.path.exists(dmg_path):
            os.remove(dmg_path)
        if os.path.exists(dmg_path_tmp):
            os.remove(dmg_path_tmp)

        # Check if a volume is already mounted
        volume_path = os.path.join('/', 'Volumes', volume_name)
        if os.path.exists(volume_path):
            raise DistutilsFileError(
                "A volume named {volume_name} is already mounted - please "
                "eject this and try again".format(volume_name=volume_name))

        shell_script = """

        # Create DMG file
        hdiutil create -volname {volume_name} -srcdir {pkg_dir} -fs HFS+ -fsargs "-c c=64,a=16,e=16" -format UDRW -size 24m {dmg_path_tmp}

        # Mount disk image, and keep reference to device
        device=$(hdiutil attach -readwrite -noverify -noautoopen {dmg_path_tmp} | egrep '^/dev/' | sed 1q | awk '{{print $1}}')

        echo '
        tell application "Finder"
            tell disk "{volume_name}"
                open
                set current view of container window to icon view
                set toolbar visible of container window to false
                set statusbar visible of container window to false
                set the bounds of container window to {{100, 100, 600, 600}}
                set theViewOptions to the icon view options of container window
                set arrangement of theViewOptions to not arranged
                set icon size of theViewOptions to 128
                set the background picture of theViewOptions to file ".background:background.png"
                set position of item "{pkg_file}" of container window to {{125, 320}}
                set position of item "Documentation" of container window to {{375, 320}}
                close
                open
                update without registering applications
                delay 5
            end tell
        end tell
        ' | osascript

        # Eject disk image
        hdiutil detach ${{device}}

        # Convert to final read-only disk image
        hdiutil convert {dmg_path_tmp} -format UDZO -imagekey zlib-level=9 -o {dmg_path}

        """.format(volume_name=volume_name, pkg_dir=pkg_dir,
                   pkg_file=pkg_file, dmg_path_tmp=dmg_path_tmp,
                   dmg_path=dmg_path)

        # Make the disk image with the above shell script
        os.system(shell_script)

        # Remove temporary disk image
        os.remove(dmg_path_tmp)
