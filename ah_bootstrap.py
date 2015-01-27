"""
This bootstrap module contains code for ensuring that the astropy_helpers
package will be importable by the time the setup.py script runs.  It also
includes some workarounds to ensure that a recent-enough version of setuptools
is being used for the installation.

This module should be the first thing imported in the setup.py of distributions
that make use of the utilities in astropy_helpers.  If the distribution ships
with its own copy of astropy_helpers, this module will first attempt to import
from the shipped copy.  However, it will also check PyPI to see if there are
any bug-fix releases on top of the current version that may be useful to get
past platform-specific bugs that have been fixed.  When running setup.py, use
the ``--offline`` command-line option to disable the auto-upgrade checks.

When this module is imported or otherwise executed it automatically calls a
main function that attempts to read the project's setup.cfg file, which it
checks for a configuration section called ``[ah_bootstrap]`` the presences of
that section, and options therein, determine the next step taken:  If it
contains an option called ``auto_use`` with a value of ``True``, it will
automatically call the main function of this module called
`use_astropy_helpers` (see that function's docstring for full details).
Otherwise no further action is taken (however,
``ah_bootstrap.use_astropy_helpers`` may be called manually from within the
setup.py script).

Additional options in the ``[ah_boostrap]`` section of setup.cfg have the same
names as the arguments to `use_astropy_helpers`, and can be used to configure
the bootstrap script when ``auto_use = True``.

See https://github.com/astropy/astropy-helpers for more details, and for the
latest version of this module.
"""

import contextlib
import errno
import imp
import io
import locale
import os
import re
import subprocess as sp
import sys

try:
    from ConfigParser import ConfigParser, RawConfigParser
except ImportError:
    from configparser import ConfigParser, RawConfigParser


if sys.version_info[0] < 3:
    _str_types = (str, unicode)
    _text_type = unicode
    PY3 = False
else:
    _str_types = (str, bytes)
    _text_type = str
    PY3 = True

# Some pre-setuptools checks to ensure that either distribute or setuptools >=
# 0.7 is used (over pre-distribute setuptools) if it is available on the path;
# otherwise the latest setuptools will be downloaded and bootstrapped with
# ``ez_setup.py``.  This used to be included in a separate file called
# setuptools_bootstrap.py; but it was combined into ah_bootstrap.py
try:
    import pkg_resources
    _setuptools_req = pkg_resources.Requirement.parse('setuptools>=0.7')
    # This may raise a DistributionNotFound in which case no version of
    # setuptools or distribute is properly installed
    _setuptools = pkg_resources.get_distribution('setuptools')
    if _setuptools not in _setuptools_req:
        # Older version of setuptools; check if we have distribute; again if
        # this results in DistributionNotFound we want to give up
        _distribute = pkg_resources.get_distribution('distribute')
        if _setuptools != _distribute:
            # It's possible on some pathological systems to have an old version
            # of setuptools and distribute on sys.path simultaneously; make
            # sure distribute is the one that's used
            sys.path.insert(1, _distribute.location)
            _distribute.activate()
            imp.reload(pkg_resources)
except:
    # There are several types of exceptions that can occur here; if all else
    # fails bootstrap and use the bootstrapped version
    from ez_setup import use_setuptools
    use_setuptools()

from distutils import log
from distutils.debug import DEBUG


# In case it didn't successfully import before the ez_setup checks
import pkg_resources

from setuptools import Distribution
from setuptools.package_index import PackageIndex
from setuptools.sandbox import run_setup

# Note: The following import is required as a workaround to
# https://github.com/astropy/astropy-helpers/issues/89; if we don't import this
# module now, it will get cleaned up after `run_setup` is called, but that will
# later cause the TemporaryDirectory class defined in it to stop working when
# used later on by setuptools
try:
    import setuptools.py31compat
except ImportError:
    pass

# TODO: Maybe enable checking for a specific version of astropy_helpers?
DIST_NAME = 'astropy-helpers'
PACKAGE_NAME = 'astropy_helpers'

# Defaults for other options
DOWNLOAD_IF_NEEDED = True
INDEX_URL = 'https://pypi.python.org/simple'
USE_GIT = True
OFFLINE = False
AUTO_UPGRADE = True

# A list of all the configuration options and their required types
CFG_OPTIONS = [
    ('auto_use', bool), ('path', str), ('download_if_needed', bool),
    ('index_url', str), ('use_git', bool), ('offline', bool),
    ('auto_upgrade', bool)
]


class _Bootstrapper(object):
    """
    Bootstrapper implementation.  See ``use_astropy_helpers`` for parameter
    documentation.
    """

    def __init__(self, path=None, index_url=None, use_git=None, offline=None,
                 download_if_needed=None, auto_upgrade=None):

        if path is None:
            path = PACKAGE_NAME

        if not (isinstance(path, _str_types) or path is False):
            raise TypeError('path must be a string or False')

        if PY3 and not isinstance(path, _text_type):
            fs_encoding = sys.getfilesystemencoding()
            path = path.decode(fs_encoding)  # path to unicode

        self.path = path

        # Set other option attributes, using defaults where necessary
        self.index_url = index_url if index_url is not None else INDEX_URL
        self.offline = offline if offline is not None else OFFLINE

        # If offline=True, override download and auto-upgrade
        if self.offline:
            download_if_needed = False
            auto_upgrade = False

        self.download = (download_if_needed
                         if download_if_needed is not None
                         else DOWNLOAD_IF_NEEDED)
        self.auto_upgrade = (auto_upgrade
                             if auto_upgrade is not None else AUTO_UPGRADE)

        # If this is a release then the .git directory will not exist so we
        # should not use git.
        git_dir_exists = os.path.exists(os.path.join(os.path.dirname(__file__), '.git'))
        if use_git is None and not git_dir_exists:
            use_git = False

        self.use_git = use_git if use_git is not None else USE_GIT
        # Declared as False by default--later we check if astropy-helpers can be
        # upgraded from PyPI, but only if not using a source distribution (as in
        # the case of import from a git submodule)
        self.is_submodule = False

    @classmethod
    def main(cls, argv=None):
        if argv is None:
            argv = sys.argv

        config = cls.parse_config()
        config.update(cls.parse_command_line(argv))

        auto_use = config.pop('auto_use', False)
        bootstrapper = cls(**config)

        if auto_use:
            # Run the bootstrapper, otherwise the setup.py is using the old
            # use_astropy_helpers() interface, in which case it will run the
            # bootstrapper manually after reconfiguring it.
            bootstrapper.run()

        return bootstrapper

    @classmethod
    def parse_config(cls):
        if not os.path.exists('setup.cfg'):
            return {}

        cfg = ConfigParser()

        try:
            cfg.read('setup.cfg')
        except Exception as e:
            if DEBUG:
                raise

            log.error(
                "Error reading setup.cfg: {0!r}\n{1} will not be "
                "automatically bootstrapped and package installation may fail."
                "\n{2}".format(e, PACKAGE_NAME, _err_help_msg))
            return {}

        if not cfg.has_section('ah_bootstrap'):
            return {}

        config = {}

        for option, type_ in CFG_OPTIONS:
            if not cfg.has_option('ah_bootstrap', option):
                continue

            if type_ is bool:
                value = cfg.getboolean('ah_bootstrap', option)
            else:
                value = cfg.get('ah_bootstrap', option)

            config[option] = value

        return config

    @classmethod
    def parse_command_line(cls, argv=None):
        if argv is None:
            argv = sys.argv

        config = {}

        # For now we just pop recognized ah_bootstrap options out of the
        # arg list.  This is imperfect; in the unlikely case that a setup.py
        # custom command or even custom Distribution class defines an argument
        # of the same name then we will break that.  However there's a catch22
        # here that we can't just do full argument parsing right here, because
        # we don't yet know *how* to parse all possible command-line arguments.
        if '--no-git' in argv:
            config['use_git'] = False
            argv.remove('--no-git')

        if '--offline' in argv:
            config['offline'] = True
            argv.remove('--offline')

        return config

    def run(self):
        strategies = ['local_directory', 'local_file', 'index']
        dist = None

        # Check to see if the path is a submodule
        self.is_submodule = self._check_submodule()

        for strategy in strategies:
            method = getattr(self, 'get_{0}_dist'.format(strategy))
            dist = method()
            if dist is not None:
                break
        else:
            raise _AHBootstrapSystemExit(
                "No source found for the {0!r} package; {0} must be "
                "available and importable as a prerequisite to building "
                "or installing this package.".format(PACKAGE_NAME))

        # Otherwise we found a version of astropy-helpers, so we're done
        # Just active the found distribution on sys.path--if we did a
        # download this usually happens automatically but it doesn't hurt to
        # do it again
        # Note: Adding the dist to the global working set also activates it
        # (makes it importable on sys.path) by default
        pkg_resources.working_set.add(dist)

    @property
    def config(self):
        """
        A `dict` containing the options this `_Bootstrapper` was configured
        with.
        """

        return dict((optname, getattr(self, optname))
                    for optname, _ in CFG_OPTIONS if hasattr(self, optname))

    def get_local_directory_dist(self):
        """
        Handle importing a vendored package from a subdirectory of the source
        distribution.
        """

        if not os.path.isdir(self.path):
            return

        log.info('Attempting to import astropy_helpers from {0} {1!r}'.format(
                 'submodule' if self.is_submodule else 'directory',
                 self.path))

        dist = self._directory_import()

        if dist is None:
            log.warn(
                'The requested path {0!r} for importing {1} does not '
                'exist, or does not contain a copy of the {1} '
                'package.'.format(self.path, PACKAGE_NAME))
        elif self.auto_upgrade and not self.is_submodule:
            # A version of astropy-helpers was found on the available path, but
            # check to see if a bugfix release is available on PyPI
            upgrade = self._do_upgrade(dist)
            if upgrade is not None:
                dist = upgrade

        return dist

    def get_local_file_dist(self):
        """
        Handle importing from a source archive; this also uses setup_requires
        but points easy_install directly to the source archive.
        """

        if not os.path.isfile(self.path):
            return

        log.info('Attempting to unpack and import astropy_helpers from '
                 '{0!r}'.format(self.path))

        try:
            dist = self._do_download(find_links=[self.path])
        except Exception as e:
            if DEBUG:
                raise

            log.warn(
                'Failed to import {0} from the specified archive {1!r}: '
                '{2}'.format(PACKAGE_NAME, self.path, str(e)))
            dist = None

        if dist is not None and self.auto_upgrade:
            # A version of astropy-helpers was found on the available path, but
            # check to see if a bugfix release is available on PyPI
            upgrade = self._do_upgrade(dist)
            if upgrade is not None:
                dist = upgrade

        return dist

    def get_index_dist(self):
        if not self.download:
            log.warn('Downloading {0!r} disabled.'.format(DIST_NAME))
            return False

        log.warn(
            "Downloading {0!r}; run setup.py with the --offline option to "
            "force offline installation.".format(DIST_NAME))

        try:
            dist = self._do_download()
        except Exception as e:
            if DEBUG:
                raise
            log.warn(
                'Failed to download and/or install {0!r} from {1!r}:\n'
                '{2}'.format(DIST_NAME, self.index_url, str(e)))
            dist = None

        # No need to run auto-upgrade here since we've already presumably
        # gotten the most up-to-date version from the package index
        return dist

    def _directory_import(self):
        """
        Import astropy_helpers from the given path, which will be added to
        sys.path.

        Must return True if the import succeeded, and False otherwise.
        """

        # Return True on success, False on failure but download is allowed, and
        # otherwise raise SystemExit
        path = os.path.abspath(self.path)

        # Use an empty WorkingSet rather than the man
        # pkg_resources.working_set, since on older versions of setuptools this
        # will invoke a VersionConflict when trying to install an upgrade
        ws = pkg_resources.WorkingSet([])
        ws.add_entry(path)
        dist = ws.by_key.get(DIST_NAME)

        if dist is None:
            # We didn't find an egg-info/dist-info in the given path, but if a
            # setup.py exists we can generate it
            setup_py = os.path.join(path, 'setup.py')
            if os.path.isfile(setup_py):
                with _silence():
                    run_setup(os.path.join(path, 'setup.py'),
                              ['egg_info'])

                for dist in pkg_resources.find_distributions(path, True):
                    # There should be only one...
                    return dist

        return dist

    def _do_download(self, version='', find_links=None):
        if find_links:
            allow_hosts = ''
            index_url = None
        else:
            allow_hosts = None
            index_url = self.index_url

        # Annoyingly, setuptools will not handle other arguments to
        # Distribution (such as options) before handling setup_requires, so it
        # is not straightforward to programmatically augment the arguments which
        # are passed to easy_install
        class _Distribution(Distribution):
            def get_option_dict(self, command_name):
                opts = Distribution.get_option_dict(self, command_name)
                if command_name == 'easy_install':
                    if find_links is not None:
                        opts['find_links'] = ('setup script', find_links)
                    if index_url is not None:
                        opts['index_url'] = ('setup script', index_url)
                    if allow_hosts is not None:
                        opts['allow_hosts'] = ('setup script', allow_hosts)
                return opts

        if version:
            req = '{0}=={1}'.format(DIST_NAME, version)
        else:
            req = DIST_NAME

        attrs = {'setup_requires': [req]}

        try:
            if DEBUG:
                _Distribution(attrs=attrs)
            else:
                with _silence():
                    _Distribution(attrs=attrs)

            # If the setup_requires succeeded it will have added the new dist to
            # the main working_set
            return pkg_resources.working_set.by_key.get(DIST_NAME)
        except Exception as e:
            if DEBUG:
                raise

            msg = 'Error retrieving {0} from {1}:\n{2}'
            if find_links:
                source = find_links[0]
            elif index_url != INDEX_URL:
                source = index_url
            else:
                source = 'PyPI'

            raise Exception(msg.format(DIST_NAME, source, repr(e)))

    def _do_upgrade(self, dist):
        # Build up a requirement for a higher bugfix release but a lower minor
        # release (so API compatibility is guaranteed)
        # sketchy version parsing--maybe come up with something a bit more
        # robust for this
        major, minor = (int(part) for part in dist.parsed_version[:2])
        next_minor = '.'.join([str(major), str(minor + 1), '0'])
        req = pkg_resources.Requirement.parse(
            '{0}>{1},<{2}'.format(DIST_NAME, dist.version, next_minor))

        package_index = PackageIndex(index_url=self.index_url)

        upgrade = package_index.obtain(req)

        if upgrade is not None:
            return self._do_download(version=upgrade.version)

    def _check_submodule(self):
        """
        Check if the given path is a git submodule.

        See the docstrings for ``_check_submodule_using_git`` and
        ``_check_submodule_no_git`` for further details.
        """

        if (self.path is None or
                (os.path.exists(self.path) and not os.path.isdir(self.path))):
            return False

        if self.use_git:
            return self._check_submodule_using_git()
        else:
            return self._check_submodule_no_git()

    def _check_submodule_using_git(self):
        """
        Check if the given path is a git submodule.  If so, attempt to initialize
        and/or update the submodule if needed.

        This function makes calls to the ``git`` command in subprocesses.  The
        ``_check_submodule_no_git`` option uses pure Python to check if the given
        path looks like a git submodule, but it cannot perform updates.
        """

        cmd = ['git', 'submodule', 'status', '--', self.path]

        try:
            log.info('Running `{0}`; use the --no-git option to disable git '
                     'commands'.format(' '.join(cmd)))
            returncode, stdout, stderr = run_cmd(cmd)
        except _CommandNotFound:
            # The git command simply wasn't found; this is most likely the
            # case on user systems that don't have git and are simply
            # trying to install the package from PyPI or a source
            # distribution.  Silently ignore this case and simply don't try
            # to use submodules
            return False

        stderr = stderr.strip()

        if returncode != 0 and stderr:
            # Unfortunately the return code alone cannot be relied on, as
            # earlier versions of git returned 0 even if the requested submodule
            # does not exist

            # This is a warning that occurs in perl (from running git submodule)
            # which only occurs with a malformatted locale setting which can
            # happen sometimes on OSX.  See again
            # https://github.com/astropy/astropy/issues/2749
            perl_warning = ('perl: warning: Falling back to the standard locale '
                            '("C").')
            if not stderr.strip().endswith(perl_warning):
                # Some other uknown error condition occurred
                log.warn('git submodule command failed '
                         'unexpectedly:\n{0}'.format(stderr))
                return False

        # Output of `git submodule status` is as follows:
        #
        # 1: Status indicator: '-' for submodule is uninitialized, '+' if
        # submodule is initialized but is not at the commit currently indicated
        # in .gitmodules (and thus needs to be updated), or 'U' if the
        # submodule is in an unstable state (i.e. has merge conflicts)
        #
        # 2. SHA-1 hash of the current commit of the submodule (we don't really
        # need this information but it's useful for checking that the output is
        # correct)
        #
        # 3. The output of `git describe` for the submodule's current commit
        # hash (this includes for example what branches the commit is on) but
        # only if the submodule is initialized.  We ignore this information for
        # now
        _git_submodule_status_re = re.compile(
            '^(?P<status>[+-U ])(?P<commit>[0-9a-f]{40}) '
            '(?P<submodule>\S+)( .*)?$')

        # The stdout should only contain one line--the status of the
        # requested submodule
        m = _git_submodule_status_re.match(stdout)
        if m:
            # Yes, the path *is* a git submodule
            self._update_submodule(m.group('submodule'), m.group('status'))
            return True
        else:
            log.warn(
                'Unexpected output from `git submodule status`:\n{0}\n'
                'Will attempt import from {1!r} regardless.'.format(
                    stdout, self.path))
            return False

    def _check_submodule_no_git(self):
        """
        Like ``_check_submodule_using_git``, but simply parses the .gitmodules file
        to determine if the supplied path is a git submodule, and does not exec any
        subprocesses.

        This can only determine if a path is a submodule--it does not perform
        updates, etc.  This function may need to be updated if the format of the
        .gitmodules file is changed between git versions.
        """

        gitmodules_path = os.path.abspath('.gitmodules')

        if not os.path.isfile(gitmodules_path):
            return False

        # This is a minimal reader for gitconfig-style files.  It handles a few of
        # the quirks that make gitconfig files incompatible with ConfigParser-style
        # files, but does not support the full gitconfig syntax (just enough
        # needed to read a .gitmodules file).
        gitmodules_fileobj = io.StringIO()

        # Must use io.open for cross-Python-compatible behavior wrt unicode
        with io.open(gitmodules_path) as f:
            for line in f:
                # gitconfig files are more flexible with leading whitespace; just
                # go ahead and remove it
                line = line.lstrip()

                # comments can start with either # or ;
                if line and line[0] in (':', ';'):
                    continue

                gitmodules_fileobj.write(line)

        gitmodules_fileobj.seek(0)

        cfg = RawConfigParser()

        try:
            cfg.readfp(gitmodules_fileobj)
        except Exception as exc:
            log.warn('Malformatted .gitmodules file: {0}\n'
                     '{1} cannot be assumed to be a git submodule.'.format(
                         exc, self.path))
            return False

        for section in cfg.sections():
            if not cfg.has_option(section, 'path'):
                continue

            submodule_path = cfg.get(section, 'path').rstrip(os.sep)

            if submodule_path == self.path.rstrip(os.sep):
                return True

        return False

    def _update_submodule(self, submodule, status):
        if status == ' ':
            # The submodule is up to date; no action necessary
            return
        elif status == '-':
            if self.offline:
                raise _AHBootstrapSystemExit(
                    "Cannot initialize the {0} submodule in --offline mode; "
                    "this requires being able to clone the submodule from an "
                    "online repository.".format(submodule))
            cmd = ['update', '--init']
            action = 'Initializing'
        elif status == '+':
            cmd = ['update']
            action = 'Updating'
            if self.offline:
                cmd.append('--no-fetch')
        elif status == 'U':
            raise _AHBoostrapSystemExit(
                'Error: Submodule {0} contains unresolved merge conflicts.  '
                'Please complete or abandon any changes in the submodule so that '
                'it is in a usable state, then try again.'.format(submodule))
        else:
            log.warn('Unknown status {0!r} for git submodule {1!r}.  Will '
                     'attempt to use the submodule as-is, but try to ensure '
                     'that the submodule is in a clean state and contains no '
                     'conflicts or errors.\n{2}'.format(status, submodule,
                                                        _err_help_msg))
            return

        err_msg = None
        cmd = ['git', 'submodule'] + cmd + ['--', submodule]
        log.warn('{0} {1} submodule with: `{2}`'.format(
            action, submodule, ' '.join(cmd)))

        try:
            log.info('Running `{0}`; use the --no-git option to disable git '
                     'commands'.format(' '.join(cmd)))
            returncode, stdout, stderr = run_cmd(cmd)
        except OSError as e:
            err_msg = str(e)
        else:
            if returncode != 0:
                err_msg = stderr

        if err_msg is not None:
            log.warn('An unexpected error occurred updating the git submodule '
                     '{0!r}:\n{1}\n{2}'.format(submodule, err_msg,
                                               _err_help_msg))

class _CommandNotFound(OSError):
    """
    An exception raised when a command run with run_cmd is not found on the
    system.
    """


def run_cmd(cmd):
    """
    Run a command in a subprocess, given as a list of command-line
    arguments.

    Returns a ``(returncode, stdout, stderr)`` tuple.
    """

    try:
        p = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE)
        # XXX: May block if either stdout or stderr fill their buffers;
        # however for the commands this is currently used for that is
        # unlikely (they should have very brief output)
        stdout, stderr = p.communicate()
    except OSError as e:
        if DEBUG:
            raise

        if e.errno == errno.ENOENT:
            msg = 'Command not found: `{0}`'.format(' '.join(cmd))
            raise _CommandNotFound(msg, cmd)
        else:
            raise _AHBoostrapSystemExit(
                'An unexpected error occurred when running the '
                '`{0}` command:\n{1}'.format(' '.join(cmd), str(e)))


    # Can fail of the default locale is not configured properly.  See
    # https://github.com/astropy/astropy/issues/2749.  For the purposes under
    # consideration 'latin1' is an acceptable fallback.
    try:
        stdio_encoding = locale.getdefaultlocale()[1] or 'latin1'
    except ValueError:
        # Due to an OSX oddity locale.getdefaultlocale() can also crash
        # depending on the user's locale/language settings.  See:
        # http://bugs.python.org/issue18378
        stdio_encoding = 'latin1'

    # Unlikely to fail at this point but even then let's be flexible
    if not isinstance(stdout, _text_type):
        stdout = stdout.decode(stdio_encoding, 'replace')
    if not isinstance(stderr, _text_type):
        stderr = stderr.decode(stdio_encoding, 'replace')

    return (p.returncode, stdout, stderr)


class _DummyFile(object):
    """A noop writeable object."""

    errors = ''  # Required for Python 3.x
    encoding = 'utf-8'

    def write(self, s):
        pass

    def flush(self):
        pass


@contextlib.contextmanager
def _silence():
    """A context manager that silences sys.stdout and sys.stderr."""

    old_stdout = sys.stdout
    old_stderr = sys.stderr
    sys.stdout = _DummyFile()
    sys.stderr = _DummyFile()
    exception_occurred = False
    try:
        yield
    except:
        exception_occurred = True
        # Go ahead and clean up so that exception handling can work normally
        sys.stdout = old_stdout
        sys.stderr = old_stderr
        raise

    if not exception_occurred:
        sys.stdout = old_stdout
        sys.stderr = old_stderr


_err_help_msg = """
If the problem persists consider installing astropy_helpers manually using pip
(`pip install astropy_helpers`) or by manually downloading the source archive,
extracting it, and installing by running `python setup.py install` from the
root of the extracted source code.
"""


class _AHBootstrapSystemExit(SystemExit):
    def __init__(self, *args):
        if not args:
            msg = 'An unknown problem occurred bootstrapping astropy_helpers.'
        else:
            msg = args[0]

        msg += '\n' + _err_help_msg

        super(_AHBootstrapSystemExit, self).__init__(msg, *args[1:])


if sys.version_info[:2] < (2, 7):
    # In Python 2.6 the distutils log does not log warnings, errors, etc. to
    # stderr so we have to wrap it to ensure consistency at least in this
    # module
    import distutils

    class log(object):
        def __getattr__(self, attr):
            return getattr(distutils.log, attr)

        def warn(self, msg, *args):
            self._log_to_stderr(distutils.log.WARN, msg, *args)

        def error(self, msg):
            self._log_to_stderr(distutils.log.ERROR, msg, *args)

        def fatal(self, msg):
            self._log_to_stderr(distutils.log.FATAL, msg, *args)

        def log(self, level, msg, *args):
            if level in (distutils.log.WARN, distutils.log.ERROR,
                         distutils.log.FATAL):
                self._log_to_stderr(level, msg, *args)
            else:
                distutils.log.log(level, msg, *args)

        def _log_to_stderr(self, level, msg, *args):
            # This is the only truly 'public' way to get the current threshold
            # of the log
            current_threshold = distutils.log.set_threshold(distutils.log.WARN)
            distutils.log.set_threshold(current_threshold)
            if level >= current_threshold:
                if args:
                    msg = msg % args
                sys.stderr.write('%s\n' % msg)
                sys.stderr.flush()

    log = log()


BOOTSTRAPPER = _Bootstrapper.main()


def use_astropy_helpers(**kwargs):
    """
    Ensure that the `astropy_helpers` module is available and is importable.
    This supports automatic submodule initialization if astropy_helpers is
    included in a project as a git submodule, or will download it from PyPI if
    necessary.

    Parameters
    ----------

    path : str or None, optional
        A filesystem path relative to the root of the project's source code
        that should be added to `sys.path` so that `astropy_helpers` can be
        imported from that path.

        If the path is a git submodule it will automatically be initialzed
        and/or updated.

        The path may also be to a ``.tar.gz`` archive of the astropy_helpers
        source distribution.  In this case the archive is automatically
        unpacked and made temporarily available on `sys.path` as a ``.egg``
        archive.

        If `None` skip straight to downloading.

    download_if_needed : bool, optional
        If the provided filesystem path is not found an attempt will be made to
        download astropy_helpers from PyPI.  It will then be made temporarily
        available on `sys.path` as a ``.egg`` archive (using the
        ``setup_requires`` feature of setuptools.  If the ``--offline`` option
        is given at the command line the value of this argument is overridden
        to `False`.

    index_url : str, optional
        If provided, use a different URL for the Python package index than the
        main PyPI server.

    use_git : bool, optional
        If `False` no git commands will be used--this effectively disables
        support for git submodules. If the ``--no-git`` option is given at the
        command line the value of this argument is overridden to `False`.

    auto_upgrade : bool, optional
        By default, when installing a package from a non-development source
        distribution ah_boostrap will try to automatically check for patch
        releases to astropy-helpers on PyPI and use the patched version over
        any bundled versions.  Setting this to `False` will disable that
        functionality. If the ``--offline`` option is given at the command line
        the value of this argument is overridden to `False`.

    offline : bool, optional
        If `False` disable all actions that require an internet connection,
        including downloading packages from the package index and fetching
        updates to any git submodule.  Defaults to `True`.
    """

    global BOOTSTRAPPER

    config = BOOTSTRAPPER.config
    config.update(**kwargs)

    # Create a new bootstrapper with the updated configuration and run it
    BOOTSTRAPPER = _Bootstrapper(**config)
    BOOTSTRAPPER.run()
