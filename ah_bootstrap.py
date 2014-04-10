import contextlib
import errno
import imp
import os
import re
import subprocess as sp
import sys

try:
    from ConfigParser import ConfigParser
except ImportError:
    from configparser import ConfigParser


if sys.version_info[0] < 3:
    _str_types = (str, unicode)
else:
    _str_types = (str, bytes)

# Some pre-setuptools checks to ensure that either distribute or setuptools >=
# 0.7 is used (over pre-distribute setuptools) if it is available on the path;
# otherwise the latest setuptools will be downloaded and bootstrapped with
# ``ez_setup.py``.  This used to be included in a separate file called
# setuptools_bootstrap.py; but it was combined into ah_bootstrap.py
try:
    import pkg_resources
    _setuptools_req = pkg_resources.Requirement.parse('setuptools>=0.7')
    # This may raise a DistributionNotFound in which case no version of
    # setuptools or distribute is properly instlaled
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

# TODO: Maybe enable checking for a specific version of astropy_helpers?
DIST_NAME = 'astropy-helpers'
PACKAGE_NAME = 'astropy_helpers'

# Defaults for other options
DOWNLOAD_IF_NEEDED = True
INDEX_URL = 'https://pypi.python.org/simple'
USE_GIT = True
AUTO_UPGRADE = True


def use_astropy_helpers(path=None, download_if_needed=None, index_url=None,
                        use_git=None, auto_upgrade=None):
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
        ``setup_requires`` feature of setuptools.

    index_url : str, optional
        If provided, use a different URL for the Python package index than the
        main PyPI server.

    use_git : bool, optional
        If `False` no git commands will be used--this effectively disables
        support for git submodules.

    auto_upgrade : bool, optional
        By default, when installing a package from a non-development source
        distribution ah_boostrap will try to automatically check for patch
        releases to astropy-helpers on PyPI and use the patched version over
        any bundled versions.  Setting this to `False` will disable that
        functionality.
    """

    # True by default, unless the --offline option was provided on the command
    # line
    if '--offline' in sys.argv:
        download_if_needed = False
        auto_upgrade = False
        sys.argv.remove('--offline')

    if path is None:
        path = PACKAGE_NAME

    if download_if_needed is None:
        download_if_needed = DOWNLOAD_IF_NEEDED

    if index_url is None:
        index_url = INDEX_URL

    if use_git is None:
        use_git = USE_GIT

    if auto_upgrade is None:
        auto_upgrade = AUTO_UPGRADE

    # Declared as False by default--later we check if astropy-helpers can be
    # upgraded from PyPI, but only if not using a source distribution (as in
    # the case of import from a git submodule)
    is_submodule = False

    if not isinstance(path, _str_types):
        if path is not None:
            raise TypeError('path must be a string or None')

        if not download_if_needed:
            log.debug('a path was not given and download from PyPI was not '
                      'allowed so this is effectively a no-op')
            return
    elif not os.path.exists(path) or os.path.isdir(path):
        # Even if the given path does not exist on the filesystem, if it *is* a
        # submodule, `git submodule init` will create it
        is_submodule = use_git and _check_submodule(path)

        if is_submodule or os.path.isdir(path):
            log.info(
                'Attempting to import astropy_helpers from {0} {1!r}'.format(
                    'submodule' if is_submodule else 'directory', path))
            dist = _directory_import(path)
        else:
            dist = None

        if dist is None:
            msg = (
                'The requested path {0!r} for importing {1} does not '
                'exist, or does not contain a copy of the {1} pacakge.  '
                'Attempting download instead.'.format(path, PACKAGE_NAME))
            if download_if_needed:
                log.warn(msg)
            else:
                raise _AHBootstrapSystemExit(msg)
    elif os.path.isfile(path):
        # Handle importing from a source archive; this also uses setup_requires
        # but points easy_install directly to the source archive
        try:
            dist = _do_download(find_links=[path])
        except Exception as e:
            if download_if_needed:
                log.warn('{0}\nWill attempt to download astropy_helpers from '
                         'PyPI instead.'.format(str(e)))
                dist = None
            else:
                raise _AHBootstrapSystemExit(e.args[0])
    else:
        msg = ('{0!r} is not a valid file or directory (it could be a '
               'symlink?)'.format(path))
        if download_if_needed:
            log.warn(msg)
            dist = None
        else:
            raise _AHBootstrapSystemExit(msg)

    if dist is not None and auto_upgrade and not is_submodule:
        # A version of astropy-helpers was found on the available path, but
        # check to see if a bugfix release is available on PyPI
        upgrade = _do_upgrade(dist, index_url)
        if upgrade is not None:
            dist = upgrade
    elif dist is None:
        # Last resort--go ahead and try to download the latest version from
        # PyPI
        try:
            if download_if_needed:
                log.warn(
                    "Downloading astropy_helpers; run setup.py with the "
                    "--offline option to force offline installation.")
                dist = _do_download(index_url=index_url)
            else:
                raise _AHBootstrapSystemExit(
                    "No source for the astropy_helpers package; "
                    "astropy_helpers must be available as a prerequisite to "
                    "installing this package.")
        except Exception as e:
            if DEBUG:
                raise
            else:
                raise _AHBootstrapSystemExit(e.args[0])

    if dist is not None:
        # Otherwise we found a version of astropy-helpers so we're done
        # Just activate the found distribibution on sys.path--if we did a
        # download this usually happens automatically but do it again just to
        # be sure
        return dist.activate()


def _do_download(version='', find_links=None, index_url=None):
    try:
        if find_links:
            allow_hosts = ''
            index_url = None
        else:
            allow_hosts = None
        # Annoyingly, setuptools will not handle other arguments to
        # Distribution (such as options) before handling setup_requires, so it
        # is not straightfoward to programmatically augment the arguments which
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
        if DEBUG:
            dist = _Distribution(attrs=attrs)
        else:
            with _silence():
                dist = _Distribution(attrs=attrs)

        return pkg_resources.working_set.by_key.get(DIST_NAME)
    except Exception as e:
        if DEBUG:
            raise

        msg = 'Error retrieving astropy helpers from {0}:\n{1}'
        if find_links:
            source = find_links[0]
        elif index_url:
            source = index_url
        else:
            source = 'PyPI'

        raise Exception(msg.format(source, repr(e)))


def _do_upgrade(dist, index_url):
    # Build up a requirement for a higher bugfix release but a lower minor
    # release (so API compatibility is guaranteed)
    # sketchy version parsing--maybe come up with something a bit more
    # robust for this
    major, minor = (int(part) for part in dist.parsed_version[:2])
    next_minor = '.'.join([str(major), str(minor + 1), '0'])
    req = pkg_resources.Requirement.parse(
        '{0}>{1},<{2}'.format(DIST_NAME, dist.version, next_minor))

    package_index = PackageIndex(index_url=index_url)

    upgrade = package_index.obtain(req)

    if upgrade is not None:
        return _do_download(version=upgrade.version, index_url=index_url)


def _directory_import(path):
    """
    Import astropy_helpers from the given path, which will be added to
    sys.path.

    Must return True if the import succeeded, and False otherwise.
    """

    # Return True on success, False on failure but download is allowed, and
    # otherwise raise SystemExit
    path = os.path.abspath(path)
    pkg_resources.working_set.add_entry(path)
    dist = pkg_resources.working_set.by_key.get(DIST_NAME)

    if dist is None:
        # We didn't find an egg-info/dist-info in the given path, but if a
        # setup.py exists we can generate it
        setup_py = os.path.join(path, 'setup.py')
        if os.path.isfile(setup_py):
            with _silence():
                run_setup(os.path.join(path, 'setup.py'), ['egg_info'])

            for dist in pkg_resources.find_distributions(path, True):
                # There should be only one...
                pkg_resources.working_set.add(dist, path, False)
                break

    return dist


def _check_submodule(path):
    try:
        p = sp.Popen(['git', 'submodule', 'status', '--', path],
                     stdout=sp.PIPE, stderr=sp.PIPE)
        stdout, stderr = p.communicate()
    except OSError as e:
        if DEBUG:
            raise

        if e.errno == errno.ENOENT:
            # The git command simply wasn't found; this is most likely the
            # case on user systems that don't have git and are simply
            # trying to install the package from PyPI or a source
            # distribution.  Silently ignore this case and simply don't try
            # to use submodules
            return False
        else:
            raise _AHBoostrapSystemExit(
                'An unexpected error occurred when running the '
                '`git submodule status` command:\n{0}'.format(str(e)))


    if p.returncode != 0 or stderr:
        # Unfortunately the return code alone cannot be relied on, as
        # earler versions of git returned 0 even if the requested submodule
        # does not exist
        log.debug('git submodule command failed '
                  'unexpectedly:\n{0}'.format(stderr))
        return False
    else:
        # The stdout should only contain one line--the status of the
        # requested submodule
        m = _git_submodule_status_re.match(stdout)
        if m:
            # Yes, the path *is* a git submodule
            _update_submodule(m.group('submodule'), m.group('status'))
            return True
        else:
            log.warn(
                'Unexected output from `git submodule status`:\n{0}\n'
                'Will attempt import from {1!r} regardless.'.format(
                    stdout, path))
            return False


def _update_submodule(submodule, status):
    if status == b' ':
        # The submodule is up to date; no action necessary
        return
    elif status == b'-':
        cmd = ['update', '--init']
        log.info('Initializing submodule {0!r}'.format(submodule))
    elif status == b'+':
        cmd = ['update']
        log.info('Updating submodule {0!r}'.format(submodule))
    elif status == b'U':
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

    try:
        p = sp.Popen(['git', 'submodule'] + cmd + ['--', submodule],
                     stdout=sp.PIPE, stderr=sp.PIPE)
        stdout, stderr = p.communicate()
    except OSError as e:
        err_msg = str(e)
    else:
        if p.returncode != 0:
            err_msg = stderr

    if err_msg:
        log.warn('An unexpected error occurred updating the git submodule '
                 '{0!r}:\n{1}\n{2}'.format(submodule, err_msg, _err_help_msg))


class _DummyFile(object):
    """A noop writeable object."""

    errors = ''  # Required for Python 3.x

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


# Output of `git submodule status` is as follows:
#
# 1: Status indicator: '-' for submodule is uninitialized, '+' if submodule is
# initialized but is not at the commit currently indicated in .gitmodules (and
# thus needs to be updated), or 'U' if the submodule is in an unstable state
# (i.e. has merge conflicts)
#
# 2. SHA-1 hash of the current commit of the submodule (we don't really need
# this information but it's useful for checking that the output is correct)
#
# 3. The output of `git describe` for the submodule's current commit hash (this
# includes for example what branches the commit is on) but only if the
# submodule is initialized.  We ignore this information for now
_git_submodule_status_re = re.compile(
    b'^(?P<status>[+-U ])(?P<commit>[0-9a-f]{40}) (?P<submodule>\S+)( .*)?$')


# Implement the auto-use feature; this allows use_astropy_helpers() to be used
# at import-time automatically so long as the correct options are specified in
# setup.cfg
_CFG_OPTIONS = [('auto_use', bool), ('path', str),
                ('download_if_needed', bool), ('index_ur', str),
                ('use_git', bool), ('auto_upgrade', bool)]

def _main():
    if not os.path.exists('setup.cfg'):
        return

    cfg = ConfigParser()

    try:
        cfg.read('setup.cfg')
    except Exception as e:
        if DEBUG:
            raise

        log.error(
            "Error reading setup.cfg: {0!r}\nastropy_helpers will not be "
            "automatically bootstrapped and package installation may fail."
            "\n{1}".format(e, _err_help_msg))
        return

    if not cfg.has_section('ah_bootstrap'):
        return

    kwargs = {}

    for option, type_ in _CFG_OPTIONS:
        if not cfg.has_option('ah_bootstrap', option):
            continue

        if type_ is bool:
            value = cfg.getboolean('ah_bootstrap', option)
        else:
            value = cfg.get('ah_bootstrap', option)

        kwargs[option] = value

    if kwargs.pop('auto_use', False):
        use_astropy_helpers(**kwargs)


_main()
