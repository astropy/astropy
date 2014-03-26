import contextlib
import errno
import imp
import os
import re
import subprocess as sp
import sys


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
from setuptools import Distribution

# TODO: Maybe enable checking for a specific version of astropy_helpers?
DIST_NAME = 'astropy-helpers'
PACKAGE_NAME = 'astropy_helpers'


def use_astropy_helpers(path=None, download_if_needed=True, index_url=None,
                        use_git=True):
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
    """

    if path is None:
        path = PACKAGE_NAME

    if not isinstance(path, _str_types):
        if path is not None:
            raise TypeError('path must be a string or None')

        if not download_if_needed:
            log.debug('a path was not given and download from PyPI was not '
                      'allowed so this is effectively a no-op')
            return
    elif not os.path.exists(path):
        # Even if the given path does not exist on the filesystem, if it *is* a
        # submodule, `git submodule init` will create it
        is_submodule = use_git and _check_submodule(path)
        if is_submodule and _directory_import(path, download_if_needed,
                                              is_submodule=is_submodule):
            # Successfully imported from submodule
            return

        if download_if_needed:
            log.warn('The requested path {0!r} for importing astropy_helpers '
                     'does not exist.  Attempting download '
                     'instead.'.format(path))
        else:
            raise _AHBootstrapSystemExit(
                'Error: The requested path {0!r} for importing '
                'astropy_helpers does not exist.'.format(path))
    elif os.path.isdir(path):
        if _directory_import(path, download_if_needed, is_submodule=use_git):
            return
    elif os.path.isfile(path):
        # Handle importing from a source archive; this also uses setup_requires
        # but points easy_install directly to the source archive
        try:
            _do_download(find_links=[path])
        except Exception as e:
            if download_if_needed:
                log.warn('{0}\nWill attempt to download astropy_helpers from '
                         'PyPI instead.'.format(str(e)))
            else:
                raise _AHBootstrapSystemExit(e.args[0])
    else:
        msg = ('{0!r} is not a valid file or directory (it could be a '
               'symlink?)'.format(path))
        if download_if_needed:
            log.warn(msg)
        else:
            raise _AHBootstrapSystemExit(msg)

    # If we made it this far, go ahead and attempt to download/activate
    try:
        _do_download(index_url=index_url)
    except Exception as e:
        if DEBUG:
            raise
        else:
            raise _AHBootstrapSystemExit(e.args[0])


def _do_download(find_links=None, index_url=None):
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

        attrs = {'setup_requires': [DIST_NAME]}
        if DEBUG:
            dist = _Distribution(attrs=attrs)
        else:
            with _silence():
                dist = _Distribution(attrs=attrs)
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


def _directory_import(path, download_if_needed, is_submodule=None):
    # Return True on success, False on failure but download is allowed, and
    # otherwise raise SystemExit
    # Check to see if the path is a git submodule
    if is_submodule is None:
        is_submodule = _check_submodule(path)

    log.info(
        'Attempting to import astropy_helpers from {0} {1!r}'.format(
            'submodule' if is_submodule else 'directory', path))
    sys.path.insert(0, path)
    try:
        __import__('astropy_helpers')
        return True
    except ImportError:
        sys.path.remove(path)

        if download_if_needed:
            log.warn(
                'Failed to import astropy_helpers from {0!r}; will '
                'attempt to download it from PyPI instead.'.format(path))
        else:
            raise _AHBoostrapSystemExit(
                'Failed to import astropy_helpers from {0!r}:\n'
                '{1}'.format(path))
    # Otherwise, success!


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
