# Licensed under a 3-clause BSD style license - see LICENSE.rst

##############################################################################
# Note: this file exists only for backward-compatibility purposes - the      #
#       contents have been moved to the separate astropy-helpers package,    #
#       located at https://github.com/astropy/astropy-helpers. Any new       #
#       development or bug fixes should be done there.                       #
##############################################################################

"""
Utilities for generating the version string for Astropy (or an affiliated
package) and the version.py module, which contains version info for the
package.

Within the generated astropy.version module, the `major`, `minor`, and `bugfix`
variables hold the respective parts of the version number (bugfix is '0' if
absent). The `release` variable is True if this is a release, and False if this
is a development version of astropy. For the actual version string, use::

    from astropy.version import version

or::

    from astropy import __version__

"""

from __future__ import division

import datetime
import imp
import os
import subprocess
import sys

from distutils import log
from warnings import warn


def _version_split(version):
    """
    Split a version string into major, minor, and bugfix numbers (with bugfix
    optional, defaulting to 0).
    """

    for prerel in ('.dev', 'a', 'b', 'rc'):
        if prerel in version:
            version = version.split(prerel)[0]

    versplit = version.split('.')
    major = int(versplit[0])
    minor = int(versplit[1])
    bugfix = 0 if len(versplit) < 3 else int(versplit[2])
    return major, minor, bugfix


def update_git_devstr(version, path=None):
    """
    Updates the git revision string if and only if the path is being imported
    directly from a git working copy.  This ensures that the revision number in
    the version string is accurate.
    """

    try:
        # Quick way to determine if we're in git or not - returns '' if not
        devstr = get_git_devstr(sha=True, show_warning=False, path=path)
    except OSError:
        return version

    if not devstr:
        # Probably not in git so just pass silently
        return version

    if 'dev' in version:  # update to the current git revision
        version_base = version.split('.dev', 1)[0]
        devstr = get_git_devstr(sha=False, show_warning=False, path=path)

        return version_base + '.dev' + devstr
    else:
        #otherwise it's already the true/release version
        return version


def get_git_devstr(sha=False, show_warning=True, path=None):
    """
    Determines the number of revisions in this repository.

    Parameters
    ----------
    sha : bool
        If True, the full SHA1 hash will be returned. Otherwise, the total
        count of commits in the repository will be used as a "revision
        number".

    show_warning : bool
        If True, issue a warning if git returns an error code, otherwise errors
        pass silently.

    path : str or None
        If a string, specifies the directory to look in to find the git
        repository.  If None, the location of the file this function is in
        is used to infer the git repository location.  If given a filename it
        uses the directory containing that file.

    Returns
    -------
    devversion : str
        Either a string with the revsion number (if `sha` is False), the
        SHA1 hash of the current commit (if `sha` is True), or an empty string
        if git version info could not be identified.

    """

    from .utils import find_current_module

    if path is None:
        try:
            mod = find_current_module(1, finddiff=True)
            path = os.path.abspath(mod.__file__)
        except (ValueError, AttributeError):
            path = __file__
    if not os.path.isdir(path):
        path = os.path.abspath(os.path.split(path)[0])

    if sha:
        cmd = 'rev-parse'  # Faster for getting just the hash of HEAD
    else:
        cmd = 'rev-list'

    try:
        p = subprocess.Popen(['git', cmd, 'HEAD'], cwd=path,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                             stdin=subprocess.PIPE)
        stdout, stderr = p.communicate()
    except OSError as e:
        if show_warning:
            warn('Error running git: ' + str(e))
        return ''

    if p.returncode == 128:
        if show_warning:
            warn('No git repository present! Using default dev version.')
        return ''
    elif p.returncode != 0:
        if show_warning:
            warn('Git failed while determining revision count: ' + stderr)
        return ''

    if sha:
        return stdout.decode('utf-8')[:40]
    else:
        nrev = stdout.decode('utf-8').count('\n')
        return  str(nrev)


# This is used by setup.py to create a new version.py - see that file for
# details. Note that the imports have to be absolute, since this is also used
# by affiliated packages.

_FROZEN_VERSION_PY_TEMPLATE = """
# Autogenerated by {packagename}'s setup.py on {timestamp}

from astropy.version_helpers import update_git_devstr, get_git_devstr

_last_generated_version = {verstr!r}

version = update_git_devstr(_last_generated_version)
githash = get_git_devstr(sha=True, show_warning=False)

major = {major}
minor = {minor}
bugfix = {bugfix}

release = {rel}
debug = {debug}

try:
    from .utils._compiler import compiler
except ImportError:
    compiler = "unknown"

try:
    from .cython_version import cython_version
except ImportError:
    cython_version = "unknown"
"""[1:]


def _get_version_py_str(packagename, version, release, debug):
    timestamp = str(datetime.datetime.now())
    major, minor, bugfix = _version_split(version)
    if packagename.lower() == 'astropy':
        packagename = 'Astropy'
    else:
        packagename = 'Astropy-affiliated package ' + packagename
    return _FROZEN_VERSION_PY_TEMPLATE.format(packagename=packagename,
                                              timestamp=timestamp,
                                              verstr=version,
                                              major=major,
                                              minor=minor,
                                              bugfix=bugfix,
                                              rel=release, debug=debug)


def generate_version_py(packagename, version, release=None, debug=None):
    """Regenerate the version.py module if necessary."""

    from .setup_helpers import is_distutils_display_option
    from .utils.compat.misc import invalidate_caches

    try:
        version_module = __import__(packagename + '.version',
                                    fromlist=['_last_generated_version',
                                              'version', 'release', 'debug'])
        try:
            last_generated_version = version_module._last_generated_version
        except AttributeError:
            # Older version.py with no _last_generated_version; this will
            # ensure a new version.py is written
            last_generated_version = None
        current_release = version_module.release
        current_debug = version_module.debug
    except ImportError:
        version_module = None
        last_generated_version = None
        current_release = None
        current_debug = None

    if release is None:
        # Keep whatever the current value is, if it exists
        release = bool(current_release)

    if debug is None:
        # Likewise, keep whatever the current value is, if it exists
        debug = bool(current_debug)

    version_py = os.path.join(packagename, 'version.py')

    if (last_generated_version != version or current_release != release or
        current_debug != debug):
        if '-q' not in sys.argv and '--quiet' not in sys.argv:
            log.set_threshold(log.INFO)

        if is_distutils_display_option():
            # Always silence unnecessary log messages when display options are
            # being used
            log.set_threshold(log.WARN)

        log.info('Freezing version number to {0}'.format(version_py))

        with open(version_py, 'w') as f:
            # This overwrites the actual version.py
            f.write(_get_version_py_str(packagename, version, release, debug))

        invalidate_caches()

        if version_module:
            imp.reload(version_module)
