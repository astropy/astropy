# Licensed under a 3-clause BSD style license - see LICENSE.rst
""" This module contains functions to determine where configuration and
data/cache files used by Astropy should be placed.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from ..extern import six

import os
import sys


__all__ = ['get_config_dir', 'get_cache_dir']


def _find_home():
    """ Locates and return the home directory (or best approximation) on this
    system.

    Raises
    ------
    OSError
        If the home directory cannot be located - usually means you are running
        Astropy on some obscure platform that doesn't have standard home
        directories.
    """


    # this is used below to make fix up encoding issues that sometimes crop up
    # in py2.x but not in py3.x
    if six.PY2:
        decodepath = lambda pth: pth.decode(sys.getfilesystemencoding())
    elif six.PY3:
        decodepath = lambda pth: pth

    # First find the home directory - this is inspired by the scheme ipython
    # uses to identify "home"
    if os.name == 'posix':
        # Linux, Unix, AIX, OS X
        if 'HOME' in os.environ:
            homedir = decodepath(os.environ['HOME'])
        else:
            raise OSError('Could not find unix home directory to search for '
                          'astropy config dir')
    elif os.name == 'nt':  # This is for all modern Windows (NT or after)
        if 'MSYSTEM' in os.environ and os.environ.get('HOME'):
            # Likely using an msys shell; use whatever it is using for its
            # $HOME directory
            homedir = decodepath(os.environ['HOME'])
        # Next try for a network home
        elif 'HOMESHARE' in os.environ:
            homedir = decodepath(os.environ['HOMESHARE'])
        # See if there's a local home
        elif 'HOMEDRIVE' in os.environ and 'HOMEPATH' in os.environ:
            homedir = os.path.join(os.environ['HOMEDRIVE'],
                                   os.environ['HOMEPATH'])
            homedir = decodepath(homedir)
        # Maybe a user profile?
        elif 'USERPROFILE' in os.environ:
            homedir = decodepath(os.path.join(os.environ['USERPROFILE']))
        else:
            try:
                from ..extern.six.moves import winreg as wreg
                shell_folders = r'Software\Microsoft\Windows\CurrentVersion\Explorer\Shell Folders'
                key = wreg.OpenKey(wreg.HKEY_CURRENT_USER, shell_folders)

                homedir = wreg.QueryValueEx(key, 'Personal')[0]
                homedir = decodepath(homedir)
                key.Close()
            except:
                # As a final possible resort, see if HOME is present
                if 'HOME' in os.environ:
                    homedir = decodepath(os.environ['HOME'])
                else:
                    raise OSError('Could not find windows home directory to '
                                  'search for astropy config dir')
    else:
        # for other platforms, try HOME, although it probably isn't there
        if 'HOME' in os.environ:
            homedir = decodepath(os.environ['HOME'])
        else:
            raise OSError('Could not find a home directory to search for '
                          'astropy config dir - are you on an unspported '
                          'platform?')
    return homedir


def get_config_dir(create=True):
    """
    Determines the Astropy configuration directory name and creates the
    directory if it doesn't exist.

    This directory is typically ``$HOME/.astropy/config``, but if the
    XDG_CONFIG_HOME environment variable is set and the
    ``$XDG_CONFIG_HOME/astropy`` directory exists, it will be that directory.
    If neither exists, the former will be created and symlinked to the latter.

    Returns
    -------
    configdir : str
        The absolute path to the configuration directory.

    """

    # symlink will be set to this if the directory is created
    linkto = None
    # first look for XDG_CONFIG_HOME
    xch = os.environ.get('XDG_CONFIG_HOME')

    if xch is not None and os.path.exists(xch):
        xchpth = os.path.join(xch, 'astropy')
        if not os.path.islink(xchpth):
            if os.path.exists(xchpth):
                return os.path.abspath(xchpth)
            else:
                linkto = xchpth
    return os.path.abspath(_find_or_create_astropy_dir('config', linkto))


def get_cache_dir():
    """
    Determines the Astropy cache directory name and creates the directory if it
    doesn't exist.

    This directory is typically ``$HOME/.astropy/cache``, but if the
    XDG_CACHE_HOME environment variable is set and the
    ``$XDG_CACHE_HOME/astropy`` directory exists, it will be that directory.
    If neither exists, the former will be created and symlinked to the latter.

    Returns
    -------
    cachedir : str
        The absolute path to the cache directory.

    """

    # symlink will be set to this if the directory is created
    linkto = None
    # first look for XDG_CACHE_HOME
    xch = os.environ.get('XDG_CACHE_HOME')

    if xch is not None and os.path.exists(xch):
        xchpth = os.path.join(xch, 'astropy')
        if not os.path.islink(xchpth):
            if os.path.exists(xchpth):
                return os.path.abspath(xchpth)
            else:
                linkto = xchpth

    return os.path.abspath(_find_or_create_astropy_dir('cache', linkto))


def _find_or_create_astropy_dir(dirnm, linkto):
    innerdir = os.path.join(_find_home(), '.astropy')
    maindir = os.path.join(_find_home(), '.astropy', dirnm)

    if not os.path.exists(maindir):
        # first create .astropy dir if needed
        if not os.path.exists(innerdir):
            try:
                os.mkdir(innerdir)
            except OSError:
                if not os.path.isdir(innerdir):
                    raise
        elif not os.path.isdir(innerdir):
            msg = 'Intended Astropy directory {0} is actually a file.'
            raise IOError(msg.format(innerdir))

        try:
            os.mkdir(maindir)
        except OSError:
            if not os.path.isdir(maindir):
                raise

        if (not sys.platform.startswith('win') and
            linkto is not None and
                not os.path.exists(linkto)):
            os.symlink(maindir, linkto)

    elif not os.path.isdir(maindir):
        msg = 'Intended Astropy {0} directory {1} is actually a file.'
        raise IOError(msg.format(dirnm, maindir))

    return os.path.abspath(maindir)
