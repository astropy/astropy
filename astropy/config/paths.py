# Licensed under a 3-clause BSD style license - see LICENSE.rst
""" This module contains functions to determine where configuration and
data/cache files used by Astropy should be placed.
"""

from __future__ import division

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
    import os
    import sys
    from os import environ as env

    # this is used below to make fix up encoding issues that sometimes crop up
    # in py2.x but not in py3.x
    if sys.version_info[0] < 3:  # pragma: py3
        decodepath = lambda pth: pth.decode(sys.getfilesystemencoding())
    else:  # pragma: py2
        decodepath = lambda pth: pth

    #First find the home directory - this is inspired by the scheme ipython
    #uses to identify "home"
    if os.name == 'posix':
        # Linux, Unix, AIX, OS X
        if 'HOME' in env:
            homedir = decodepath(env['HOME'])
        else:
            raise OSError('Could not find unix home directory to search for '
                          'astropy config dir')
    elif os.name == 'nt':  # This is for all modern Windows (NT or after)
        #Try for a network home first
        if 'HOMESHARE' in env:
            homedir = decodepath(env['HOMESHARE'])
        #See if there's a local home
        elif 'HOMEDRIVE' in env and 'HOMEPATH' in env:
            homedir = os.path.join(env['HOMEDRIVE'], env['HOMEPATH'])
            homedir = decodepath(homedir)
        #Maybe a user profile?
        elif 'USERPROFILE' in env:
            homedir = decodepath(os.path.join(env['USERPROFILE']))
        else:
            try:
                import _winreg as wreg
                key = wreg.OpenKey(wreg.HKEY_CURRENT_USER,
            r'Software\Microsoft\Windows\CurrentVersion\Explorer\Shell Folders')

                homedir = wreg.QueryValueEx(key, 'Personal')[0]
                homedir = decodepath(homedir)
                key.Close()
            except:
                #As a final possible resort, see if HOME is present
                if 'HOME' in env:
                    homedir = decodepath(env['HOME'])
                else:
                    raise OSError('Could not find windows home directory to '
                                  'search for astropy config dir')
    else:
        #for other platforms, try HOME, although it probably isn't there
        if 'HOME' in env:
            homedir = decodepath(env['HOME'])
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

    from os import path, environ

    #symlink will be set to this if the directory is created
    linkto = None
    #first look for XDG_CONFIG_HOME
    xch = environ.get('XDG_CONFIG_HOME')

    if xch is not None and path.exists(xch):
        xchpth = path.join(xch, 'astropy')
        if not path.islink(xchpth):
            if path.exists(xchpth):
                return path.abspath(xchpth)
            else:
                linkto = xchpth
    return path.abspath(_find_or_create_astropy_dir('config', linkto))


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
    from os import path, environ

    #symlink will be set to this if the directory is created
    linkto = None
    #first look for XDG_CACHE_HOME
    xch = environ.get('XDG_CACHE_HOME')

    if xch is not None and path.exists(xch):
        xchpth = path.join(xch, 'astropy')
        if not path.islink(xchpth):
            if path.exists(xchpth):
                return path.abspath(xchpth)
            else:
                linkto = xchpth

    return path.abspath(_find_or_create_astropy_dir('cache', linkto))


def _find_or_create_astropy_dir(dirnm, linkto):
    from os import path, mkdir
    import sys

    innerdir = path.join(_find_home(), '.astropy')
    maindir = path.join(_find_home(), '.astropy', dirnm)

    if not path.exists(maindir):
        #first create .astropy dir if needed
        if not path.exists(innerdir):
            try:
                mkdir(innerdir)
            except OSError:
                if not path.isdir(innerdir):
                    raise
        elif not path.isdir(innerdir):
            msg = 'Intended Astropy directory {0} is actually a file.'
            raise IOError(msg.format(innerdir))

        try:
            mkdir(maindir)
        except OSError:
            if not path.isdir(maindir):
                raise

        if (not sys.platform.startswith('win') and
            linkto is not None and
            not path.exists(linkto)):
            from os import symlink
            symlink(maindir, linkto)

    elif not path.isdir(maindir):
        msg = 'Intended Astropy {0} directory {1} is actually a file.'
        raise IOError(msg.format(dirnm, maindir))

    return path.abspath(maindir)
