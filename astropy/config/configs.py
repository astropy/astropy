# Licensed under a 3-clause BSD style license - see LICENSE.rst  
from __future__ import division

"""This module contains functions to standardize access to configuration files
for astropy and affiliated modules.
"""

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

    #First find the home directory - this is inspired by the scheme ipython
    #uses to identify "home"
    if os.name == 'posix':
        # Linux, Unix, AIX, OS X
        if 'HOME' in env:
            homedir = env['HOME'].decode(sys.getfilesystemencoding())
        else:
            raise OSError('Could not find unix home directory to search for' +\
                          ' astropy config dir')
    elif os.name == 'nt':  # This is for all modern Windows (NT or after)
        #Try for a network home first
        if 'HOMESHARE' in env: 
            homedir = env['HOMESHARE'].decode(sys.getfilesystemencoding())
        #See if there's a local home
        elif 'HOMEDRIVE' in env and 'HOMEPATH' in env: 
            homedir = os.path.join(env['HOMEDRIVE'], env['HOMEPATH'])
            homedir = homedir.decode(sys.getfilesystemencoding())
        #Maybe a user profile?
        elif 'USERPROFILE' in env:
            homedir = os.path.join(env['USERPROFILE']).decode(sys.getfilesystemencoding())
        else:            
            try:
                import _winreg as wreg
                key = wreg.OpenKey(wreg.HKEY_CURRENT_USER,
            'Software\Microsoft\Windows\CurrentVersion\Explorer\Shell Folders')
                
                homedir = wreg.QueryValueEx(key, 'Personal')[0]
                homedir = homedir.decode(sys.getfilesystemencoding())
                key.Close()
            except:
                #As a final possible resort, see if HOME is present for some reason
                if 'HOME' in env:
                    homedir = env['HOME'].decode(sys.getfilesystemencoding())
                else:
                    raise OSError('Could not find windows home directory to' +\
                                  ' search for astropy config dir')
    else:
        #for other platforms, try HOME, although it probably isn't there
        if 'HOME' in env:
            homedir = env['HOME'].decode(sys.getfilesystemencoding())
        else:
            raise OSError('Could not find a home directory to search for ' +\
                 'astropy config dir - are you on an unspported platform?')
    return homedir

def get_config_dir(create=True):
    """
    Determines the Astropy configuration directory name and creates the 
    directory if it doesn't exist.
    
    This directory is typically ``$HOME/.astropy/config``, but if the 
    XDG_CONFIG_HOME environment variable is set and the  
    ``$XDG_CONFIG_HOME/astropy`` directory exists, it will be that directory.
    If it neither exists, the former will be created and symlinked to the 
    latter.
    
    Returns
    -------
    configdir : str
        The absolute path to the configuration directory.
        
    """
    from os import path,environ
    
    #symlink will be set to this if the directory is created
    linkto = None 
    #first look for XDG_CONFIG_HOME
    xch = environ.get('XDG_CONFIG_HOME')

    if xch is not None:
        xchpth =  path.join(xch,'astropy')
        if not path.islink(xchpth):
            if path.exists(xchpth):
                return xchpth
            else:
                linkto = xchpth
                
    return _find_or_create_astropy_dir('config',linkto)
    
def get_cache_dir():
    """
    Determines the Astropy cache directory name and creates the directory if it
    doesn't exist.
    
    This directory is typically ``$HOME/.astropy/cache``, but if the 
    XDG_CACHE_HOME environment variable is set and the  
    ``$XDG_CACHE_HOME/astropy`` directory exists, it will be that directory.
    If it neither exists, the former will be created and symlinked to the 
    latter.
    
    Returns
    -------
    cachedir : str
        The absolute path to the cache directory.
        
    """
    from os import path,environ
    
    #symlink will be set to this if the directory is created
    linkto = None 
    #first look for XDG_CACHE_HOME
    xch = environ.get('XDG_CACHE_HOME')
    
    if xch is not None:
        xchpth =  path.join(xch,'astropy')
        if not path.islink(xchpth):
            if path.exists(xchpth):
                return xchpth
            else:
                linkto = xchpth
                
    return _find_or_create_astropy_dir('cache',linkto)
    
def _find_or_create_astropy_dir(dirnm,linkto):
    from os import path,mkdir,symlink
    
    innerdir = path.join(_find_home(),'.astropy')
    dir = path.join(_find_home(),'.astropy',dirnm)
    
    if not path.exists(dir):
        #first create .astropy dir if needed
        if not path.exists(innerdir):
            mkdir(innerdir)
        elif not path.isdir(innerdir):
            msg = 'Intended Astropy directory {0} is actually a file.'
            raise IOError(msg.format(innerdir))
            
        mkdir(dir)

        if linkto is not None and not path.exists(linkto):
            symlink(dir,linkto)
    elif not path.isdir(dir):
        msg = 'Intended Astropy {0} directory {1} is actually a file.'
        raise IOError(msg.format(dirnm,dir))
    
    return path.abspath(dir)


def get_config(subpackage, mainpackage='astropy'):
    """ Retrieves the object associated with a configuration file for a 
    subpackage.
    
    This function takes care of locating and opening the configuration file for
    an astropy subpackage.  The file itself is stored as a simple text-based
    format like that used in the standard library `ConfigParser`.  However, the
    interface is more fully featured, making use of the `configobj` package, 
    for which documentation can be found at 
    `http://www.voidspace.org.uk/python/configobj.html`_ . See these docs for 
    further information regarding the definition of the `ConfigObj` class or 
    similar.
    
    .. note::
        There is no automatic mechanism for writing `ConfigObj` objects, so if
        you want any configuration changes to be updated in the file, you  
        *must* call the `ConfigObj.write` method after you have made your 
        changes.
    
    
    Parameters
    ----------
    subpackage : str
        The name of the subpackage that the configuration file is desired for.
    
    Returns
    -------
    cfgobj : ConfigObj
        A ConfigObj object for the file storing the requested subpackage's 
        configuration info.
    
    Raises
    ------
    IOError
        If the config directory is a file, or some other problem occurs when 
        accessing the configuration file.
    """
    
    from ..extern.configobj import configobj
    from os import mkdir
    from os.path import join, exists, isfile
    
    cfgdir = join(get_config_dir(True), mainpackage)
    if exists(cfgdir):
        if isfile(cfgdir):
            msg = 'Tried to generate configuration directory {0} but a file' +\
                  ' exists with the same name'
            raise IOError(msg.format(cfgdir))
    else:
        mkdir(cfgdir)
        
    cfgpath = join(cfgdir, subpackage + '.cfg')
    return configobj.ConfigObj(cfgpath)
