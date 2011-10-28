# Licensed under a 3-clause BSD style license - see LICENSE.rst  
from __future__ import division

"""This module contains functions to standardize access to configuration files
for astropy and affiliated modules.
"""

__all__ = ['get_config_dir', 'get_config']


def get_config_dir(create=True):
    """
    Determines the Astropy configuration directory name.
    
    If a $HOME environment variable is defined, this is used as the base. If
    not, this function goes through various convoluted methods to try to
    identify the user's home directory (inspired by ipython's scheme for
    finding these).
    
    Parameters
    ----------
    create : bool
        If True, the directory will be created if it doesn't exist.
    
    Returns
    -------
    configdir : str
        The absolute path to the configuration directory.
        
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
                          ' astropysics config dir')
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
                                  ' search for astropysics config dir')
    else:
        #for other platforms, try HOME, although it probably isn't there
        if 'HOME' in env:
            homedir = env['HOME'].decode(sys.getfilesystemencoding())
        else:
            raise OSError('Could not find a home directory to search for ' +\
                 'astropysics config dir - are you on an unspported platform?')
        
    configdir = os.path.realpath(os.path.join(homedir, '.astropy'))
        
    if create and not os.path.isdir(configdir):
        #try to create it
        os.mkdir(configdir)
    return configdir


def get_config(subpackage, mainpackage='astropy'):
    """ Retrieves the object associated with a configuration file for a 
    subpackage.
    
    This function takes care of locating and opening the configuration file for
    an astropy subclass.  The file itself is stored as a simple text-based
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
