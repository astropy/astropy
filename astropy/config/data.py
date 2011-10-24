#!/usr/bin/env python
from __future__ import division

"""This module contains helper functions for accessing, downloading, and 
caching data files.
"""

__all__ = []

#TODO: replace this with configobj config setting
DATAURL = 'http://data.astropy.org/'

def get_data_fileobj(dataname,cache=True):
    """
    Retrieves a data file from the standard locations and provides the file as
    a file-like object.
    
    """

def get_data_filename(dataname):
    """
    Retrieves a data file from the standard locations and provides the local 
    name of the file.
    
    This function is similar to `get_data_fileobj` but returns the file *name* 
    instead of a readable file-like object.  This means that this function must
    always cache remote files locally, unlike `get_data_fileobj`.
    
    Parameters
    ----------
    dataname : str
        The name 
       
    
    
    Returns
    -------
    filename : str
        A file path on the local file system corresponding to the data requested
        in `dataname`.
    """
    from urlparse import urlparse
    from .configs import get_config_dir
    
    url = urlparse(dataname)
    if url.scheme!='':
        #it's actually a url for a net location
        raise NotImplementedError
    else:
        datafn = _find_pkg_data_fn(dataname)
        if datafn is None:
            #no local - need to get remote data
            return _cache_remote_fn(DATAURL+datafn)
        else:
            return datafn
    
def _find_pkg_data_fn(dataname):
    """
    Look for data in the source-included data directory and return the filename
    or None if the file is not found.
    """
    
    from .. import __file__ as rootfile
    from os import sep,path
    
    #TBD: should this look through all possible data dirs, or just the root data dir?
    
    path = path.dirname(rootfile)+sep+'data'+sep+dataname
    if path.isdir(path):
        raise IOError("Tried to access a data file that's actually a package data directory")
    elif path.isfile(path):
        return path
    else:
        return None
    
def _cache_remote_fn(remoteurl):
    raise NotImplementedError