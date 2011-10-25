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
    
def compute_hash(localfn):
    """ Computes the MD5 hash for a file.
    
    The hash for a data file is used for looking up data files in a unique 
    fashion. This is of particular use for tests; a test may require a 
    particular version of a particular file, in which case it can be accessed 
    via hash to get the appropriate version. 
    
    Typically, if you wish to write a test that requires a particular data file,
    you will want to submit that file to the astropy data servers, and use
    e.g. ``get_data_filename('hash/a725fa6ba642587436612c2df0451956')``, but 
    with the hash for your file in place of the hash in the example.
    
    Parameters
    ----------
    localfn : str
        The path to the file for which the hash should be generated.
        
    Returns
    -------
    md5hash : str
        The hex digest of the MD5 hash for the contents of the `localfn` file.
    
    """
    import hashlib
    
    with open(localfn) as f:
        h = hashlib.md5(f.read())
    
    return h.hexdigest()
    
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