#!/usr/bin/env python
from __future__ import division

"""This module contains helper functions for accessing, downloading, and 
caching data files.
"""

__all__ = ['get_data_fileobj','get_data_filename','compute_hash',
           'clear_data_cache']

#TODO: replace this with configobj config setting
DATAURL = 'http://data.astropy.org/'

def get_data_fileobj(dataname,cache=True):
    """
    Retrieves a data file from the standard locations and provides the file as
    a file-like object.
    
    Parameters
    ----------
    dataname : str
        The locator for the data file.  One of the following:
        
            * The name of a data file included in the source distribution.  This
              can be in a directory, in which case the directory is expected to
              be inside the source code 'data' directory.
            * The name of a data file stored on the astropy data server.
            * A hash referencing a particular version of a file on the Astropy
              data server, e.g. 'hash/395dd6493cc584df1e78b474fb150840'
            * A URL to some other file.
            
    cache : bool
        If True, the file will be downloaded and saved locally.  If False, the 
        file-like object will directly access the resource (e.g. if a remote URL
        is accessed, an objectect like that from `urllib2.urlopen` is returned).
    
    Returns
    -------
    fileobj : file-like
        An object with the contents of the data file available via :func:`read`.
        
    """
    from urlparse import urlparse
    from urllib2 import urlopen
    
    if cache:
        return open(get_data_filename(dataname),'r')
    else:
        url = urlparse(dataname)
        if url.scheme!='':
            #it's actually a url for a net location
            return urlopen(datafn)
        else:
            datafn = _find_pkg_data_fn(dataname)
            if datafn is None:
                #not local file - need to get remote data
                return urlopen(DATAURL+datafn)
            else:
                return open(datafn,'r')

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
        The locator for the data file.  One of the following:
        
            * The name of a data file included in the source distribution.  This
              can be in a directory, in which case the directory is expected to
              be inside the source code 'data' directory.
            * The name of a data file stored on the astropy data server.
            * A hash referencing a particular version of a file on the Astropy
              data server, e.g. 'hash/395dd6493cc584df1e78b474fb150840'
            * A URL to some other file.
    
    Returns
    -------
    filename : str
        A file path on the local file system corresponding to the data requested
        in `dataname`.
    """
    from urlparse import urlparse
    
    url = urlparse(dataname)
    if url.scheme!='':
        #it's actually a url for a net location
        return _cache_remote(datafn)
    else:
        datafn = _find_pkg_data_fn(dataname)
        if datafn is None:
            #no local - need to get remote data
            return _cache_remote(DATAURL+datafn)
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
    
def _cache_remote(remoteurl):
    """
    Accepts a URL, downloads and caches the result returning the resulting 
    filename. If present in the cache, just returns the filename.
    """
    import urllib
    import shelve
    from os import mkdir
    from os.path import join,exists
    from urllib2 import urlopen
    from urlparse import urlsplit
    from contextlib import closing
    from .configs import get_config_dir
    
    cfgdir = get_config_dir()
    dldir = join(cfgdir,'datacache')
    if not exists(dldir):
        mkdir(dldir)
    
   
    
    #use a special mapping file to determine if this url is already downloaded
    with closing(shelve.open(join(cfgdir,'datacache_urlmap'))) as url2fn:
        if remoteurl in url2fn:
            localpath =  url2fn[remoteurl]
        else:
            #if not in the mapping file, download the file to the cache
            with closing(urlopen(remoteurl)) as remote:
                #determine the proper local name for the file from the 
                #headers if possible
                #TODO: make use of actual hash file name for hash/... - need data server info for this
                if 'Content-Disposition' in rinfo:
                    #often URLs that redirect to a download provide the fielname
                    #in the header info
                    localfn = rinfo['Content-Disposition'].split('filename=')[1]
                else:
                    #otherwise fallback on the url filename
                    localfn = urlsplit(dataurl)[2].split('/')[-1]
                
                    
                localpath = join(dldir,localfn)
                
                with open(localpath,'w') as f:
                    #TODO: add in download reporter that sends a message to the
                    # log when the log is implemented
                    f.write(remote.read())
            url2fn[remoteurl] = localpath
            
    return localpath

def clear_data_cache(filename=None):
    """ Clears the data file cache by deleting the local version.
    
    Parameters
    ----------
    filename : str or None
        If a specific file or URL, the associated file will be deleted.  
        If None, all cached data files will be removed.
        
    Raises
    ------
    OSEerror
        If the requested filename is not present in the data directory.
    """
    from os import unlink
    from os.path import join,exists
    from shutil import rmtree
    from .configs import get_config_dir
    
    dldir = join(get_config_dir(),'datacache')
    if filename is None:
        if exists(dldir):
            rmtree(dldir)
    else:
        unlink(join(dldir,filename))
    