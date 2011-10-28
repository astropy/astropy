# Licensed under a 3-clause BSD style license - see LICENSE.rst  
from __future__ import division

""" This module contains helper functions for accessing, downloading, and
caching data files.
"""

__all__ = ['get_data_fileobj', 'get_data_filename', 'compute_hash',
           'clear_data_cache']

#TODO: replace this with configobj config setting
DATAURL = 'http://data.astropy.org/'
REMOTE_TIMEOUT = 3.  

def get_data_fileobj(dataname, cache=True):
    """
    Retrieves a data file from the standard locations and provides the file as
    a file-like object.

    Parameters
    ----------
    dataname : str
        The locator for the data file.  One of the following:

            * The name of a data file included in the source distribution.
              This can be in a directory, in which case the directory is
              expected to be inside the source code 'data' directory.
            * The name of a data file stored on the astropy data server.
            * A hash referencing a particular version of a file on the Astropy
              data server, e.g. 'hash/395dd6493cc584df1e78b474fb150840'
            * A URL to some other file.

    cache : bool or str
        If True, the file will be downloaded and saved locally.  If False, the
        file-like object will directly access the resource (e.g. if a remote
        URL is accessed, an object like that from `urllib2.urlopen` is
        returned). If a string, that string will be used as the local file name
        in the cache directory.
    
    Returns
    -------
    fileobj : file-like
        An object with the contents of the data file available via
        :func:`read`.
    
    Raises
    ------
    ValueError
        If `cache` is specified, but the name does not match an already
        present local cached version of the file.
    urllib2.URLError
        If a remote file cannot be found.
    IOError
        If problems occur writing or reading a local file.
    
    """
    from urlparse import urlparse
    from urllib2 import urlopen
    
    if cache:
        return open(get_data_filename(dataname, cache), 'r')
    else:
        url = urlparse(dataname)
        if url.scheme != '':
            #it's actually a url for a net location
            return urlopen(dataname,timeout=REMOTE_TIMEOUT)
        else:
            datafn = _find_pkg_data_fn(dataname)
            if datafn is None:
                #not local file - need to get remote data
                return urlopen(DATAURL + datafn,timeout=REMOTE_TIMEOUT)
            else:
                return open(datafn, 'r')


def get_data_filename(dataname, cachename=None):
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
        
            * The name of a data file included in the source distribution.
              This can be in a directory, in which case the directory is 
              expected to be inside the source code 'data' directory.
            * The name of a data file stored on the astropy data server.
            * A hash referencing a particular version of a file on the Astropy
              data server, e.g. 'hash/395dd6493cc584df1e78b474fb150840'
            * A URL to some other file.
    cachename : str of None
        Specifies the local name to be used for the cached file.  If None, the
        name will be automatically determined from the requested file's name.
        
    Raises
    ------
    ValueError
        If `cachename` is specified, but the name does not match an already 
        present local cached version of the file.
    
    Returns
    -------
    filename : str
        A file path on the local file system corresponding to the data 
        requested in `dataname`.
    """
    from urlparse import urlparse
    
    url = urlparse(dataname)
    if url.scheme != '':
        #it's actually a url for a net location
        return _cache_remote(dataname, cachename)
    else:
        datafn = _find_pkg_data_fn(dataname)
        if datafn is None:
            #no local - need to get remote data
            return _cache_remote(DATAURL + datafn, cachename)
        else:
            return datafn

_COMPUTE_HASH_BLOCK_SIZE = 1024 #TBD: replace with a configuration option?

def compute_hash(localfn):
    """ Computes the MD5 hash for a file.
    
    The hash for a data file is used for looking up data files in a unique 
    fashion. This is of particular use for tests; a test may require a 
    particular version of a particular file, in which case it can be accessed 
    via hash to get the appropriate version. 
    
    Typically, if you wish to write a test that requires a particular data 
    file, you will want to submit that file to the astropy data servers, and 
    use e.g. ``get_data_filename('hash/a725fa6ba642587436612c2df0451956')``, 
    but with the hash for your file in place of the hash in the example.
    
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
        h = hashlib.md5()
        block = f.read(_COMPUTE_HASH_BLOCK_SIZE)
        while block!='':
            h.update(block)
            block = f.read(_COMPUTE_HASH_BLOCK_SIZE)
    
    return h.hexdigest()


def _find_pkg_data_fn(dataname):
    """
    Look for data in the source-included data directory and return the filename
    or None if the file is not found.
    """
    
    from .. import __file__ as rootfile
    from os.path import dirname, join
    
    #TBD: should this look through all possible data dirs, or just the 
    
    path = join(dirname(rootfile),'data',dataname)
    if path.isdir(path):
        raise IOError("Tried to access a data file that's actually a package data directory")
    elif path.isfile(path):
        return path
    else:
        return None


def _cache_remote(remoteurl, localname=None):
    """
    Accepts a URL, downloads and caches the result returning the resulting 
    filename. If present in the cache, just returns the filename.
    
    If `localname` is given, it specifies the name to be used for the cached
    version of the file.
    """
    import urllib
    import shelve
    from os import mkdir
    from os.path import join, split, exists
    from urllib2 import urlopen
    from urlparse import urlsplit
    from contextlib import closing
    
    dldir = get_data_cache_dir()
    if not exists(dldir):
        mkdir(dldir)   
    
    #use a special mapping file to determine if this url is already downloaded
    urlmapfn = join(split(dldir)[0], 'datacache_urlmap')
    with closing(shelve.open(urlmapfn)) as url2fn:
        if remoteurl in url2fn:
            localpath = url2fn[remoteurl]
            if localname is not None and join(dldir, localname) != localpath:
                msgstr = 'Requested localname {0} does not match cached data ' +\
                         'file name {1}'
                raise ValueError(msgstr.format(localname, localpath))
        else:
            #if not in the mapping file, download the file to the cache
            with closing(urlopen(remoteurl,timeout=REMOTE_TIMEOUT)) as remote:
                if localname is None:
                    #determine the proper local name for the file from the 
                    #headers if possible
                    #TODO: make use of actual hash file name for hash/... - need data server info for this
                    rinfo = remote.info()
                    if 'Content-Disposition' in rinfo:
                        #often URLs that redirect to a download provide the 
                        #fielname in the header info
                        localname = rinfo['Content-Disposition'].split('filename=')[1]
                    else:
                        #otherwise fallback on the url filename
                        localname = urlsplit(remoteurl)[2].split('/')[-1]
                    
                localpath = join(dldir, localname)
                
                with open(localpath, 'w') as f:
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
    import shelve
    from os import unlink
    from os.path import join, split, exists
    from shutil import rmtree
    from contextlib import closing
    
    dldir = get_data_cache_dir()
    urlmapfn = join(split(dldir)[0], 'datacache_urlmap')
    if filename is None:
        if exists(dldir):
            rmtree(dldir)
        if exists(urlmapfn):
            unlink(urlmapfn)
    else:
        filepath = join(dldir, filename)
        urlmapfn = join(split(dldir)[0], 'datacache_urlmap')
        with closing(shelve.open(urlmapfn)) as url2fn:
            for k,v in url2fn.items():
                if v==filepath:
                    del url2fn[k]
        unlink(filepath)


def get_data_cache_dir():
    """ Finds the path to the data cache directory.
    
    Returns
    -------
    datadir : str
        The path to the data cache directory.
    """
    from .configs import get_config_dir
    from os.path import join
    
    return join(get_config_dir(), 'datacache')
