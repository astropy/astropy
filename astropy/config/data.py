# Licensed under a 3-clause BSD style license - see LICENSE.rst  
""" This module contains helper functions for accessing, downloading, and
caching data files.
"""

from __future__ import division
from .configs import ConfigurationItem,save_config

__all__ = ['get_data_fileobj', 'get_data_filename', 'compute_hash',
           'clear_data_cache']

DATAURL = ConfigurationItem('dataurl','http://data.astropy.org/',
                            'URL for astropy remote data site.')
REMOTE_TIMEOUT = ConfigurationItem('remote_timeout',3.,
                  'Time to wait for remote data query (in seconds).')  
COMPUTE_HASH_BLOCK_SIZE = ConfigurationItem('hash_block_size',2**16, #64K
                            'Block size for computing MD5 file hashes.')  
DATA_CACHE_DL_BLOCK_SIZE = ConfigurationItem('data_dl_block_size',2**16, #64K
                        'Number of bytes of remote data to download per step.')

DATA_CACHE_LOCK_ATTEMPTS = ConfigurationItem('data_cache_lock_attempts',3,
                            'Number of times to try to get the lock while ' + 
                            'accessing the data cache before giving up.')

#used for supporting with statements in get_data_fileobj
def _fake_enter(self):
    return self
def _fake_exit(self, type, value, traceback):
    self.close()

def get_data_fileobj(dataname, cache=True):
    """
    Retrieves a data file from the standard locations and provides the file as
    a file-like object.
    
    .. note::
        If a file is requested using the form used to search inside the source
        distribution (the first form given in the `dataname` parameter 
        description), the subpackage that stores the data will be imported.

    Parameters
    ----------
    dataname : str
        Name/location of the desired data file.  One of the following:

            * The name of a data file included in the source distribution.
              Data files are specified as ``astropy/pkgname/data/file.dat`` .
            * If a matching local file does not exist, the Astropy data server
              will be queried for the file.
            * A hash like that produced by `compute_hash` can be requested, 
              prefixed by 'hash/' e.g. 'hash/395dd6493cc584df1e78b474fb150840'.
              The hash will first be searched for locally, and if not found, the
              Astropy data server will be queried.
            * A URL to some other file.

    cache : bool
        If True, the file will be downloaded and saved locally or the
        already-cached local copy will be accessed. If False, the file-like
        object will directly access the resource (e.g. if a remote URL is
        accessed, an object like that from `urllib2.urlopen` is returned).
    
    Returns
    -------
    fileobj : file-like
        An object with the contents of the data file available via :func:`read`.
        Can be used as part of a ``with`` statement, automatically closing
        itself after the ``with`` block.
    
    Raises
    ------
    urllib2.URLError
        If a remote file cannot be found.
    IOError
        If problems occur writing or reading a local file.
        
    Examples
    --------
    
    This will retrieve a data file and its contents for the `astropy.wcs` 
    tests:
    
        from astropy.config import get_data_fileobj
        
        fobj = get_data_fileobj('astropy/wcs/tests/data/3d_cd.hdr')
        
        try:
            fcontents = fobj.read()
        finally:
            fobj.close()
            
            
    This downloads a data file and its contents from a specified URL, and does
    *not* cache it remotely:
    
        from astropy.config import get_data_fileobj
        
        vegaurl = 'ftp://ftp.stsci.edu/cdbs/grid/k93models/standards/vega.fits'
        with get_data_fileobj(vegaurl,False) as fobj:
            fcontents = fobj.read()
            
    """
    from urlparse import urlparse
    from urllib2 import urlopen
    from types import MethodType
    
    if cache:
        return open(get_data_filename(dataname), 'r')
    else:
        url = urlparse(dataname)
        if url.scheme != '':
            #it's actually a url for a net location
            urlres = urlopen(dataname,timeout=REMOTE_TIMEOUT())
        else:
            datafn = _find_pkg_data_fn(dataname)
            if datafn is None:
                #not local file - need to get remote data
                urlres = urlopen(DATAURL() + datafn,timeout=REMOTE_TIMEOUT())
            else:
                return open(datafn, 'r')
        
        #need to add in context managers to support with urlopen
        urlres.__enter__ = MethodType(_fake_enter, urlres)
        urlres.__exit__ = MethodType(_fake_exit, urlres)
        
        return urlres

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
        Name/location of the desired data file.  One of the following:

            * The name of a data file included in the source distribution.
              Data files are specified as ``astropy/pkgname/data/file.dat`` .
            * If a matching local file does not exist, the Astropy data server
              will be queried for the file.
            * A hash like that produced by `compute_hash` can be requested, 
              prefixed by 'hash/' e.g. 'hash/395dd6493cc584df1e78b474fb150840'.
              The hash will first be searched for locally, and if not found, the
              Astropy data server will be queried.
            * A URL to some other file.
        
    Raises
    ------
    urllib2.URLError
        If a remote file cannot be found.
    IOError
        If problems occur writing or reading a local file.
    
    Returns
    -------
    filename : str
        A file path on the local file system corresponding to the data 
        requested in `dataname`.
        
    Examples
    --------
    
    This will retrieve the contents of the data file for the `astropy.wcs` 
    tests:
    
        from astropy.config import get_data_filename
        
        fn = get_data_filename('astropy/wcs/tests/data/3d_cd.hdr')
        with open(fn) as f:
            fcontents = f.read()
            
            
    This retrieves a data file by hash either locally or from the astropy data
    server:
    
        from astropy.config import get_data_filename
        
        fn = get_data_filename('hash/da34a7b07ef153eede67387bf950bb32')
        with open(fn) as f:
            fcontents = f.read()
    """
    from urlparse import urlparse
    
    url = urlparse(dataname)
    if url.scheme != '':
        #it's actually a url for a net location
        return _cache_remote(dataname)
    elif dataname.startswith('hash/'):
            #first try looking for a local version if a hash is specified
            hashfn = _find_hash_fn(dataname[5:])
            if hashfn is None:
                return _cache_remote(DATAURL() + dataname)
            else:
                return hashfn
    else:
        datafn = _find_pkg_data_fn(dataname)
        if datafn is None:
            #look for the file on the data server
            return _cache_remote(DATAURL() + dataname)
        else:
            return datafn



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
        block = f.read(COMPUTE_HASH_BLOCK_SIZE())
        while block!='':
            h.update(block)
            block = f.read(COMPUTE_HASH_BLOCK_SIZE())
    
    return h.hexdigest()


def _find_pkg_data_fn(dataname):
    """
    Look for data in the source-included data directories and return the 
    filename or None if the file is not found.
    
    Note that this imports any intermediate packages it finds.
    """
    from os import sep
    from os.path import split,isfile,isdir
    from pkgutil import find_loader
    
    #first split off the file itself
    pkgloc = split(dataname)[0]
    
    #now try to locate the module/package until one is found
    while find_loader(pkgloc) is None:
        pkgloc = split(pkgloc)[0]
        
    #if it's a file, go back one more because it means there was something like
    #/data.py present in the same place as a /data/ dir
    if isfile(pkgloc):
        pkgloc = split(dataname)[0]
    
    path = find_loader(pkgloc).filename + dataname[len(pkgloc):]
    
    if isdir(path):
        raise IOError("Tried to access a data file that's actually a package data directory")
    elif isfile(path):
        return path
    else:
        return None

def _find_hash_fn(hash):
    """
    Looks for a local file by hash - returns file name if found, otherwise
    returns None.
    """
    from os.path import exists,join
    
    dldir,urlmapfn = _get_data_cache_locs()
    hashfn = join(dldir,hash)
    if exists(hashfn):
        return hashfn
    else:
        return None

def _cache_remote(remoteurl):
    """
    Accepts a URL, downloads and caches the result returning the filename, with
    a name determined by the file's MD5 hash. If present in the cache, just
    returns the filename.
    """
    from contextlib import closing
    from os.path import join
    from shutil import move
    from urllib2 import urlopen
    import hashlib
    import shelve
    
    dldir,urlmapfn = _get_data_cache_locs()
    _acquire_data_cache_lock()
    try:
        with closing(shelve.open(urlmapfn)) as url2hash:
            if str(remoteurl) in url2hash:
                localpath = url2hash[str(remoteurl)]
            else:
                with closing(urlopen(remoteurl,timeout=REMOTE_TIMEOUT())) as remote:
                    #save the file to a temporary file
                    tmpfn = join(dldir,'cachedl.tmp')
                    hash = hashlib.md5()
                    
                    with open(tmpfn,'w') as f:
                        block = remote.read(DATA_CACHE_DL_BLOCK_SIZE())
                        while block!='':
                            #TODO: add in something to update a progress bar 
                            #when the download occurs.  Or use the log. either
                            #requires stuff that isn't yet in master
                            f.write(block)
                            hash.update(block)
                            block = remote.read(DATA_CACHE_DL_BLOCK_SIZE())
                    
                    localpath = join(dldir,hash.hexdigest())
                    move(tmpfn,localpath)
                    url2hash[str(remoteurl)] = localpath
    finally:
        _release_data_cache_lock()
    
    return localpath


def clear_data_cache(hashorurl=None):
    """ Clears the data file cache by deleting the local file(s).
    
    Parameters
    ----------
    hashorurl : str or None
        If None, the whole cache is cleared.  Otherwise, either specifies a hash
        for the cached file that is supposed to be deleted, or a URL that has
        previously been downloaded to the cache.
        
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
    
    dldir,urlmapfn = _get_data_cache_locs()
    _acquire_data_cache_lock()
    try:
        
        if hashorurl is None:
            if exists(dldir):
                rmtree(dldir)
            if exists(urlmapfn):
                unlink(urlmapfn)
        else:
            with closing(shelve.open(urlmapfn)) as url2hash:
                filepath = join(dldir, hashorurl)
                
                if exists(filepath):
                    for k,v in url2hash.items():
                        if v==filepath:
                            del url2hash[k]
                    unlink(filepath)
                    
                elif hashorurl in url2hash:
                    filepath = url2hash[hashorurl]
                    del url2hash[hashorurl]
                    unlink(filepath)
                    
                else:
                    msg = 'Could not find file or url {0}'
                    raise OSError(msg.format(hashorurl))
    finally:
        _release_data_cache_lock()

def _get_data_cache_locs():
    """ Finds the path to the data cache directory.
    
    Returns
    -------
    datadir : str
        The path to the data cache directory.
    shelveloc : str
        The path to the shelve object that stores the cache info.
    """
    import shelve
    from .paths import get_cache_dir
    from os.path import exists,isdir,join
    from os import mkdir
    
    datadir = join(get_cache_dir(), 'data')
    shelveloc = join(get_cache_dir(), 'data_urlmap')
    if not exists(datadir):
        mkdir(datadir)
    elif not isdir(datadir):
        msg = 'Data cache directory {0} is not a directory'
        raise IOError(msg.format(datadir))
        
    if isdir(shelveloc):
        msg = 'Data cache shelve object location {0} is a directory'
        raise IOError(msg.format(shelveloc))
    
    return datadir,shelveloc

#the cache directory must be locked before any writes are performed.  Same for
#the hash shelve, so this should be used for both.  Note
def _acquire_data_cache_lock():
    from os.path import join
    from os import mkdir
    from time import sleep
    
    lockdir = join(_get_data_cache_locs()[0],'lock')
    for i in range(DATA_CACHE_LOCK_ATTEMPTS()):
        try:
            mkdir(lockdir)
        except OSError:
            sleep(1)
        else:
            return
    msg = 'Unable to acquire lock for cache directory ({0} exists)'
    raise RuntimeError(msg.format(lockdir))

def _release_data_cache_lock():
    from os.path import join,exists,isdir
    from os import rmdir
    
    lockdir = join(_get_data_cache_locs()[0],'lock')
    
    if exists(lockdir) and isdir(lockdir):
        rmdir(lockdir)
    else:
        msg = 'Error releasing lock. "{0}" either does not exist or is not '+\
              'a directory.'
        raise RuntimeError(msg.format(lockdir))
