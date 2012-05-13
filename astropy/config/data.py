# Licensed under a 3-clause BSD style license - see LICENSE.rst
""" This module contains helper functions for accessing, downloading, and
caching data files.
"""

from __future__ import division

import sys
import atexit

from .configuration import ConfigurationItem


__all__ = ['get_data_fileobj', 'get_data_filename', 'get_data_contents',
           'get_data_fileobjs', 'get_data_filenames', 'compute_hash',
           'clear_data_cache', 'CacheMissingWarning']

DATAURL = ConfigurationItem(
    'dataurl', 'http://data.astropy.org/', 'URL for astropy remote data site.')
REMOTE_TIMEOUT = ConfigurationItem(
    'remote_timeout', 3., 'Time to wait for remote data query (in seconds).')
COMPUTE_HASH_BLOCK_SIZE = ConfigurationItem(
    'hash_block_size', 2 ** 16,  # 64K
    'Block size for computing MD5 file hashes.')
DATA_CACHE_DL_BLOCK_SIZE = ConfigurationItem(
    'data_dl_block_size', 2 ** 16,  # 64K
    'Number of bytes of remote data to download per step.')
DATA_CACHE_LOCK_ATTEMPTS = ConfigurationItem(
    'data_cache_lock_attempts', 3, 'Number of times to try to get the lock ' +
    'while accessing the data cache before giving up.')
DELETE_TEMPORARY_DOWNLOADS_AT_EXIT = ConfigurationItem(
    'delete_temporary_downloads_at_exit', True, 'If True, temporary download' +
    ' files created when the cache is inacessible will be deleted at the end' +
    ' of the python session.')


if sys.version_info[0] < 3:  # pragma: py2
    #used for supporting with statements in get_data_fileobj
    def _fake_enter(self):
        return self

    def _fake_exit(self, type, value, traceback):
        self.close()


class CacheMissingWarning(Warning):
    """
    This warning indicates the standard cache directory is not accessible, with
    the first argument providing the warning message. If args[1] is present, it
    is a filename indicating the path to a temporary file that was created to
    store a remote data download in the absence of the cache.
    """


def get_data_fileobj(dataname, cache=True):
    """
    Retrieves a data file from the standard locations and provides the file as
    a file-like object.

    Parameters
    ----------
    dataname : str
        Name/location of the desired data file.  One of the following:

            * The name of a data file included in the source
              distribution.  The path is relative to the module
              calling this function.  For example, if calling from
              `astropy.pkname`, use ``'data/file.dat'`` to get the
              file in ``astropy/pkgname/data/file.dat``.  Double-dots
              can be used to go up a level.  In the same example, use
              ``'../data/file.dat'`` to get ``astropy/data/file.dat``.
            * If a matching local file does not exist, the Astropy
              data server will be queried for the file.
            * A hash like that produced by `compute_hash` can be
              requested, prefixed by 'hash/'
              e.g. 'hash/395dd6493cc584df1e78b474fb150840'.  The hash
              will first be searched for locally, and if not found,
              the Astropy data server will be queried.
            * A URL to some other file.

    cache : bool
        If True, the file will be downloaded and saved locally or the
        already-cached local copy will be accessed. If False, the file-like
        object will directly access the resource (e.g. if a remote URL is
        accessed, an object like that from `urllib2.urlopen` is returned).

    Returns
    -------
    fileobj : file-like
        An object with the contents of the data file available via
        :func:`read`.  Can be used as part of a ``with`` statement,
        automatically closing itself after the ``with`` block.

    Raises
    ------
    urllib2.URLError
        If a remote file cannot be found.
    IOError
        If problems occur writing or reading a local file.

    Examples
    --------

    This will retrieve a data file and its contents for the `astropy.wcs`
    tests::

        from astropy.config import get_data_fileobj

        with get_data_fileobj('data/3d_cd.hdr') as fobj:
            fcontents = fobj.read()


    This downloads a data file and its contents from a specified URL, and does
    *not* cache it remotely::

        from astropy.config import get_data_fileobj

        vegaurl = 'ftp://ftp.stsci.edu/cdbs/grid/k93models/standards/vega.fits'
        with get_data_fileobj(vegaurl,False) as fobj:
            fcontents = fobj.read()

    See Also
    --------
    get_data_contents : returns the contents of a file or url as a bytes object
    get_data_filename : returns a local name for a file containing the data
    """
    from os.path import isdir, isfile
    from urlparse import urlparse
    from urllib2 import urlopen
    from types import MethodType

    if cache:
        return open(get_data_filename(dataname), 'rb')
    else:
        url = urlparse(dataname)
        if url[0] != '':  # url[0]==url.scheme, but url[0] is py 2.6-compat
            #it's actually a url for a net location
            urlres = urlopen(dataname, timeout=REMOTE_TIMEOUT())
        else:
            datafn = _find_pkg_data_path(dataname)
            if isdir(datafn):
                raise IOError(
                    "Tried to access a data file that's actually "
                    "a package data directory")
            elif isfile(datafn):
                return open(datafn, 'rb')
            else:
                #not local file - need to get remote data
                urlres = urlopen(DATAURL() + datafn, timeout=REMOTE_TIMEOUT())

        if sys.version_info[0] < 3:  # pragma: py2
            #need to add in context managers to support with urlopen for <3.x
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

            * The name of a data file included in the source
              distribution.  The path is relative to the module
              calling this function.  For example, if calling from
              `astropy.pkname`, use ``'data/file.dat'`` to get the
              file in ``astropy/pkgname/data/file.dat``.  Double-dots
              can be used to go up a level.  In the same example, use
              ``'../data/file.dat'`` to get ``astropy/data/file.dat``.
            * If a matching local file does not exist, the Astropy
              data server will be queried for the file.
            * A hash like that produced by `compute_hash` can be
              requested, prefixed by 'hash/'
              e.g. 'hash/395dd6493cc584df1e78b474fb150840'.  The hash
              will first be searched for locally, and if not found,
              the Astropy data server will be queried.
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
    tests::

        from astropy.config import get_data_filename

        fn = get_data_filename('data/3d_cd.hdr')
        with open(fn) as f:
            fcontents = f.read()


    This retrieves a data file by hash either locally or from the astropy data
    server::

        from astropy.config import get_data_filename

        fn = get_data_filename('hash/da34a7b07ef153eede67387bf950bb32')
        with open(fn) as f:
            fcontents = f.read()

    See Also
    --------
    get_data_contents : returns the contents of a file or url as a bytes object
    get_data_fileobj : returns a file-like object with the data
    """
    from os.path import isdir, isfile
    from urlparse import urlparse

    url = urlparse(dataname)
    if url[0] != '':  # url[0]==url.scheme, but url[0] is py 2.6-compat
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
        datafn = _find_pkg_data_path(dataname)
        if isdir(datafn):
            raise IOError(
                "Tried to access a data file that's actually "
                "a package data directory")
        elif isfile(datafn):
            return datafn
        else:
            #look for the file on the data server
            return _cache_remote(DATAURL() + dataname)


def get_data_contents(dataname, cache=True):
    """
    Retrieves a data file from the standard locations and returns its
    contents as a bytes object.

    Parameters
    ----------
    dataname : str
        Name/location of the desired data file.  One of the following:

            * The name of a data file included in the source
              distribution.  The path is relative to the module
              calling this function.  For example, if calling from
              `astropy.pkname`, use ``'data/file.dat'`` to get the
              file in ``astropy/pkgname/data/file.dat``.  Double-dots
              can be used to go up a level.  In the same example, use
              ``'../data/file.dat'`` to get ``astropy/data/file.dat``.
            * If a matching local file does not exist, the Astropy
              data server will be queried for the file.
            * A hash like that produced by `compute_hash` can be
              requested, prefixed by 'hash/'
              e.g. 'hash/395dd6493cc584df1e78b474fb150840'.  The hash
              will first be searched for locally, and if not found,
              the Astropy data server will be queried.
            * A URL to some other file.

    cache : bool
        If True, the file will be downloaded and saved locally or the
        already-cached local copy will be accessed. If False, the file-like
        object will directly access the resource (e.g. if a remote URL is
        accessed, an object like that from `urllib2.urlopen` is returned).

    Returns
    -------
    contents : bytes
        The complete contents of the file as a bytes object.

    Raises
    ------
    urllib2.URLError
        If a remote file cannot be found.
    IOError
        If problems occur writing or reading a local file.

    See Also
    --------
    get_data_fileobj : returns a file-like object with the data
    get_data_filename : returns a local name for a file containing the data
    """
    with get_data_fileobj(dataname, cache=cache) as fd:
        contents = fd.read()
    return contents


def get_data_filenames(datadir, pattern='*'):
    """
    Returns the path of all of the data files in a given directory
    that match a given glob pattern.

    Parameters
    ----------
    datadir : str
        Name/location of the desired data files.  One of the following:

            * The name of a directory included in the source
              distribution.  The path is relative to the module
              calling this function.  For example, if calling from
              `astropy.pkname`, use ``'data'`` to get the
              files in ``astropy/pkgname/data``.
            * Remote URLs are not currently supported.

    pattern : str, optional
        A UNIX-style filename glob pattern to match files.  See the
        `glob` module in the standard library for more information.
        By default, matches all files.

    Returns
    -------
    filenames : iterator of str
        Paths on the local filesystem in *datadir* matching *pattern*.

    Examples
    --------
    This will retrieve the contents of the data file for the `astropy.wcs`
    tests::

        from astropy.config import get_data_filenames

        for fn in get_data_filename('maps', '*.hdr'):
            with open(fn) as f:
                fcontents = f.read()
    """
    import fnmatch
    from os.path import isdir, isfile, join
    from os import listdir

    path = _find_pkg_data_path(datadir)
    if isfile(path):
        raise IOError(
            "Tried to access a data directory that's actually "
            "a package data file")
    elif isdir(path):
        for filename in listdir(path):
            if fnmatch.fnmatch(filename, pattern):
                yield join(path, filename)
    else:
        raise IOError("Path not found")


def get_data_fileobjs(datadir, pattern='*'):
    """
    Returns readable file objects for all of the data files in a given
    directory that match a given glob pattern.

    Parameters
    ----------
    datadir : str
        Name/location of the desired data files.  One of the following:

            * The name of a directory included in the source
              distribution.  The path is relative to the module
              calling this function.  For example, if calling from
              `astropy.pkname`, use ``'data'`` to get the
              files in ``astropy/pkgname/data``
            * Remote URLs are not currently supported

    pattern : str, optional
        A UNIX-style filename glob pattern to match files.  See the
        `glob` module in the standard library for more information.
        By default, matches all files.

    Returns
    -------
    fileobjs : iterator of file objects
        File objects for each of the files on the local filesystem in
        *datadir* matching *pattern*.

    Examples
    --------
    This will retrieve the contents of the data file for the `astropy.wcs`
    tests::

        from astropy.config import get_data_filenames

        for fd in get_data_filename('maps', '*.hdr'):
            fcontents = fd.read()
    """
    for fn in get_data_filenames(datadir, pattern):
        with open(fn, 'rb') as fd:
            yield fd


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

    with open(localfn, 'rb') as f:
        h = hashlib.md5()
        block = f.read(COMPUTE_HASH_BLOCK_SIZE())
        while block:
            h.update(block)
            block = f.read(COMPUTE_HASH_BLOCK_SIZE())

    return h.hexdigest()


def _find_pkg_data_path(dataname):
    """
    Look for data in the source-included data directories and return the
    path.
    """
    from os.path import abspath, dirname, join
    from ..utils.misc import find_current_module

    module = find_current_module(1, True)
    rootpkgname = module.__package__.split('.')[0]
    rootpkg = __import__(rootpkgname)

    module_path = dirname(module.__file__)
    path = join(module_path, dataname)

    root_dir = dirname(rootpkg.__file__)
    assert abspath(path).startswith(abspath(root_dir)), \
           ("attempted to get a local data file outside "
            "of the " + rootpkgname + " tree")

    return path


def _find_hash_fn(hash):
    """
    Looks for a local file by hash - returns file name if found and a valid
    file, otherwise returns None.
    """
    from os.path import isfile, join
    from warnings import warn

    try:
        dldir, urlmapfn = _get_data_cache_locs()
    except (IOError, OSError) as e:
        msg = 'Could not access cache directory to search for data file: '
        warn(CacheMissingWarning(msg + str(e)))
        return None
    hashfn = join(dldir, hash)
    if isfile(hashfn):
        return hashfn
    else:
        return None


def _cache_remote(remoteurl):
    """
    Accepts a URL, downloads and caches the result returning the filename, with
    a name determined by the file's MD5 hash. If present in the cache, just
    returns the filename.
    """
    import hashlib
    from contextlib import closing
    from os.path import join
    from tempfile import NamedTemporaryFile
    from shutil import move
    from urllib2 import urlopen
    from warnings import warn

    from ..utils.console import ProgressBarOrSpinner

    try:
        dldir, urlmapfn = _get_data_cache_locs()
        docache = True
    except (IOError, OSError) as e:
        msg = 'Remote data cache could not be accessed due to '
        estr = '' if len(e.args) < 1 else (': ' + str(e))
        warn(CacheMissingWarning(msg + e.__class__.__name__ + estr))
        docache = False

    if docache:
        _acquire_data_cache_lock()
    try:
        if docache:
            with _open_shelve(urlmapfn, True) as url2hash:
                if str(remoteurl) in url2hash:
                    return url2hash[str(remoteurl)]

        with closing(urlopen(remoteurl, timeout=REMOTE_TIMEOUT())) as remote:
            #keep a hash to rename the local file to the hashed name
            hash = hashlib.md5()

            info = remote.info()
            if 'Content-Length' in info:
                try:
                    size = int(info['Content-Length'])
                except ValueError:
                    size = None
            else:
                size = None

            dlmsg = "Downloading {0}".format(remoteurl)
            with ProgressBarOrSpinner(size, dlmsg) as p:
                with NamedTemporaryFile(delete=False) as f:
                    bytes_read = 0
                    block = remote.read(DATA_CACHE_DL_BLOCK_SIZE())
                    while block:
                        f.write(block)
                        hash.update(block)
                        bytes_read += len(block)
                        p.update(bytes_read)
                        block = remote.read(DATA_CACHE_DL_BLOCK_SIZE())
        if docache:
            with _open_shelve(urlmapfn, True) as url2hash:
                localpath = join(dldir, hash.hexdigest())
                move(f.name, localpath)
                url2hash[str(remoteurl)] = localpath
        else:
            localpath = f.name
            msg = 'File downloaded to temp file due to lack of cache access.'
            warn(CacheMissingWarning(msg, localpath))
            if DELETE_TEMPORARY_DOWNLOADS_AT_EXIT():
                global _tempfilestodel
                _tempfilestodel.append(localpath)

    finally:
        if docache:
            _release_data_cache_lock()

    return localpath


#this is used by _cache_remote and _deltemps to determine the files to delete
# when the interpreter exits
_tempfilestodel = []


@atexit.register
def _deltemps():
    import os

    global _tempfilestodel

    while len(_tempfilestodel) > 0:
        fn = _tempfilestodel.pop()
        if os.path.isfile(fn):
            os.remove(fn)


def clear_data_cache(hashorurl=None):
    """ Clears the data file cache by deleting the local file(s).

    Parameters
    ----------
    hashorurl : str or None
        If None, the whole cache is cleared.  Otherwise, either specifies a
        hash for the cached file that is supposed to be deleted, or a URL that
        has previously been downloaded to the cache.

    Raises
    ------
    OSEerror
        If the requested filename is not present in the data directory.

    """
    from os import unlink
    from os.path import join, exists, abspath
    from shutil import rmtree
    from warnings import warn

    try:
        dldir, urlmapfn = _get_data_cache_locs()
    except (IOError, OSError) as e:
        msg = 'Not clearing data cache - cache inacessable due to '
        estr = '' if len(e.args) < 1 else (': ' + str(e))
        warn(CacheMissingWarning(msg + e.__class__.__name__ + estr))
        return

    _acquire_data_cache_lock()
    try:

        if hashorurl is None:
            if exists(dldir):
                rmtree(dldir)
            if exists(urlmapfn):
                unlink(urlmapfn)
        else:
            with _open_shelve(urlmapfn, True) as url2hash:
                filepath = join(dldir, hashorurl)
                assert abspath(filepath).startswith(abspath(dldir)), \
                       ("attempted to use clear_data_cache on a location" +
                        " that's not inside the data cache directory")

                if exists(filepath):
                    for k, v in url2hash.items():
                        if v == filepath:
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
    from .paths import get_cache_dir
    from os.path import exists, isdir, join
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

    return datadir, shelveloc


def _open_shelve(shelffn, withclosing=False):
    """
    opens a shelf in a way that is py3.x and py2.x compatible.  If
    `withclosing` is  True, it will be opened with closing, allowing use like:

        with _open_shelve('somefile',True) as s:
            ...
    """
    import shelve
    from contextlib import closing

    if sys.version_info[0] > 2:  # pragma: py3
        shelf = shelve.open(shelffn, protocol=2)
    else:  # pragma: py2
        shelf = shelve.open(shelffn + '.db', protocol=2)

    if withclosing:
        return closing(shelf)
    else:
        return shelf


#the cache directory must be locked before any writes are performed.  Same for
#the hash shelve, so this should be used for both.  Note
def _acquire_data_cache_lock():
    from os.path import join
    from os import mkdir
    from time import sleep

    lockdir = join(_get_data_cache_locs()[0], 'lock')
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
    from os.path import join, exists, isdir
    from os import rmdir

    lockdir = join(_get_data_cache_locs()[0], 'lock')

    if exists(lockdir) and isdir(lockdir):
        rmdir(lockdir)
    else:
        msg = 'Error releasing lock. "{0}" either does not exist or is not ' +\
              'a directory.'
        raise RuntimeError(msg.format(lockdir))
