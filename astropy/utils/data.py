# Licensed under a 3-clause BSD style license - see LICENSE.rst
""" This module contains helper functions for accessing, downloading, and
caching data files.
"""

from __future__ import division

import os
import io
import sys
import atexit
import contextlib
import urllib2

from ..config.configuration import ConfigurationItem

__all__ = ['get_fileobj', 'get_pkg_data_fileobj', 'get_pkg_data_filename',
           'get_pkg_data_contents', 'get_pkg_data_fileobjs',
           'get_pkg_data_filenames', 'compute_hash',
           'clear_data_cache', 'CacheMissingWarning', 'cache_remote']

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


PY3K = sys.version_info[0] >= 3


if not PY3K:  # pragma: py2
    #used for supporting with statements in get_pkg_data_fileobj
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


def _is_url(string):
    """
    Test whether a string is a valid URL

    Parameters
    ----------
    string : str
        The string to test
    """
    from urlparse import urlparse
    url = urlparse(string)
    return url[0] != '' # url[0]==url.scheme, but url[0] is py 2.6-compat


@contextlib.contextmanager
def get_fileobj(name_or_obj, encoding=None, cache=False):
    """
    Given a filename or a readable file-like object, return a readable
    file-like object.

    This supports passing filenames, URLs, and readable file-like
    objects, any of which can be compressed in gzip or bzip2.

    Parameters
    ----------
    name_or_obj : str or file-like object
        The filename of the file to access (if given as a string), or
        the file-like object to access.

        If a file-like object, it must be opened in binary mode.

    encoding : str, optional
        When `None` (default), returns a file-like object with a
        `read` method that on Python 2.x returns `bytes` objects and
        on Python 3.x returns `str` (`unicode`) objects, using
        `locale.getpreferredencoding()` as an encoding.  This matches
        the default behavior of the built-in `open` when no `mode`
        argument is provided.

        When `'binary'`, returns a file-like object where its `read`
        method returns `bytes` objects.

        When another string, it is the name of an encoding, and the
        file-like object's `read` method will return `str` (`unicode`)
        objects, decoded from binary using the given encoding.

    cache : bool, optional
        Whether to cache the contents of remote URLs
    """
    close_fds = []

    # Get a file object to the content
    if isinstance(name_or_obj, basestring):
        if _is_url(name_or_obj):
            if cache:
                fileobj = open(cache_remote(name_or_obj))
            else:
                fileobj = urllib2.urlopen(name_or_obj, timeout=REMOTE_TIMEOUT())
                close_fds.append(fileobj)
                from types import MethodType
                if not PY3K:  # pragma: py2
                    # Need to add in context managers to support with urlopen for <3.x
                    fileobj.__enter__ = MethodType(_fake_enter, fileobj)
                    fileobj.__exit__ = MethodType(_fake_exit, fileobj)
        else:
            if PY3K:
                fileobj = io.FileIO(name_or_obj, 'r')
            else:
                fileobj = open(name_or_obj, 'rb')
            close_fds.append(fileobj)
    else:
        fileobj = name_or_obj

    # Check if the file object supports random access, and if not, then wrap
    # it in a StringIO buffer.
    if not hasattr(fileobj, 'seek'):
        from StringIO import StringIO
        fileobj = StringIO(fileobj.read())

    # Now read enough bytes to look at signature
    signature = fileobj.read(4)
    fileobj.seek(0)

    if signature[:3] == b'\x1f\x8b\x08':  # gzip
        import struct
        try:
            from .compat import gzip
            fileobj_new = gzip.GzipFile(fileobj=fileobj, mode='rb')
            fileobj_new.read(1)  # need to check that the file is really gzip
        except IOError:  # invalid gzip file
            fileobj.seek(0)
        except struct.error:  # invalid gzip file on Python 3
            fileobj.seek(0)
        else:
            fileobj_new.seek(0)
            fileobj = fileobj_new
    elif signature[:3] == b'BZh':  # bzip2
        try:
            # bz2.BZ2File does not support file objects, only filenames, so we
            # need to write the data to a temporary file
            import tempfile
            tmp = tempfile.NamedTemporaryFile()
            tmp.write(fileobj.read())
            tmp.flush()
            close_fds.append(tmp)
            import bz2
            fileobj_new = bz2.BZ2File(tmp.name, mode='rb')
            fileobj_new.read(1)  # need to check that the file is really bzip2
        except IOError:  # invalid bzip2 file
            fileobj.seek(0)
        else:
            fileobj_new.seek(0)
            fileobj = fileobj_new

    # By this point, we have a file, io.FileIO, gzip.GzipFile, or
    # bz2.BZ2File instance opened in binary mode (that is, read
    # returns bytes).  Now we need to, if requested, wrap it in a
    # io.TextIOWrapper so read will return unicode based on the
    # encoding parameter.

    if PY3K:
        needs_textio_wrapper = encoding != 'binary'
    else:
        needs_textio_wrapper = encoding != 'binary' and encoding is not None

    if needs_textio_wrapper:
        import bz2
        # FIXME: A bz2.BZ2File can not be wrapped by a TextIOWrapper,
        # so on Python 3 the user will get back bytes from the file
        # rather than Unicode as expected.
        if not isinstance(fileobj, bz2.BZ2File):
            # On Python 2.x, we need to first wrap the regular `file`
            # instance in a `io.FileIO` object before it can be
            # wrapped in a `TextIOWrapper`.  We don't just create an
            # `io.FileIO` object in the first place, because we can't
            # get a raw file descriptor out of it on Python 2.x, which
            # is required for the XML iterparser.
            if not PY3K and isinstance(fileobj, file):
                fileobj = io.FileIO(fileobj.fileno())

            fileobj = io.BufferedReader(fileobj)
            fileobj = io.TextIOWrapper(fileobj, encoding=encoding)

    yield fileobj

    for fd in close_fds:
        fd.close()


def get_pkg_data_fileobj(data_name, cache=True, encoding=None):
    """
    Retrieves a data file from the standard locations for the package and
    provides the file as a file-like object that reads bytes.

    Parameters
    ----------
    data_name : str
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

    cache : bool
        If True, the file will be downloaded and saved locally or the
        already-cached local copy will be accessed. If False, the file-like
        object will directly access the resource (e.g. if a remote URL is
        accessed, an object like that from `urllib2.urlopen` is returned).

    encoding : str, optional
        When `None` (default), returns a file-like object with a
        `read` method that on Python 2.x returns `bytes` objects and
        on Python 3.x returns `str` (`unicode`) objects, using
        `locale.getpreferredencoding()` as an encoding.  This matches
        the default behavior of the built-in `open` when no `mode`
        argument is provided.

        When `'binary'`, returns a file-like object where its `read`
        method returns `bytes` objects.

        When another string, it is the name of an encoding, and the
        file-like object's `read` method will return `str` (`unicode`)
        objects, decoded from binary using the given encoding.


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

        from astropy.config import get_pkg_data_fileobj

        with get_pkg_data_fileobj('data/3d_cd.hdr') as fobj:
            fcontents = fobj.read()


    This downloads a data file and its contents from a specified URL, and does
    *not* cache it remotely::

        from astropy.config import get_pkg_data_fileobj

        vegaurl = 'ftp://ftp.stsci.edu/cdbs/grid/k93models/standards/vega.fits'
        with get_pkg_data_fileobj(vegaurl,False) as fobj:
            fcontents = fobj.read()

    See Also
    --------
    get_pkg_data_contents : returns the contents of a file or url as a bytes object
    get_pkg_data_filename : returns a local name for a file containing the data
    """
    datafn = _find_pkg_data_path(data_name)
    if os.path.isdir(datafn):
        raise IOError("Tried to access a data file that's actually "
                      "a package data directory")
    elif os.path.isfile(datafn):  # local file
        return get_fileobj(datafn, encoding=encoding)
    else:  # remote file
        return get_fileobj(DATAURL() + datafn, encoding=encoding)


def get_pkg_data_filename(data_name):
    """
    Retrieves a data file from the standard locations for the package and
    provides a local filename for the data.

    This function is similar to `get_pkg_data_fileobj` but returns the file *name*
    instead of a readable file-like object.  This means that this function must
    always cache remote files locally, unlike `get_pkg_data_fileobj`.

    Parameters
    ----------
    data_name : str
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
        requested in `data_name`.

    Examples
    --------

    This will retrieve the contents of the data file for the `astropy.wcs`
    tests::

        from astropy.config import get_pkg_data_filename

        fn = get_pkg_data_filename('data/3d_cd.hdr')
        with open(fn) as f:
            fcontents = f.read()


    This retrieves a data file by hash either locally or from the astropy data
    server::

        from astropy.config import get_pkg_data_filename

        fn = get_pkg_data_filename('hash/da34a7b07ef153eede67387bf950bb32')
        with open(fn) as f:
            fcontents = f.read()

    See Also
    --------
    get_pkg_data_contents : returns the contents of a file or url as a bytes object
    get_pkg_data_fileobj : returns a file-like object with the data
    """

    if data_name.startswith('hash/'):
        # first try looking for a local version if a hash is specified
        hashfn = _find_hash_fn(data_name[5:])
        if hashfn is None:
            return cache_remote(DATAURL() + data_name)
        else:
            return hashfn
    else:
        datafn = _find_pkg_data_path(data_name)
        if os.path.isdir(datafn):
            raise IOError("Tried to access a data file that's actually "
                          "a package data directory")
        elif os.path.isfile(datafn):  # local file
            return datafn
        else:  # remote file
            return cache_remote(DATAURL() + data_name)


def get_pkg_data_contents(data_name, cache=True, encoding=None):
    """
    Retrieves a data file from the standard locations and returns its
    contents as a bytes object.

    Parameters
    ----------
    data_name : str
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

    encoding : str, optional
        When `None` (default), returns a file-like object with a
        `read` method that on Python 2.x returns `bytes` objects and
        on Python 3.x returns `str` (`unicode`) objects, using
        `locale.getpreferredencoding()` as an encoding.  This matches
        the default behavior of the built-in `open` when no `mode`
        argument is provided.

        When `'binary'`, returns a file-like object where its `read`
        method returns `bytes` objects.

        When another string, it is the name of an encoding, and the
        file-like object's `read` method will return `str` (`unicode`)
        objects, decoded from binary using the given encoding.


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
    get_pkg_data_fileobj : returns a file-like object with the data
    get_pkg_data_filename : returns a local name for a file containing the data
    """
    with get_pkg_data_fileobj(data_name, cache=cache, encoding=encoding) as fd:
        contents = fd.read()
    return contents


def get_pkg_data_filenames(datadir, pattern='*'):
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

        from astropy.config import get_pkg_data_filenames

        for fn in get_pkg_data_filename('maps', '*.hdr'):
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


def get_pkg_data_fileobjs(datadir, pattern='*', encoding=None):
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

    encoding : str, optional
        When `None` (default), returns a file-like object with a
        `read` method that on Python 2.x returns `bytes` objects and
        on Python 3.x returns `str` (`unicode`) objects, using
        `locale.getpreferredencoding()` as an encoding.  This matches
        the default behavior of the built-in `open` when no `mode`
        argument is provided.

        When `'binary'`, returns a file-like object where its `read`
        method returns `bytes` objects.

        When another string, it is the name of an encoding, and the
        file-like object's `read` method will return `str` (`unicode`)
        objects, decoded from binary using the given encoding.


    Returns
    -------
    fileobjs : iterator of file objects
        File objects for each of the files on the local filesystem in
        *datadir* matching *pattern*.

    Examples
    --------
    This will retrieve the contents of the data file for the `astropy.wcs`
    tests::

        from astropy.config import get_pkg_data_filenames

        for fd in get_pkg_data_filename('maps', '*.hdr'):
            fcontents = fd.read()
    """
    for fn in get_pkg_data_filenames(datadir, pattern):
        with get_pkg_data_fileobj(fn, encoding=encoding) as fd:
            yield fd


def compute_hash(localfn):
    """ Computes the MD5 hash for a file.

    The hash for a data file is used for looking up data files in a unique
    fashion. This is of particular use for tests; a test may require a
    particular version of a particular file, in which case it can be accessed
    via hash to get the appropriate version.

    Typically, if you wish to write a test that requires a particular data
    file, you will want to submit that file to the astropy data servers, and
    use e.g. ``get_pkg_data_filename('hash/a725fa6ba642587436612c2df0451956')``,
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


def _find_pkg_data_path(data_name):
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
    path = join(module_path, data_name)

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


def cache_remote(remoteurl):
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

        with closing(urllib2.urlopen(remoteurl, timeout=REMOTE_TIMEOUT())) as remote:
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


#this is used by cache_remote and _deltemps to determine the files to delete
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
    from ..config.paths import get_cache_dir
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

    if not PY3K:  # pragma: py3
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
