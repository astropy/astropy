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

__all__ = ['get_readable_fileobj', 'get_file_contents', 'get_pkg_data_fileobj',
           'get_pkg_data_filename', 'get_pkg_data_contents',
           'get_pkg_data_fileobjs', 'get_pkg_data_filenames', 'compute_hash',
           'clear_download_cache', 'CacheMissingWarning', 'download_file',
           'download_files_in_parallel']

DATAURL = ConfigurationItem(
    'dataurl', 'http://data.astropy.org/', 'URL for astropy remote data site.')
REMOTE_TIMEOUT = ConfigurationItem(
    'remote_timeout', 3., 'Time to wait for remote data query (in seconds).')
COMPUTE_HASH_BLOCK_SIZE = ConfigurationItem(
    'hash_block_size', 2 ** 16,  # 64K
    'Block size for computing MD5 file hashes.')
DOWNLOAD_CACHE_BLOCK_SIZE = ConfigurationItem(
    'download_block_size', 2 ** 16,  # 64K
    'Number of bytes of remote data to download per step.')
DOWNLOAD_CACHE_LOCK_ATTEMPTS = ConfigurationItem(
    'download_cache_lock_attempts', 5, 'Number of times to try to get the lock ' +
    'while accessing the data cache before giving up.')
DELETE_TEMPORARY_DOWNLOADS_AT_EXIT = ConfigurationItem(
    'delete_temporary_downloads_at_exit', True, 'If True, temporary download' +
    ' files created when the cache is inacessible will be deleted at the end' +
    ' of the python session.')


PY3K = sys.version_info[0] >= 3


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
    # url[0]==url.scheme, but url[0] is py 2.6-compat
    # we can't just check that url[0] is not an empty string, because
    # file paths in windows would return a non-empty scheme (e.g. e:\\
    # returns 'e').
    return url[0].lower() in ['http', 'https', 'ftp']


def _is_inside(path, parent_path):
    # We have to try realpath too to avoid issues with symlinks, but we leave
    # abspath because some systems like debian have the absolute path (with no
    # symlinks followed) match, but the real directories in different
    # locations, so need to try both cases.
    return os.path.abspath(path).startswith(os.path.abspath(parent_path)) \
        or os.path.realpath(path).startswith(os.path.realpath(parent_path))

@contextlib.contextmanager
def get_readable_fileobj(name_or_obj, encoding=None, cache=False):
    """
    Given a filename or a readable file-like object, return a context
    manager that yields a readable file-like object.

    This supports passing filenames, URLs, and readable file-like
    objects, any of which can be compressed in gzip or bzip2.

    Notes
    -----

    This function is a context manager, and should be used for example
    as::

        with get_readable_fileobj('file.dat') as f:
            contents = f.read()

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
        Whether to cache the contents of remote URLs.
    """
    import tempfile

    # close_fds is a list of file handles created by this function
    # that need to be closed.  We don't want to always just close the
    # returned file handle, because it may simply be the file handle
    # passed in.  In that case it is not the responsibility of this
    # function to close it: doing so could result in a "double close"
    # and an "invalid file descriptor" exception.
    close_fds = []
    delete_fds = []

    # Get a file object to the content
    if isinstance(name_or_obj, basestring):
        if _is_url(name_or_obj):
            name_or_obj = download_file(name_or_obj, cache=cache)
        if PY3K:
            fileobj = io.FileIO(name_or_obj, 'r')
        else:
            fileobj = open(name_or_obj, 'rb')
        close_fds.append(fileobj)
    else:
        fileobj = name_or_obj

    # Check if the file object supports random access, and if not,
    # then wrap it in a BytesIO buffer.  It would be nicer to use a
    # BufferedReader to avoid reading loading the whole file first,
    # but that is not compatible with streams or urllib2.urlopen
    # objects on Python 2.x.
    if not hasattr(fileobj, 'seek'):
        fileobj = io.BytesIO(fileobj.read())

    # Now read enough bytes to look at signature
    signature = fileobj.read(4)
    fileobj.seek(0)

    if signature[:3] == b'\x1f\x8b\x08':  # gzip
        import struct
        try:
            from .compat import gzip
            fileobj_new = gzip.GzipFile(fileobj=fileobj, mode='rb')
            fileobj_new.read(1)  # need to check that the file is really gzip
        except (IOError, EOFError):  # invalid gzip file
            fileobj.seek(0)
            fileobj_new.close()
        except struct.error:  # invalid gzip file on Python 3
            fileobj.seek(0)
            fileobj_new.close()
        else:
            fileobj_new.seek(0)
            fileobj = fileobj_new
    elif signature[:3] == b'BZh':  # bzip2
        try:
            # bz2.BZ2File does not support file objects, only filenames, so we
            # need to write the data to a temporary file
            tmp = tempfile.NamedTemporaryFile("wb", delete=False)
            tmp.write(fileobj.read())
            tmp.close()
            delete_fds.append(tmp)
            import bz2
            fileobj_new = bz2.BZ2File(tmp.name, mode='rb')
            fileobj_new.read(1)  # need to check that the file is really bzip2
        except IOError:  # invalid bzip2 file
            fileobj.seek(0)
            fileobj_new.close()
        else:
            fileobj_new.seek(0)
            close_fds.append(fileobj_new)
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
        # A bz2.BZ2File can not be wrapped by a TextIOWrapper,
        # so we decompress it to a temporary file and then
        # return a handle to that.
        import bz2
        if isinstance(fileobj, bz2.BZ2File):
            tmp = tempfile.NamedTemporaryFile("wb", delete=False)
            data = fileobj.read()
            tmp.write(data)
            tmp.close()
            delete_fds.append(tmp)
            if PY3K:
                fileobj = io.FileIO(tmp.name, 'r')
            else:
                fileobj = open(tmp.name, 'rb')
            close_fds.append(fileobj)

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

        # Ensure that file is at the start - io.FileIO will for
        # example not always be at the start:
        # >>> import io
        # >>> f = open('test.fits', 'rb')
        # >>> f.read(4)
        # 'SIMP'
        # >>> f.seek(0)
        # >>> fileobj = io.FileIO(f.fileno())
        # >>> fileobj.tell()
        # 4096L

        fileobj.seek(0)

    try:
        yield fileobj
    finally:
        for fd in close_fds:
            fd.close()
        for fd in delete_fds:
            os.remove(fd.name)


def get_file_contents(name_or_obj, encoding=None, cache=False):
    """
    Retrieves the contents of a filename or file-like object.

    See  the `get_readable_fileobj` docstring for details on parameters.

    Returns
    -------
    content
        The content of the file (as requested by `encoding`).

    """
    with get_readable_fileobj(name_or_obj, encoding, cache) as f:
        return f.read()


def get_pkg_data_fileobj(data_name, encoding=None, cache=True):
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

        from astropy.config import get_pkg_data_fileobj

        with get_pkg_data_fileobj('data/3d_cd.hdr') as fobj:
            fcontents = fobj.read()

    This would download a data file from the astropy data server
    because the ``standards/vega.fits`` file is not present in the
    source distribution.  It will also save the file locally so the
    next time it is accessed it won't need to be downloaded.::

        from astropy.config import get_pkg_data_fileobj

        with get_pkg_data_fileobj('standards/vega.fits') as fobj:
            fcontents = fobj.read()

    This does the same thing but does *not* cache it locally::

        with get_pkg_data_fileobj('standards/vega.fits', cache=False) as fobj:
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
        return get_readable_fileobj(datafn, encoding=encoding)
    else:  # remote file
        return get_readable_fileobj(DATAURL() + datafn, encoding=encoding,
                                    cache=cache)


def get_pkg_data_filename(data_name):
    """
    Retrieves a data file from the standard locations for the package and
    provides a local filename for the data.

    This function is similar to `get_pkg_data_fileobj` but returns the
    file *name* instead of a readable file-like object.  This means
    that this function must always cache remote files locally, unlike
    `get_pkg_data_fileobj`.

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
            return download_file(DATAURL() + data_name, cache=True)
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
            return download_file(DATAURL() + data_name, cache=True)


def get_pkg_data_contents(data_name, encoding=None, cache=True):
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
    get_pkg_data_fileobj : returns a file-like object with the data
    get_pkg_data_filename : returns a local name for a file containing the data
    """
    with get_pkg_data_fileobj(data_name, encoding=encoding, cache=cache) as fd:
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
    from os.path import dirname, join
    from ..utils.misc import find_current_module

    module = find_current_module(1, True)
    if module is None:
        # not called from inside an astropy package.  So just pass name through
        return data_name
    rootpkgname = module.__package__.split('.')[0]

    rootpkg = __import__(rootpkgname)

    module_path = dirname(module.__file__)
    path = join(module_path, data_name)

    root_dir = dirname(rootpkg.__file__)
    assert _is_inside(path, root_dir), \
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
        dldir, urlmapfn = _get_download_cache_locs()
    except (IOError, OSError) as e:
        msg = 'Could not access cache directory to search for data file: '
        warn(CacheMissingWarning(msg + str(e)))
        return None
    hashfn = join(dldir, hash)
    if isfile(hashfn):
        return hashfn
    else:
        return None


def download_file(remote_url, cache=False):
    """
    Accepts a URL, downloads and optionally caches the result
    returning the filename, with a name determined by the file's MD5
    hash. If ``cache=True`` and the file is present in the cache, just
    returns the filename.

    Parameters
    ----------
    remote_url : str
        The URL of the file to download
    cache : bool, optional
        Whether to use the cache
    """

    import hashlib
    from contextlib import closing
    from tempfile import NamedTemporaryFile
    from shutil import move
    from warnings import warn

    from ..utils.console import ProgressBarOrSpinner

    missing_cache = False

    if cache:
        try:
            dldir, urlmapfn = _get_download_cache_locs()
        except (IOError, OSError) as e:
            msg = 'Remote data cache could not be accessed due to '
            estr = '' if len(e.args) < 1 else (': ' + str(e))
            warn(CacheMissingWarning(msg + e.__class__.__name__ + estr))
            cache = False
            missing_cache = True  # indicates that the cache is missing to raise a warning later
    try:
        if cache:
            # We don't need to acquire the lock here, since we are only reading
            with _open_shelve(urlmapfn, True) as url2hash:
                if str(remote_url) in url2hash:
                    return url2hash[str(remote_url)]

        with closing(urllib2.urlopen(remote_url, timeout=REMOTE_TIMEOUT())) as remote:
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

            dlmsg = "Downloading {0}".format(remote_url)
            with ProgressBarOrSpinner(size, dlmsg) as p:
                try:
                    with NamedTemporaryFile(delete=False) as f:
                        bytes_read = 0
                        block = remote.read(DOWNLOAD_CACHE_BLOCK_SIZE())
                        while block:
                            f.write(block)
                            hash.update(block)
                            bytes_read += len(block)
                            p.update(bytes_read)
                            block = remote.read(DOWNLOAD_CACHE_BLOCK_SIZE())
                except:
                    if os.path.exists(f.name):
                        os.remove(f.name)
                    raise

        if cache:
            _acquire_download_cache_lock()
            try:
                with _open_shelve(urlmapfn, True) as url2hash:
                    # We check now to see if another process has
                    # inadvertently written the file underneath us
                    # already
                    if str(remote_url) in url2hash:
                        return url2hash[str(remote_url)]
                    local_path = os.path.join(dldir, hash.hexdigest())
                    move(f.name, local_path)
                    url2hash[str(remote_url)] = local_path
            finally:
                _release_download_cache_lock()
        else:
            local_path = f.name
            if missing_cache:
                msg = ('File downloaded to temporary location due to problem '
                       'with cache directory and will not be cached.')
                warn(CacheMissingWarning(msg, local_path))
            if DELETE_TEMPORARY_DOWNLOADS_AT_EXIT():
                global _tempfilestodel
                _tempfilestodel.append(local_path)
    except urllib2.URLError as e:
        if hasattr(e, 'reason') and hasattr(e.reason, 'errno') and e.reason.errno == 8:
            e.reason.strerror = e.reason.strerror + '. requested URL: ' + remote_url
            e.reason.args = (e.reason.errno, e.reason.strerror)
        raise e

    return local_path


def _do_download_files_in_parallel(args):
    orig_stdout = sys.stdout
    sys.stdout = io.BytesIO()
    try:
        return download_file(*args)
    finally:
        sys.stdout = orig_stdout


def download_files_in_parallel(urls, cache=False):
    """
    Downloads multiple files in parallel from the given URLs.  Blocks until
    all files have downloaded.  The result is a list of local file paths
    corresponding to the given urls.

    Parameters
    ----------
    urls : list of str
        The URLs to retrieve.

    cache : bool, optional
        Whether to use the cache

    Returns
    -------
    paths : list of str
        The local file paths corresponding to the downloaded URLs.
    """
    from .console import ProgressBar

    # Combine duplicate URLs
    combined_urls = list(set(urls))
    combined_paths = ProgressBar.map(
        _do_download_files_in_parallel,
        [(x, cache) for x in combined_urls],
        multiprocess=True)
    paths = []
    for url in urls:
        paths.append(combined_paths[combined_urls.index(url)])
    return paths


# This is used by download_file and _deltemps to determine the files to delete
# when the interpreter exits
_tempfilestodel = []


@atexit.register
def _deltemps():

    global _tempfilestodel

    while len(_tempfilestodel) > 0:
        fn = _tempfilestodel.pop()
        if os.path.isfile(fn):
            os.remove(fn)


def clear_download_cache(hashorurl=None):
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
    from os.path import join, exists
    from shutil import rmtree
    from warnings import warn

    try:
        dldir, urlmapfn = _get_download_cache_locs()
    except (IOError, OSError) as e:
        msg = 'Not clearing data cache - cache inacessable due to '
        estr = '' if len(e.args) < 1 else (': ' + str(e))
        warn(CacheMissingWarning(msg + e.__class__.__name__ + estr))
        return

    _acquire_download_cache_lock()
    try:
        if hashorurl is None:
            if exists(dldir):
                rmtree(dldir)
            if exists(urlmapfn):
                unlink(urlmapfn)
        else:
            with _open_shelve(urlmapfn, True) as url2hash:
                filepath = join(dldir, hashorurl)
                assert _is_inside(filepath, dldir), \
                       ("attempted to use clear_download_cache on a location" +
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
        # the lock will be gone if rmtree was used above, but release otherwise
        if exists(join(_get_download_cache_locs()[0], 'lock')):
            _release_download_cache_lock()


def _get_download_cache_locs():
    """ Finds the path to the data cache directory and makes them if
    they don't exist.

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

    datadir = join(get_cache_dir(), 'download')
    shelveloc = join(get_cache_dir(), 'download_urlmap')

    if not exists(datadir):
        try:
            mkdir(datadir)
        except OSError as e:
            if not exists(datadir):
                raise
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
#the hash shelve, so this should be used for both.
def _acquire_download_cache_lock():
    """
    Uses the lock directory method.  This is good because `mkdir` is
    atomic at the system call level, so it's thread-safe.
    """
    from os.path import join
    from os import mkdir
    from time import sleep

    lockdir = os.path.join(_get_download_cache_locs()[0], 'lock')
    for i in range(DOWNLOAD_CACHE_LOCK_ATTEMPTS()):
        try:
            mkdir(lockdir)
            #write the pid of this process for informational purposes
            with open(join(lockdir, 'pid'), 'w') as f:
                f.write(str(os.getpid()))

        except OSError:
            sleep(1)
        else:
            return
    msg = 'Unable to acquire lock for cache directory ({0} exists)'
    raise RuntimeError(msg.format(lockdir))


def _release_download_cache_lock():
    from os.path import join, exists, isdir
    from os import rmdir

    lockdir = join(_get_download_cache_locs()[0], 'lock')

    if isdir(lockdir):
        #if the pid file is present, be sure to remove it
        pidfn = join(lockdir, 'pid')
        if exists(pidfn):
            os.remove(pidfn)
        rmdir(lockdir)
    else:
        msg = 'Error releasing lock. "{0}" either does not exist or is not ' +\
              'a directory.'
        raise RuntimeError(msg.format(lockdir))
