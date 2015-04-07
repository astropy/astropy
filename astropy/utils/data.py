# Licensed under a 3-clause BSD style license - see LICENSE.rst

""" This module contains helper functions for accessing, downloading, and
caching data files.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from ..extern import six
from ..extern.six.moves import urllib

import atexit
import contextlib
import fnmatch
import hashlib
import os
import io
import shutil
import socket
import sys
import time

from tempfile import NamedTemporaryFile, gettempdir
from warnings import warn

from .. import config as _config
from ..utils.exceptions import AstropyWarning


__all__ = [
    'Conf', 'conf', 'get_readable_fileobj', 'get_file_contents',
    'get_pkg_data_fileobj', 'get_pkg_data_filename',
    'get_pkg_data_contents', 'get_pkg_data_fileobjs',
    'get_pkg_data_filenames', 'compute_hash', 'clear_download_cache',
    'CacheMissingWarning', 'get_free_space_in_dir',
    'check_free_space_in_dir', 'download_file',
    'download_files_in_parallel']


class Conf(_config.ConfigNamespace):
    """
    Configuration parameters for `astropy.utils.data`.
    """

    dataurl = _config.ConfigItem(
        'http://data.astropy.org/',
        'URL for astropy remote data site.')
    remote_timeout = _config.ConfigItem(
        3.,
        'Time to wait for remote data queries (in seconds).',
        aliases=['astropy.coordinates.name_resolve.name_resolve_timeout'])
    compute_hash_block_size = _config.ConfigItem(
        2 ** 16,  # 64K
        'Block size for computing MD5 file hashes.')
    download_block_size = _config.ConfigItem(
        2 ** 16,  # 64K
        'Number of bytes of remote data to download per step.')
    download_cache_lock_attempts = _config.ConfigItem(
        5,
        'Number of times to try to get the lock ' +
        'while accessing the data cache before giving up.')
    delete_temporary_downloads_at_exit = _config.ConfigItem(
        True,
        'If True, temporary download files created when the cache is '
        'inaccessible will be deleted at the end of the python session.')
conf = Conf()


DATAURL = _config.ConfigAlias(
    '0.4', 'DATAURL', 'dataurl')
REMOTE_TIMEOUT = _config.ConfigAlias(
    '0.4', 'REMOTE_TIMEOUT', 'remote_timeout')
COMPUTE_HASH_BLOCK_SIZE = _config.ConfigAlias(
    '0.4', 'COMPUTE_HASH_BLOCK_SIZE', 'compute_hash_block_size')
DOWNLOAD_CACHE_BLOCK_SIZE = _config.ConfigAlias(
    '0.4', 'DOWNLOAD_CACHE_BLOCK_SIZE', 'download_block_size')
DOWNLOAD_CACHE_LOCK_ATTEMPTS = _config.ConfigAlias(
    '0.4', 'DOWNLOAD_CACHE_LOCK_ATTEMPTS', 'download_cache_lock_attempts')
DELETE_TEMPORARY_DOWNLOADS_AT_EXIT = _config.ConfigAlias(
    '0.4', 'DELETE_TEMPORARY_DOWNLOADS_AT_EXIT', 'delete_temporary_downloads_at_exit')


class CacheMissingWarning(AstropyWarning):
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
    url = urllib.parse.urlparse(string)
    # url[0]==url.scheme, but url[0] is py 2.6-compat
    # we can't just check that url[0] is not an empty string, because
    # file paths in windows would return a non-empty scheme (e.g. e:\\
    # returns 'e').
    return url[0].lower() in ['http', 'https', 'ftp', 'sftp', 'ssh', 'file']


def _is_inside(path, parent_path):
    # We have to try realpath too to avoid issues with symlinks, but we leave
    # abspath because some systems like debian have the absolute path (with no
    # symlinks followed) match, but the real directories in different
    # locations, so need to try both cases.
    return os.path.abspath(path).startswith(os.path.abspath(parent_path)) \
        or os.path.realpath(path).startswith(os.path.realpath(parent_path))


@contextlib.contextmanager
def get_readable_fileobj(name_or_obj, encoding=None, cache=False,
                         show_progress=True, remote_timeout=None):
    """
    Given a filename or a readable file-like object, return a context
    manager that yields a readable file-like object.

    This supports passing filenames, URLs, and readable file-like
    objects, any of which can be compressed in gzip, bzip2 or xz
    (the latter on systems with lzma support).

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
        ``read`` method that on Python 2.x returns `bytes` objects and
        on Python 3.x returns `str` (`unicode`) objects, using
        `locale.getpreferredencoding` as an encoding.  This matches
        the default behavior of the built-in `open` when no ``mode``
        argument is provided.

        When ``'binary'``, returns a file-like object where its ``read``
        method returns `bytes` objects.

        When another string, it is the name of an encoding, and the
        file-like object's ``read`` method will return `str` (`unicode`)
        objects, decoded from binary using the given encoding.

    cache : bool, optional
        Whether to cache the contents of remote URLs.

    show_progress : bool, optional
        Whether to display a progress bar if the file is downloaded
        from a remote server.  Default is `True`.

    remote_timeout : float
        Timeout for remote requests in seconds (default is the configurable
        REMOTE_TIMEOUT, which is 3s by default)

    Returns
    -------
    file : readable file-like object
    """

    # close_fds is a list of file handles created by this function
    # that need to be closed.  We don't want to always just close the
    # returned file handle, because it may simply be the file handle
    # passed in.  In that case it is not the responsibility of this
    # function to close it: doing so could result in a "double close"
    # and an "invalid file descriptor" exception.
    close_fds = []
    delete_fds = []

    if remote_timeout is None:
        # use configfile default
        remote_timeout = conf.remote_timeout

    # Get a file object to the content
    if isinstance(name_or_obj, six.string_types):
        if _is_url(name_or_obj):
            name_or_obj = download_file(
                name_or_obj, cache=cache, show_progress=show_progress,
                timeout=remote_timeout)
        if six.PY3:
            fileobj = io.FileIO(name_or_obj, 'r')
        elif six.PY2:
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
            import bz2
        except ImportError:
            for fd in close_fds:
                fd.close()
            raise ValueError(
                ".bz2 format files are not supported since the Python "
                "interpreter does not include the bz2 module")
        try:
            # bz2.BZ2File does not support file objects, only filenames, so we
            # need to write the data to a temporary file
            tmp = NamedTemporaryFile("wb", delete=False)
            tmp.write(fileobj.read())
            tmp.close()
            delete_fds.append(tmp)
            fileobj_new = bz2.BZ2File(tmp.name, mode='rb')
            fileobj_new.read(1)  # need to check that the file is really bzip2
        except (IOError) as e:  # invalid bzip2 file
            fileobj.seek(0)
            fileobj_new.close()
            # raise
        else:
            fileobj_new.seek(0)
            close_fds.append(fileobj_new)
            fileobj = fileobj_new
    elif signature[:3] == b'\xfd7z':  # xz
        try:
            # for Python < 3.3 try backports.lzma; pyliblzma installs as lzma,
            # but does not support TextIOWrapper
            if sys.version_info >= (3,3,0):
                import lzma
                fileobj_new = lzma.LZMAFile(fileobj, mode='rb')
            else:
                from backports import lzma
                from backports.lzma import LZMAFile
                # when called with file object, returns a non-seekable instance
                # need a filename here, too, so have to write the data to a 
                # temporary file
                tmp = NamedTemporaryFile("wb", delete=False)
                tmp.write(fileobj.read())
                tmp.close()
                delete_fds.append(tmp)
                fileobj_new = LZMAFile(tmp.name, mode='rb')
            fileobj_new.read(1)  # need to check that the file is really xz
        except ImportError:
            for fd in close_fds:
                fd.close()
            raise ValueError(
                ".xz format files are not supported since the Python "
                "interpreter does not include the lzma module. "
                "On Python versions < 3.3 consider installing backports.lzma")
        except (IOError, EOFError) as e:  # invalid xz file
            fileobj.seek(0)
            fileobj_new.close()
            # raise
        else:
            fileobj_new.seek(0)
            fileobj = fileobj_new

    # By this point, we have a file, io.FileIO, gzip.GzipFile, bz2.BZ2File
    # or lzma.LZMAFile instance opened in binary mode (that is, read
    # returns bytes).  Now we need to, if requested, wrap it in a
    # io.TextIOWrapper so read will return unicode based on the
    # encoding parameter.

    if six.PY3:
        needs_textio_wrapper = encoding != 'binary'
    elif six.PY2:
        needs_textio_wrapper = encoding != 'binary' and encoding is not None

    if needs_textio_wrapper:
        # A bz2.BZ2File can not be wrapped by a TextIOWrapper,
        # so we decompress it to a temporary file and then
        # return a handle to that.
        try:
            import bz2
        except ImportError:
            pass
        else:
            if isinstance(fileobj, bz2.BZ2File):
                tmp = NamedTemporaryFile("wb", delete=False)
                data = fileobj.read()
                tmp.write(data)
                tmp.close()
                delete_fds.append(tmp)
                if six.PY3:
                    fileobj = io.FileIO(tmp.name, 'r')
                elif six.PY2:
                    fileobj = open(tmp.name, 'rb')
                close_fds.append(fileobj)

        # On Python 2.x, we need to first wrap the regular `file`
        # instance in a `io.FileIO` object before it can be
        # wrapped in a `TextIOWrapper`.  We don't just create an
        # `io.FileIO` object in the first place, because we can't
        # get a raw file descriptor out of it on Python 2.x, which
        # is required for the XML iterparser.
        if six.PY2 and isinstance(fileobj, file):
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
        The content of the file (as requested by ``encoding``).

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
              ``astropy.pkname``, use ``'data/file.dat'`` to get the
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
        ``read`` method that on Python 2.x returns `bytes` objects and
        on Python 3.x returns `str` (`unicode`) objects, using
        `locale.getpreferredencoding` as an encoding.  This matches
        the default behavior of the built-in `open` when no ``mode``
        argument is provided.

        When ``'binary'``, returns a file-like object where its ``read``
        method returns `bytes` objects.

        When another string, it is the name of an encoding, and the
        file-like object's ``read`` method will return `str` (`unicode`)
        objects, decoded from binary using the given encoding.

    cache : bool
        If True, the file will be downloaded and saved locally or the
        already-cached local copy will be accessed. If False, the
        file-like object will directly access the resource (e.g. if a
        remote URL is accessed, an object like that from
        `urllib2.urlopen` on Python 2 or `urllib.request.urlopen` on
        Python 3 is returned).

    Returns
    -------
    fileobj : file-like
        An object with the contents of the data file available via
        ``read`` function.  Can be used as part of a ``with`` statement,
        automatically closing itself after the ``with`` block.

    Raises
    ------
    urllib2.URLError, urllib.error.URLError
        If a remote file cannot be found.
    IOError
        If problems occur writing or reading a local file.

    Examples
    --------

    This will retrieve a data file and its contents for the `astropy.wcs`
    tests::

        from astropy.utils.data import get_pkg_data_fileobj

        with get_pkg_data_fileobj('data/3d_cd.hdr') as fobj:
            fcontents = fobj.read()

    This would download a data file from the astropy data server
    because the ``standards/vega.fits`` file is not present in the
    source distribution.  It will also save the file locally so the
    next time it is accessed it won't need to be downloaded.::

        from astropy.utils.data import get_pkg_data_fileobj

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
        return get_readable_fileobj(conf.dataurl + datafn, encoding=encoding,
                                    cache=cache)


def get_pkg_data_filename(data_name, show_progress=True, remote_timeout=None):
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
              ``astropy.pkname``, use ``'data/file.dat'`` to get the
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

    show_progress : bool, optional
        Whether to display a progress bar if the file is downloaded
        from a remote server.  Default is `True`.

    remote_timeout : float
        Timeout for the requests in seconds (default is the
        configurable `astropy.utils.data.Conf.remote_timeout`, which
        is 3s by default)

    Raises
    ------
    urllib2.URLError, urllib.error.URLError
        If a remote file cannot be found.
    IOError
        If problems occur writing or reading a local file.

    Returns
    -------
    filename : str
        A file path on the local file system corresponding to the data
        requested in ``data_name``.

    Examples
    --------

    This will retrieve the contents of the data file for the `astropy.wcs`
    tests::

        from astropy.utils.data import get_pkg_data_filename

        fn = get_pkg_data_filename('data/3d_cd.hdr')
        with open(fn) as f:
            fcontents = f.read()


    This retrieves a data file by hash either locally or from the astropy data
    server::

        from astropy.utils.data import get_pkg_data_filename

        fn = get_pkg_data_filename('hash/da34a7b07ef153eede67387bf950bb32')
        with open(fn) as f:
            fcontents = f.read()

    See Also
    --------
    get_pkg_data_contents : returns the contents of a file or url as a bytes object
    get_pkg_data_fileobj : returns a file-like object with the data
    """

    data_name = os.path.normpath(data_name)

    if remote_timeout is None:
        # use configfile default
        remote_timeout = conf.remote_timeout

    if data_name.startswith('hash/'):
        # first try looking for a local version if a hash is specified
        hashfn = _find_hash_fn(data_name[5:])
        if hashfn is None:
            return download_file(
                conf.dataurl + data_name, cache=True,
                show_progress=show_progress,
                timeout=remote_timeout)
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
            return download_file(
                conf.dataurl + data_name, cache=True,
                show_progress=show_progress,
                timeout=remote_timeout)


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
              ``astropy.pkname``, use ``'data/file.dat'`` to get the
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
        ``read`` method that on Python 2.x returns `bytes` objects and
        on Python 3.x returns `str` (`unicode`) objects, using
        `locale.getpreferredencoding` as an encoding.  This matches
        the default behavior of the built-in `open` when no ``mode``
        argument is provided.

        When ``'binary'``, returns a file-like object where its ``read``
        method returns `bytes` objects.

        When another string, it is the name of an encoding, and the
        file-like object's ``read`` method will return `str` (`unicode`)
        objects, decoded from binary using the given encoding.

    cache : bool
        If True, the file will be downloaded and saved locally or the
        already-cached local copy will be accessed. If False, the
        file-like object will directly access the resource (e.g. if a
        remote URL is accessed, an object like that from
        `urllib2.urlopen` on Python 2 or `urllib.request.urlopen` on
        Python 3 is returned).

    Returns
    -------
    contents : bytes
        The complete contents of the file as a bytes object.

    Raises
    ------
    urllib2.URLError, urllib.error.URLError
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
              ``astropy.pkname``, use ``'data'`` to get the
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

        from astropy.utils.data import get_pkg_data_filenames

        for fn in get_pkg_data_filename('maps', '*.hdr'):
            with open(fn) as f:
                fcontents = f.read()
    """

    path = _find_pkg_data_path(datadir)
    if os.path.isfile(path):
        raise IOError(
            "Tried to access a data directory that's actually "
            "a package data file")
    elif os.path.isdir(path):
        for filename in os.listdir(path):
            if fnmatch.fnmatch(filename, pattern):
                yield os.path.join(path, filename)
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
              ``astropy.pkname``, use ``'data'`` to get the
              files in ``astropy/pkgname/data``
            * Remote URLs are not currently supported

    pattern : str, optional
        A UNIX-style filename glob pattern to match files.  See the
        `glob` module in the standard library for more information.
        By default, matches all files.

    encoding : str, optional
        When `None` (default), returns a file-like object with a
        ``read`` method that on Python 2.x returns `bytes` objects and
        on Python 3.x returns `str` (`unicode`) objects, using
        `locale.getpreferredencoding` as an encoding.  This matches
        the default behavior of the built-in `open` when no ``mode``
        argument is provided.

        When ``'binary'``, returns a file-like object where its ``read``
        method returns `bytes` objects.

        When another string, it is the name of an encoding, and the
        file-like object's ``read`` method will return `str` (`unicode`)
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

        from astropy.utils.data import get_pkg_data_filenames

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
    use
    e.g. ``get_pkg_data_filename('hash/a725fa6ba642587436612c2df0451956')``,
    but with the hash for your file in place of the hash in the example.

    Parameters
    ----------
    localfn : str
        The path to the file for which the hash should be generated.

    Returns
    -------
    md5hash : str
        The hex digest of the MD5 hash for the contents of the ``localfn``
        file.

    """

    with open(localfn, 'rb') as f:
        h = hashlib.md5()
        block = f.read(conf.compute_hash_block_size)
        while block:
            h.update(block)
            block = f.read(conf.compute_hash_block_size)

    return h.hexdigest()


def _find_pkg_data_path(data_name):
    """
    Look for data in the source-included data directories and return the
    path.
    """

    from ..utils import find_current_module

    module = find_current_module(1, True)
    if module is None:
        # not called from inside an astropy package.  So just pass name through
        return data_name

    if not hasattr(module, '__package__') or not module.__package__:
        # The __package__ attribute may be missing or set to None; see PEP-366,
        # also astropy issue #1256
        if '.' in module.__name__:
            pkgname = module.__name__.rpartition('.')[0]
        else:
            pkgname = module.__name__
    else:
        pkgname = module.__package__

    rootpkgname = pkgname.partition('.')[0]

    rootpkg = __import__(rootpkgname)

    module_path = os.path.dirname(module.__file__)
    path = os.path.join(module_path, data_name)

    root_dir = os.path.dirname(rootpkg.__file__)
    assert _is_inside(path, root_dir), \
           ("attempted to get a local data file outside "
            "of the " + rootpkgname + " tree")

    return path


def _find_hash_fn(hash):
    """
    Looks for a local file by hash - returns file name if found and a valid
    file, otherwise returns None.
    """

    try:
        dldir, urlmapfn = _get_download_cache_locs()
    except (IOError, OSError) as e:
        msg = 'Could not access cache directory to search for data file: '
        warn(CacheMissingWarning(msg + str(e)))
        return None
    hashfn = os.path.join(dldir, hash)
    if os.path.isfile(hashfn):
        return hashfn
    else:
        return None


def get_free_space_in_dir(path):
    """
    Given a path to a directory, returns the amount of free space (in
    bytes) on that filesystem.

    Parameters
    ----------
    path : str
        The path to a directory

    Returns
    -------
    bytes : int
        The amount of free space on the partition that the directory
        is on.
    """

    if sys.platform.startswith('win'):
        import ctypes
        free_bytes = ctypes.c_ulonglong(0)
        retval = ctypes.windll.kernel32.GetDiskFreeSpaceExW(
                ctypes.c_wchar_p(path), None, None, ctypes.pointer(free_bytes))
        if retval == 0:
            raise IOError('Checking free space on %r failed unexpectedly.' %
                          path)
        return free_bytes.value
    else:
        stat = os.statvfs(path)
        return stat.f_bavail * stat.f_frsize


def check_free_space_in_dir(path, size):
    """
    Determines if a given directory has enough space to hold a file of
    a given size.  Raises an IOError if the file would be too large.

    Parameters
    ----------
    path : str
        The path to a directory

    size : int
        A proposed filesize (in bytes)

    Raises
    -------
    IOError : There is not enough room on the filesystem
    """
    from ..utils.console import human_file_size

    space = get_free_space_in_dir(path)
    if space < size:
        raise IOError(
            "Not enough free space in '{0}' "
            "to download a {1} file".format(
                path, human_file_size(size)))


def download_file(remote_url, cache=False, show_progress=True, timeout=None):
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

    show_progress : bool, optional
        Whether to display a progress bar during the download (default
        is `True`)

    timeout : float, optional
        The timeout, in seconds.  Otherwise, use
        `astropy.utils.data.Conf.remote_timeout`.

    Returns
    -------
    local_path : str
        Returns the local path that the file was download to.

    Raises
    ------
    urllib2.URLError, urllib.error.URLError
        Whenever there's a problem getting the remote file.
    """

    from ..utils.console import ProgressBarOrSpinner

    if timeout is None:
        timeout = conf.remote_timeout

    missing_cache = False

    if timeout is None:
        # use configfile default
        timeout = REMOTE_TIMEOUT()

    if cache:
        try:
            dldir, urlmapfn = _get_download_cache_locs()
        except (IOError, OSError) as e:
            msg = 'Remote data cache could not be accessed due to '
            estr = '' if len(e.args) < 1 else (': ' + str(e))
            warn(CacheMissingWarning(msg + e.__class__.__name__ + estr))
            cache = False
            missing_cache = True  # indicates that the cache is missing to raise a warning later

    if six.PY2 and isinstance(remote_url, six.text_type):
        # shelve DBs don't accept unicode strings in Python 2
        url_key = remote_url.encode('utf-8')
    else:
        url_key = remote_url

    try:
        if cache:
            # We don't need to acquire the lock here, since we are only reading
            with _open_shelve(urlmapfn, True) as url2hash:
                if url_key in url2hash:
                    return url2hash[url_key]

        with contextlib.closing(urllib.request.urlopen(
                remote_url, timeout=timeout)) as remote:
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

            if size is not None:
                check_free_space_in_dir(gettempdir(), size)
                if cache:
                    check_free_space_in_dir(dldir, size)

            if show_progress:
                progress_stream = sys.stdout
            else:
                progress_stream = io.StringIO()

            dlmsg = "Downloading {0}".format(remote_url)
            with ProgressBarOrSpinner(size, dlmsg, file=progress_stream) as p:
                with NamedTemporaryFile(delete=False) as f:
                    try:
                        bytes_read = 0
                        block = remote.read(conf.download_block_size)
                        while block:
                            f.write(block)
                            hash.update(block)
                            bytes_read += len(block)
                            p.update(bytes_read)
                            block = remote.read(conf.download_block_size)
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
                    if url_key in url2hash:
                        return url2hash[url_key]
                    local_path = os.path.join(dldir, hash.hexdigest())
                    shutil.move(f.name, local_path)
                    url2hash[url_key] = local_path
            finally:
                _release_download_cache_lock()
        else:
            local_path = f.name
            if missing_cache:
                msg = ('File downloaded to temporary location due to problem '
                       'with cache directory and will not be cached.')
                warn(CacheMissingWarning(msg, local_path))
            if conf.delete_temporary_downloads_at_exit:
                global _tempfilestodel
                _tempfilestodel.append(local_path)
    except urllib.error.URLError as e:
        if hasattr(e, 'reason') and hasattr(e.reason, 'errno') and e.reason.errno == 8:
            e.reason.strerror = e.reason.strerror + '. requested URL: ' + remote_url
            e.reason.args = (e.reason.errno, e.reason.strerror)
        raise e
    except socket.timeout as e:
        # this isn't supposed to happen, but occasionally a socket.timeout gets
        # through.  It's supposed to be caught in `urrlib2` and raised in this
        # way, but for some reason in mysterious circumstances it doesn't. So
        # we'll just re-raise it here instead
        raise urllib.error.URLError(e)

    return local_path


def _do_download_files_in_parallel(args):
    return download_file(*args, show_progress=False)


def download_files_in_parallel(urls, cache=False, show_progress=True,
                               timeout=None):
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

    show_progress : bool, optional
        Whether to display a progress bar during the download (default
        is `True`)

    timeout : float, optional
        Timeout for the requests in seconds (default is the
        configurable `astropy.utils.data.Conf.remote_timeout`).

    Returns
    -------
    paths : list of str
        The local file paths corresponding to the downloaded URLs.
    """
    from .console import ProgressBar

    if timeout is None:
        timeout = conf.remote_timeout

    if show_progress:
        progress = sys.stdout
    else:
        progress = io.BytesIO()

    if timeout is None:
        # use configfile default
        timeout = REMOTE_TIMEOUT()

    # Combine duplicate URLs
    combined_urls = list(set(urls))
    combined_paths = ProgressBar.map(
        _do_download_files_in_parallel,
        [(x, cache) for x in combined_urls],
        file=progress,
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

    if _tempfilestodel is not None:
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
            if os.path.exists(dldir):
                shutil.rmtree(dldir)
            if os.path.exists(urlmapfn):
                os.unlink(urlmapfn)
        else:
            with _open_shelve(urlmapfn, True) as url2hash:
                filepath = os.path.join(dldir, hashorurl)
                assert _is_inside(filepath, dldir), \
                       ("attempted to use clear_download_cache on a path "
                        "outside the data cache directory")

                # shelve DBs don't accept unicode strings as keys in Python 2
                if six.PY2 and isinstance(hashorurl, six.text_type):
                    hash_key = hashorurl.encode('utf-8')
                else:
                    hash_key = hashorurl

                if os.path.exists(filepath):
                    for k, v in list(six.iteritems(url2hash)):
                        if v == filepath:
                            del url2hash[k]
                    os.unlink(filepath)
                elif hash_key in url2hash:
                    filepath = url2hash[hash_key]
                    del url2hash[hash_key]
                    os.unlink(filepath)
                else:
                    msg = 'Could not find file or url {0}'
                    raise OSError(msg.format(hashorurl))
    finally:
        # the lock will be gone if rmtree was used above, but release otherwise
        if os.path.exists(os.path.join(_get_download_cache_locs()[0], 'lock')):
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

    datadir = os.path.join(get_cache_dir(), 'download')
    shelveloc = os.path.join(get_cache_dir(), 'download_urlmap')

    if not os.path.exists(datadir):
        try:
            os.mkdir(datadir)
        except OSError as e:
            if not os.path.exists(datadir):
                raise
    elif not os.path.isdir(datadir):
        msg = 'Data cache directory {0} is not a directory'
        raise IOError(msg.format(datadir))

    if os.path.isdir(shelveloc):
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

    if six.PY2:
        shelf = shelve.open(shelffn, protocol=2)
    elif six.PY3:
        shelf = shelve.open(shelffn + '.db', protocol=2)

    if withclosing:
        return contextlib.closing(shelf)
    else:
        return shelf


#the cache directory must be locked before any writes are performed.  Same for
#the hash shelve, so this should be used for both.
def _acquire_download_cache_lock():
    """
    Uses the lock directory method.  This is good because `mkdir` is
    atomic at the system call level, so it's thread-safe.
    """

    lockdir = os.path.join(_get_download_cache_locs()[0], 'lock')
    for i in range(conf.download_cache_lock_attempts):
        try:
            os.mkdir(lockdir)
            #write the pid of this process for informational purposes
            with open(os.path.join(lockdir, 'pid'), 'w') as f:
                f.write(str(os.getpid()))

        except OSError:
            time.sleep(1)
        else:
            return
    msg = 'Unable to acquire lock for cache directory ({0} exists)'
    raise RuntimeError(msg.format(lockdir))


def _release_download_cache_lock():
    lockdir = os.path.join(_get_download_cache_locs()[0], 'lock')

    if os.path.isdir(lockdir):
        #if the pid file is present, be sure to remove it
        pidfn = os.path.join(lockdir, 'pid')
        if os.path.exists(pidfn):
            os.remove(pidfn)
        os.rmdir(lockdir)
    else:
        msg = 'Error releasing lock. "{0}" either does not exist or is not ' +\
              'a directory.'
        raise RuntimeError(msg.format(lockdir))
