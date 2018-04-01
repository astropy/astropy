# Licensed under a 3-clause BSD style license - see LICENSE.rst

""" This module contains helper functions for accessing, downloading, and
caching data files.
"""

import atexit
import contextlib
import fnmatch
import hashlib
import os
import io
import pathlib
import shutil
import socket
import sys
import time
import urllib.request
import urllib.error
import urllib.parse
import shelve

from tempfile import NamedTemporaryFile, gettempdir
from warnings import warn

from .. import config as _config
from ..utils.exceptions import AstropyWarning
from ..utils.introspection import find_current_module, resolve_name

__all__ = [
    'Conf', 'conf', 'get_readable_fileobj', 'get_file_contents',
    'get_pkg_data_fileobj', 'get_pkg_data_filename',
    'get_pkg_data_contents', 'get_pkg_data_fileobjs',
    'get_pkg_data_filenames', 'compute_hash', 'clear_download_cache',
    'CacheMissingWarning', 'get_free_space_in_dir',
    'check_free_space_in_dir', 'download_file',
    'download_files_in_parallel', 'is_url_in_cache', 'get_cached_urls']

_dataurls_to_alias = {}


class Conf(_config.ConfigNamespace):
    """
    Configuration parameters for `astropy.utils.data`.
    """

    dataurl = _config.ConfigItem(
        'http://data.astropy.org/',
        'Primary URL for astropy remote data site.')
    dataurl_mirror = _config.ConfigItem(
        'http://www.astropy.org/astropy-data/',
        'Mirror URL for astropy remote data site.')
    remote_timeout = _config.ConfigItem(
        10.,
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
    # we can't just check that url.scheme is not an empty string, because
    # file paths in windows would return a non-empty scheme (e.g. e:\\
    # returns 'e').
    return url.scheme.lower() in ['http', 'https', 'ftp', 'sftp', 'ssh', 'file']


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
    Given a filename, pathlib.Path object or a readable file-like object, return a context
    manager that yields a readable file-like object.

    This supports passing filenames, URLs, and readable file-like objects,
    any of which can be compressed in gzip, bzip2 or lzma (xz) if the
    appropriate compression libraries are provided by the Python installation.

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
        ``read`` method that returns `str` (``unicode``) objects, using
        `locale.getpreferredencoding` as an encoding.  This matches
        the default behavior of the built-in `open` when no ``mode``
        argument is provided.

        When ``'binary'``, returns a file-like object where its ``read``
        method returns `bytes` objects.

        When another string, it is the name of an encoding, and the
        file-like object's ``read`` method will return `str` (``unicode``)
        objects, decoded from binary using the given encoding.

    cache : bool, optional
        Whether to cache the contents of remote URLs.

    show_progress : bool, optional
        Whether to display a progress bar if the file is downloaded
        from a remote server.  Default is `True`.

    remote_timeout : float
        Timeout for remote requests in seconds (default is the configurable
        `astropy.utils.data.Conf.remote_timeout`, which is 3s by default)

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
    PATH_TYPES = (str, pathlib.Path)

    close_fds = []
    delete_fds = []

    if remote_timeout is None:
        # use configfile default
        remote_timeout = conf.remote_timeout

    # Get a file object to the content
    if isinstance(name_or_obj, PATH_TYPES):
        # name_or_obj could be a Path object if pathlib is available
        name_or_obj = str(name_or_obj)

        is_url = _is_url(name_or_obj)
        if is_url:
            name_or_obj = download_file(
                name_or_obj, cache=cache, show_progress=show_progress,
                timeout=remote_timeout)
        fileobj = io.FileIO(name_or_obj, 'r')
        if is_url and not cache:
            delete_fds.append(fileobj)
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
            import gzip
            fileobj_new = gzip.GzipFile(fileobj=fileobj, mode='rb')
            fileobj_new.read(1)  # need to check that the file is really gzip
        except (OSError, EOFError, struct.error):  # invalid gzip file
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
            with NamedTemporaryFile("wb", delete=False) as tmp:
                tmp.write(fileobj.read())
                tmp.close()
                fileobj_new = bz2.BZ2File(tmp.name, mode='rb')
            fileobj_new.read(1)  # need to check that the file is really bzip2
        except OSError:  # invalid bzip2 file
            fileobj.seek(0)
            fileobj_new.close()
            # raise
        else:
            fileobj_new.seek(0)
            close_fds.append(fileobj_new)
            fileobj = fileobj_new
    elif signature[:3] == b'\xfd7z':  # xz
        try:
            import lzma
            fileobj_new = lzma.LZMAFile(fileobj, mode='rb')
            fileobj_new.read(1)  # need to check that the file is really xz
        except ImportError:
            for fd in close_fds:
                fd.close()
            raise ValueError(
                ".xz format files are not supported since the Python "
                "interpreter does not include the lzma module.")
        except (OSError, EOFError) as e:  # invalid xz file
            fileobj.seek(0)
            fileobj_new.close()
            # should we propagate this to the caller to signal bad content?
            # raise ValueError(e)
        else:
            fileobj_new.seek(0)
            fileobj = fileobj_new

    # By this point, we have a file, io.FileIO, gzip.GzipFile, bz2.BZ2File
    # or lzma.LZMAFile instance opened in binary mode (that is, read
    # returns bytes).  Now we need to, if requested, wrap it in a
    # io.TextIOWrapper so read will return unicode based on the
    # encoding parameter.

    needs_textio_wrapper = encoding != 'binary'

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

                fileobj = io.FileIO(tmp.name, 'r')
                close_fds.append(fileobj)

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


def get_file_contents(*args, **kwargs):
    """
    Retrieves the contents of a filename or file-like object.

    See  the `get_readable_fileobj` docstring for details on parameters.

    Returns
    -------
    content
        The content of the file (as requested by ``encoding``).

    """
    with get_readable_fileobj(*args, **kwargs) as f:
        return f.read()


@contextlib.contextmanager
def get_pkg_data_fileobj(data_name, package=None, encoding=None, cache=True):
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
              e.g. 'hash/34c33b3eb0d56eb9462003af249eff28'.  The hash
              will first be searched for locally, and if not found,
              the Astropy data server will be queried.

    package : str, optional
        If specified, look for a file relative to the given package, rather
        than the default of looking relative to the calling module's package.

    encoding : str, optional
        When `None` (default), returns a file-like object with a
        ``read`` method returns `str` (``unicode``) objects, using
        `locale.getpreferredencoding` as an encoding.  This matches
        the default behavior of the built-in `open` when no ``mode``
        argument is provided.

        When ``'binary'``, returns a file-like object where its ``read``
        method returns `bytes` objects.

        When another string, it is the name of an encoding, and the
        file-like object's ``read`` method will return `str` (``unicode``)
        objects, decoded from binary using the given encoding.

    cache : bool
        If True, the file will be downloaded and saved locally or the
        already-cached local copy will be accessed. If False, the
        file-like object will directly access the resource (e.g. if a
        remote URL is accessed, an object like that from
        `urllib.request.urlopen` is returned).

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
    OSError
        If problems occur writing or reading a local file.

    Examples
    --------

    This will retrieve a data file and its contents for the `astropy.wcs`
    tests::

        >>> from astropy.utils.data import get_pkg_data_fileobj
        >>> with get_pkg_data_fileobj('data/3d_cd.hdr',
        ...                           package='astropy.wcs.tests') as fobj:
        ...     fcontents = fobj.read()
        ...

    This next example would download a data file from the astropy data server
    because the ``allsky/allsky_rosat.fits`` file is not present in the
    source distribution.  It will also save the file locally so the
    next time it is accessed it won't need to be downloaded.::

        >>> from astropy.utils.data import get_pkg_data_fileobj
        >>> with get_pkg_data_fileobj('allsky/allsky_rosat.fits',
        ...                           encoding='binary') as fobj:  # doctest: +REMOTE_DATA
        ...     fcontents = fobj.read()
        ...
        Downloading http://data.astropy.org/allsky/allsky_rosat.fits [Done]

    This does the same thing but does *not* cache it locally::

        >>> with get_pkg_data_fileobj('allsky/allsky_rosat.fits',
        ...                           encoding='binary', cache=False) as fobj:  # doctest: +REMOTE_DATA
        ...     fcontents = fobj.read()
        ...
        Downloading http://data.astropy.org/allsky/allsky_rosat.fits [Done]

    See Also
    --------
    get_pkg_data_contents : returns the contents of a file or url as a bytes object
    get_pkg_data_filename : returns a local name for a file containing the data
    """

    datafn = _find_pkg_data_path(data_name, package=package)
    if os.path.isdir(datafn):
        raise OSError("Tried to access a data file that's actually "
                      "a package data directory")
    elif os.path.isfile(datafn):  # local file
        with get_readable_fileobj(datafn, encoding=encoding) as fileobj:
            yield fileobj
    else:  # remote file
        all_urls = (conf.dataurl, conf.dataurl_mirror)
        for url in all_urls:
            try:
                with get_readable_fileobj(url + data_name, encoding=encoding,
                                          cache=cache) as fileobj:
                    # We read a byte to trigger any URLErrors
                    fileobj.read(1)
                    fileobj.seek(0)
                    yield fileobj
                    break
            except urllib.error.URLError:
                pass
        else:
            urls = '\n'.join('  - {0}'.format(url) for url in all_urls)
            raise urllib.error.URLError("Failed to download {0} from the following "
                                        "repositories:\n\n{1}".format(data_name, urls))


def get_pkg_data_filename(data_name, package=None, show_progress=True,
                          remote_timeout=None):
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
              e.g. 'hash/34c33b3eb0d56eb9462003af249eff28'.  The hash
              will first be searched for locally, and if not found,
              the Astropy data server will be queried.

    package : str, optional
        If specified, look for a file relative to the given package, rather
        than the default of looking relative to the calling module's package.

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
    OSError
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

        >>> from astropy.utils.data import get_pkg_data_filename
        >>> fn = get_pkg_data_filename('data/3d_cd.hdr',
        ...                            package='astropy.wcs.tests')
        >>> with open(fn) as f:
        ...     fcontents = f.read()
        ...

    This retrieves a data file by hash either locally or from the astropy data
    server::

        >>> from astropy.utils.data import get_pkg_data_filename
        >>> fn = get_pkg_data_filename('hash/34c33b3eb0d56eb9462003af249eff28')  # doctest: +SKIP
        >>> with open(fn) as f:
        ...     fcontents = f.read()
        ...

    See Also
    --------
    get_pkg_data_contents : returns the contents of a file or url as a bytes object
    get_pkg_data_fileobj : returns a file-like object with the data
    """

    if remote_timeout is None:
        # use configfile default
        remote_timeout = conf.remote_timeout

    if data_name.startswith('hash/'):
        # first try looking for a local version if a hash is specified
        hashfn = _find_hash_fn(data_name[5:])

        if hashfn is None:
            all_urls = (conf.dataurl, conf.dataurl_mirror)
            for url in all_urls:
                try:
                    return download_file(url + data_name, cache=True,
                                         show_progress=show_progress,
                                         timeout=remote_timeout)
                except urllib.error.URLError:
                    pass
            urls = '\n'.join('  - {0}'.format(url) for url in all_urls)
            raise urllib.error.URLError("Failed to download {0} from the following "
                                        "repositories:\n\n{1}\n\n".format(data_name, urls))

        else:
            return hashfn
    else:
        fs_path = os.path.normpath(data_name)
        datafn = _find_pkg_data_path(fs_path, package=package)
        if os.path.isdir(datafn):
            raise OSError("Tried to access a data file that's actually "
                          "a package data directory")
        elif os.path.isfile(datafn):  # local file
            return datafn
        else:  # remote file
            all_urls = (conf.dataurl, conf.dataurl_mirror)
            for url in all_urls:
                try:
                    return download_file(url + data_name, cache=True,
                                         show_progress=show_progress,
                                         timeout=remote_timeout)
                except urllib.error.URLError:
                    pass
            urls = '\n'.join('  - {0}'.format(url) for url in all_urls)
            raise urllib.error.URLError("Failed to download {0} from the following "
                                        "repositories:\n\n{1}".format(data_name, urls))


def get_pkg_data_contents(data_name, package=None, encoding=None, cache=True):
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
              e.g. 'hash/34c33b3eb0d56eb9462003af249eff28'.  The hash
              will first be searched for locally, and if not found,
              the Astropy data server will be queried.
            * A URL to some other file.

    package : str, optional
        If specified, look for a file relative to the given package, rather
        than the default of looking relative to the calling module's package.


    encoding : str, optional
        When `None` (default), returns a file-like object with a
        ``read`` method that returns `str` (``unicode``) objects, using
        `locale.getpreferredencoding` as an encoding.  This matches
        the default behavior of the built-in `open` when no ``mode``
        argument is provided.

        When ``'binary'``, returns a file-like object where its ``read``
        method returns `bytes` objects.

        When another string, it is the name of an encoding, and the
        file-like object's ``read`` method will return `str` (``unicode``)
        objects, decoded from binary using the given encoding.

    cache : bool
        If True, the file will be downloaded and saved locally or the
        already-cached local copy will be accessed. If False, the
        file-like object will directly access the resource (e.g. if a
        remote URL is accessed, an object like that from
        `urllib.request.urlopen` is returned).

    Returns
    -------
    contents : bytes
        The complete contents of the file as a bytes object.

    Raises
    ------
    urllib2.URLError, urllib.error.URLError
        If a remote file cannot be found.
    OSError
        If problems occur writing or reading a local file.

    See Also
    --------
    get_pkg_data_fileobj : returns a file-like object with the data
    get_pkg_data_filename : returns a local name for a file containing the data
    """

    with get_pkg_data_fileobj(data_name, package=package, encoding=encoding,
                              cache=cache) as fd:
        contents = fd.read()
    return contents


def get_pkg_data_filenames(datadir, package=None, pattern='*'):
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

    package : str, optional
        If specified, look for a file relative to the given package, rather
        than the default of looking relative to the calling module's package.

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

        >>> from astropy.utils.data import get_pkg_data_filenames
        >>> for fn in get_pkg_data_filenames('maps', 'astropy.wcs.tests',
        ...                                  '*.hdr'):
        ...     with open(fn) as f:
        ...         fcontents = f.read()
        ...
    """

    path = _find_pkg_data_path(datadir, package=package)
    if os.path.isfile(path):
        raise OSError(
            "Tried to access a data directory that's actually "
            "a package data file")
    elif os.path.isdir(path):
        for filename in os.listdir(path):
            if fnmatch.fnmatch(filename, pattern):
                yield os.path.join(path, filename)
    else:
        raise OSError("Path not found")


def get_pkg_data_fileobjs(datadir, package=None, pattern='*', encoding=None):
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

    package : str, optional
        If specified, look for a file relative to the given package, rather
        than the default of looking relative to the calling module's package.

    pattern : str, optional
        A UNIX-style filename glob pattern to match files.  See the
        `glob` module in the standard library for more information.
        By default, matches all files.

    encoding : str, optional
        When `None` (default), returns a file-like object with a
        ``read`` method that returns `str` (``unicode``) objects, using
        `locale.getpreferredencoding` as an encoding.  This matches
        the default behavior of the built-in `open` when no ``mode``
        argument is provided.

        When ``'binary'``, returns a file-like object where its ``read``
        method returns `bytes` objects.

        When another string, it is the name of an encoding, and the
        file-like object's ``read`` method will return `str` (``unicode``)
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

        >>> from astropy.utils.data import get_pkg_data_filenames
        >>> for fd in get_pkg_data_fileobjs('maps', 'astropy.wcs.tests',
        ...                                 '*.hdr'):
        ...     fcontents = fd.read()
        ...
    """

    for fn in get_pkg_data_filenames(datadir, package=package,
                                     pattern=pattern):
        with get_readable_fileobj(fn, encoding=encoding) as fd:
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
    e.g. ``get_pkg_data_filename('hash/34c33b3eb0d56eb9462003af249eff28')``,
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


def _find_pkg_data_path(data_name, package=None):
    """
    Look for data in the source-included data directories and return the
    path.
    """

    if package is None:
        module = find_current_module(1, finddiff=['astropy.utils.data', 'contextlib'])
        if module is None:
            # not called from inside an astropy package.  So just pass name
            # through
            return data_name

        if not hasattr(module, '__package__') or not module.__package__:
            # The __package__ attribute may be missing or set to None; see
            # PEP-366, also astropy issue #1256
            if '.' in module.__name__:
                package = module.__name__.rpartition('.')[0]
            else:
                package = module.__name__
        else:
            package = module.__package__
    else:
        module = resolve_name(package)

    rootpkgname = package.partition('.')[0]

    rootpkg = resolve_name(rootpkgname)

    module_path = os.path.dirname(module.__file__)
    path = os.path.join(module_path, data_name)

    root_dir = os.path.dirname(rootpkg.__file__)
    if not _is_inside(path, root_dir):
        raise RuntimeError("attempted to get a local data file outside "
                           "of the {} tree.".format(rootpkgname))

    return path


def _find_hash_fn(hash):
    """
    Looks for a local file by hash - returns file name if found and a valid
    file, otherwise returns None.
    """

    try:
        dldir, urlmapfn = _get_download_cache_locs()
    except OSError as e:
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
            raise OSError('Checking free space on {!r} failed '
                          'unexpectedly.'.format(path))
        return free_bytes.value
    else:
        stat = os.statvfs(path)
        return stat.f_bavail * stat.f_frsize


def check_free_space_in_dir(path, size):
    """
    Determines if a given directory has enough space to hold a file of
    a given size.  Raises an OSError if the file would be too large.

    Parameters
    ----------
    path : str
        The path to a directory

    size : int
        A proposed filesize (in bytes)

    Raises
    -------
    OSError : There is not enough room on the filesystem
    """
    from ..utils.console import human_file_size

    space = get_free_space_in_dir(path)
    if space < size:
        raise OSError(
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

    if cache:
        try:
            dldir, urlmapfn = _get_download_cache_locs()
        except OSError as e:
            msg = 'Remote data cache could not be accessed due to '
            estr = '' if len(e.args) < 1 else (': ' + str(e))
            warn(CacheMissingWarning(msg + e.__class__.__name__ + estr))
            cache = False
            missing_cache = True  # indicates that the cache is missing to raise a warning later

    url_key = remote_url

    # Check if URL is Astropy data server, which has alias, and cache it.
    if (url_key.startswith(conf.dataurl) and
            conf.dataurl not in _dataurls_to_alias):
        with urllib.request.urlopen(conf.dataurl, timeout=timeout) as remote:
            _dataurls_to_alias[conf.dataurl] = [conf.dataurl, remote.geturl()]

    try:
        if cache:
            # We don't need to acquire the lock here, since we are only reading
            with shelve.open(urlmapfn) as url2hash:
                if url_key in url2hash:
                    return url2hash[url_key]
                # If there is a cached copy from mirror, use it.
                else:
                    for cur_url in _dataurls_to_alias.get(conf.dataurl, []):
                        if url_key.startswith(cur_url):
                            url_mirror = url_key.replace(cur_url,
                                                         conf.dataurl_mirror)
                            if url_mirror in url2hash:
                                return url2hash[url_mirror]

        with urllib.request.urlopen(remote_url, timeout=timeout) as remote:
            # keep a hash to rename the local file to the hashed name
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
                    except BaseException:
                        if os.path.exists(f.name):
                            os.remove(f.name)
                        raise

        if cache:
            _acquire_download_cache_lock()
            try:
                with shelve.open(urlmapfn) as url2hash:
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


def is_url_in_cache(url_key):
    """
    Check if a download from ``url_key`` is in the cache.

    Parameters
    ----------
    url_key : string
        The URL retrieved

    Returns
    -------
    in_cache : bool
        `True` if a download from ``url_key`` is in the cache
    """
    # The code below is modified from astropy.utils.data.download_file()
    try:
        dldir, urlmapfn = _get_download_cache_locs()
    except OSError as e:
        msg = 'Remote data cache could not be accessed due to '
        estr = '' if len(e.args) < 1 else (': ' + str(e))
        warn(CacheMissingWarning(msg + e.__class__.__name__ + estr))
        return False

    with shelve.open(urlmapfn) as url2hash:
        if url_key in url2hash:
            return True
    return False


def _do_download_files_in_parallel(args):
    return download_file(*args)


def download_files_in_parallel(urls, cache=True, show_progress=True,
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
        Whether to use the cache (default is `True`).

        .. versionchanged:: 3.0
            The default was changed to ``True`` and setting it to ``False`` will
            print a Warning and set it to ``True`` again, because the function
            will not work properly without cache.

    show_progress : bool, optional
        Whether to display a progress bar during the download (default
        is `True`)

    timeout : float, optional
        Timeout for each individual requests in seconds (default is the
        configurable `astropy.utils.data.Conf.remote_timeout`).

    Returns
    -------
    paths : list of str
        The local file paths corresponding to the downloaded URLs.
    """
    from .console import ProgressBar

    if timeout is None:
        timeout = conf.remote_timeout

    if not cache:
        # See issue #6662, on windows won't work because the files are removed
        # again before they can be used. On *NIX systems it will behave as if
        # cache was set to True because multiprocessing cannot insert the items
        # in the list of to-be-removed files.
        warn("Disabling the cache does not work because of multiprocessing, it "
             "will be set to ``True``. You may need to manually remove the "
             "cached files afterwards.", AstropyWarning)
        cache = True

    if show_progress:
        progress = sys.stdout
    else:
        progress = io.BytesIO()

    # Combine duplicate URLs
    combined_urls = list(set(urls))
    combined_paths = ProgressBar.map(
        _do_download_files_in_parallel,
        [(x, cache, False, timeout) for x in combined_urls],
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
        should be removed from the cache if present.
    """

    try:
        dldir, urlmapfn = _get_download_cache_locs()
    except OSError as e:
        msg = 'Not clearing data cache - cache inacessable due to '
        estr = '' if len(e.args) < 1 else (': ' + str(e))
        warn(CacheMissingWarning(msg + e.__class__.__name__ + estr))
        return

    _acquire_download_cache_lock()
    try:
        if hashorurl is None:
            # dldir includes both the download files and the urlmapfn.  This structure
            # is required since we cannot know a priori the actual file name corresponding
            # to the shelve map named urlmapfn.
            if os.path.exists(dldir):
                shutil.rmtree(dldir)
        else:
            with shelve.open(urlmapfn) as url2hash:
                filepath = os.path.join(dldir, hashorurl)
                if not _is_inside(filepath, dldir):
                    raise RuntimeError("attempted to use clear_download_cache on"
                                       " a path outside the data cache directory")

                hash_key = hashorurl

                if os.path.exists(filepath):
                    for k, v in url2hash.items():
                        if v == filepath:
                            del url2hash[k]
                    os.unlink(filepath)
                elif hash_key in url2hash:
                    filepath = url2hash[hash_key]
                    del url2hash[hash_key]
                    if os.path.exists(filepath):
                        # Make sure the filepath still actually exists (perhaps user removed it)
                        os.unlink(filepath)
                # Otherwise could not find file or url, but no worries.
                # Clearing download cache just makes sure that the file or url
                # is no longer in the cache regardless of starting condition.

    finally:
        # the lock will be gone if rmtree was used above, but release otherwise
        if os.path.exists(os.path.join(dldir, 'lock')):
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

    # datadir includes both the download files and the shelveloc.  This structure
    # is required since we cannot know a priori the actual file name corresponding
    # to the shelve map named shelveloc.  (The backend can vary and is allowed to
    # do whatever it wants with the filename.  Filename munging can and does happen
    # in practice).
    py_version = 'py' + str(sys.version_info.major)
    datadir = os.path.join(get_cache_dir(), 'download', py_version)
    shelveloc = os.path.join(datadir, 'urlmap')

    if not os.path.exists(datadir):
        try:
            os.makedirs(datadir)
        except OSError as e:
            if not os.path.exists(datadir):
                raise
    elif not os.path.isdir(datadir):
        msg = 'Data cache directory {0} is not a directory'
        raise OSError(msg.format(datadir))

    if os.path.isdir(shelveloc):
        msg = 'Data cache shelve object location {0} is a directory'
        raise OSError(msg.format(shelveloc))

    return datadir, shelveloc


# the cache directory must be locked before any writes are performed.  Same for
# the hash shelve, so this should be used for both.
def _acquire_download_cache_lock():
    """
    Uses the lock directory method.  This is good because `mkdir` is
    atomic at the system call level, so it's thread-safe.
    """

    lockdir = os.path.join(_get_download_cache_locs()[0], 'lock')
    for i in range(conf.download_cache_lock_attempts):
        try:
            os.mkdir(lockdir)
            # write the pid of this process for informational purposes
            with open(os.path.join(lockdir, 'pid'), 'w') as f:
                f.write(str(os.getpid()))

        except OSError:
            time.sleep(1)
        else:
            return
    msg = ("Unable to acquire lock for cache directory ({0} exists). "
           "You may need to delete the lock if the python interpreter wasn't "
           "shut down properly.")
    raise RuntimeError(msg.format(lockdir))


def _release_download_cache_lock():
    lockdir = os.path.join(_get_download_cache_locs()[0], 'lock')

    if os.path.isdir(lockdir):
        # if the pid file is present, be sure to remove it
        pidfn = os.path.join(lockdir, 'pid')
        if os.path.exists(pidfn):
            os.remove(pidfn)
        os.rmdir(lockdir)
    else:
        msg = 'Error releasing lock. "{0}" either does not exist or is not ' +\
              'a directory.'
        raise RuntimeError(msg.format(lockdir))


def get_cached_urls():
    """
    Get the list of URLs in the cache. Especially useful for looking up what
    files are stored in your cache when you don't have internet access.

    Returns
    -------
    cached_urls : list
        List of cached URLs.
    """
    # The code below is modified from astropy.utils.data.download_file()
    try:
        dldir, urlmapfn = _get_download_cache_locs()
    except OSError as e:
        msg = 'Remote data cache could not be accessed due to '
        estr = '' if len(e.args) < 1 else (': ' + str(e))
        warn(CacheMissingWarning(msg + e.__class__.__name__ + estr))
        return False

    with shelve.open(urlmapfn) as url2hash:
        return list(url2hash.keys())
