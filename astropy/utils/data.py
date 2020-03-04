# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Functions for accessing, downloading, and caching data files."""

import atexit
import contextlib
import dbm
import errno
import fnmatch
import hashlib
import os
import io
import pathlib
import re
import shutil
import socket
import sys
import time
import urllib.request
import urllib.error
import urllib.parse
import shelve
import zipfile

from tempfile import NamedTemporaryFile, gettempdir, TemporaryDirectory
from warnings import warn

import astropy.config.paths
from astropy import config as _config
from astropy.utils.exceptions import AstropyWarning
from astropy.utils.introspection import find_current_module, resolve_name

# Order here determines order in the autosummary
__all__ = [
    'Conf', 'conf',
    'download_file', 'download_files_in_parallel',
    'get_readable_fileobj',
    'get_pkg_data_fileobj', 'get_pkg_data_filename',
    'get_pkg_data_contents', 'get_pkg_data_fileobjs',
    'get_pkg_data_filenames',
    'is_url_in_cache', 'get_cached_urls',
    'cache_total_size', 'cache_contents',
    'export_download_cache', 'import_download_cache', 'import_file_to_cache',
    'check_download_cache',
    'clear_download_cache',
    'compute_hash',
    'get_free_space_in_dir',
    'check_free_space_in_dir',
    'get_file_contents',
    'CacheMissingWarning',
]

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
    default_http_user_agent = _config.ConfigItem(
        'astropy',
        'Default User-Agent for HTTP request headers. This can be overwritten'
        'for a particular call via http_headers option, where available.'
        'This only provides the default value when not set by https_headers.')
    remote_timeout = _config.ConfigItem(
        10.,
        'Time to wait for remote data queries (in seconds).',
        aliases=['astropy.coordinates.name_resolve.name_resolve_timeout'])
    compute_hash_block_size = _config.ConfigItem(
        2 ** 16,  # 64K
        'Block size for computing file hashes.')
    download_block_size = _config.ConfigItem(
        2 ** 16,  # 64K
        'Number of bytes of remote data to download per step.')
    download_cache_lock_attempts = _config.ConfigItem(
        5,
        'Number of seconds to wait for the cache lock to be free. It should '
        'normally only ever be held long enough to copy an already-downloaded '
        'file into the cache, so this will normally only run over if '
        'something goes wrong and the lock is left held by a dead process; '
        'the exception raised should indicate this and what to do to fix it.')
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
                         show_progress=True, remote_timeout=None,
                         sources=None, http_headers=None):
    """Yield a readable, seekable file-like object from a file or URL.

    This supports passing filenames, URLs, and readable file-like objects,
    any of which can be compressed in gzip, bzip2 or lzma (xz) if the
    appropriate compression libraries are provided by the Python installation.

    Notes
    -----

    This function is a context manager, and should be used for example
    as::

        with get_readable_fileobj('file.dat') as f:
            contents = f.read()

    If a URL is provided and the cache is in use, the provided URL will be the
    name used in the cache. The contents may already be stored in the cache
    under this URL provided, they may be downloaded from this URL, or they may
    be downloaded from one of the locations listed in ``sources``. See
    `~download_file` for details.

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

    cache : bool or "update", optional
        Whether to cache the contents of remote URLs. If "update",
        check the remote URL for a new version but store the result
        in the cache.

    show_progress : bool, optional
        Whether to display a progress bar if the file is downloaded
        from a remote server.  Default is `True`.

    remote_timeout : float
        Timeout for remote requests in seconds (default is the configurable
        `astropy.utils.data.Conf.remote_timeout`, which is 3s by default)

    sources : list of str, optional
        If provided, a list of URLs to try to obtain the file from. The
        result will be stored under the original URL. The original URL
        will *not* be tried unless it is in this list; this is to prevent
        long waits for a primary server that is known to be inaccessible
        at the moment.

    http_headers : dict or None
        HTTP request headers to pass into ``urlopen`` if needed. (These headers
        are ignored if the protocol for the ``name_or_obj``/``sources`` entry
        is not a remote HTTP URL.) In the default case (None), the headers are
        ``User-Agent: some_value`` and ``Accept: */*``, where ``some_value``
        is set by ``astropy.utils.data.conf.default_http_user_agent``.

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
                timeout=remote_timeout, sources=sources,
                http_headers=http_headers)
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
        try:
            # py.path.LocalPath objects have .read() method but it uses
            # text mode, which won't work. .read_binary() does, and
            # surely other ducks would return binary contents when
            # called like this.
            # py.path.LocalPath is what comes from the tmpdir fixture
            # in pytest.
            fileobj = io.BytesIO(fileobj.read_binary())
        except AttributeError:
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
        except (OSError, EOFError):  # invalid xz file
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
    urllib.error.URLError
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
        ...                           encoding='binary') as fobj:  # doctest: +REMOTE_DATA +IGNORE_OUTPUT
        ...     fcontents = fobj.read()
        ...
        Downloading http://data.astropy.org/allsky/allsky_rosat.fits [Done]

    This does the same thing but does *not* cache it locally::

        >>> with get_pkg_data_fileobj('allsky/allsky_rosat.fits',
        ...                           encoding='binary', cache=False) as fobj:  # doctest: +REMOTE_DATA +IGNORE_OUTPUT
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
        with get_readable_fileobj(
            conf.dataurl + data_name,
            encoding=encoding,
            cache=cache,
            sources=[conf.dataurl + data_name,
                     conf.dataurl_mirror + data_name],
        ) as fileobj:
            # We read a byte to trigger any URLErrors
            fileobj.read(1)
            fileobj.seek(0)
            yield fileobj


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
    urllib.error.URLError
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
            return download_file(conf.dataurl + data_name, cache=True,
                                 show_progress=show_progress,
                                 timeout=remote_timeout,
                                 sources=[conf.dataurl + data_name,
                                          conf.dataurl_mirror + data_name])
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
            return download_file(conf.dataurl + data_name, cache=True,
                                 show_progress=show_progress,
                                 timeout=remote_timeout,
                                 sources=[conf.dataurl + data_name,
                                          conf.dataurl_mirror + data_name])


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
    urllib.error.URLError
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
        >>> for fn in get_pkg_data_filenames('data/maps', 'astropy.wcs.tests',
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
        >>> for fd in get_pkg_data_fileobjs('data/maps', 'astropy.wcs.tests',
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
    hash : str
        The hex digest of the cryptographic hash for the contents of the
        ``localfn`` file.
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
        raise RuntimeError(f"attempted to get a local data file outside "
                           f"of the {rootpkgname} tree.")

    return path


def _find_hash_fn(hexdigest, pkgname='astropy'):
    """
    Looks for a local file by hash - returns file name if found and a valid
    file, otherwise returns None.
    """

    with _cache(pkgname) as (dldir, url2hash):
        if dldir is None:
            return None
        hashfn = os.path.join(dldir, hexdigest)
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
    if not os.path.isdir(path):
        raise OSError(
            "Can only determine free space associated with directories, "
            "not files.")
        # Actually you can on Linux but I want to avoid code that fails
        # on Windows only.
    return shutil.disk_usage(path).free


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
    from astropy.utils.console import human_file_size

    space = get_free_space_in_dir(path)
    if space < size:
        raise OSError(f"Not enough free space in {path} "
                      f"to download a {human_file_size(size)} file")


def _download_file_from_source(source_url, show_progress=True, timeout=None,
                               remote_url=None, cache=False, pkgname='astropy',
                               http_headers=None):
    from astropy.utils.console import ProgressBarOrSpinner

    if remote_url is None:
        remote_url = source_url
    if http_headers is None:
        http_headers = {}

    req = urllib.request.Request(source_url, headers=http_headers)
    with urllib.request.urlopen(req, timeout=timeout) as remote:
        # keep a hash to rename the local file to the hashed name
        hasher = hashlib.md5()

        info = remote.info()
        try:
            size = int(info['Content-Length'])
        except (KeyError, ValueError, TypeError):
            size = None

        if size is not None:
            check_free_space_in_dir(gettempdir(), size)
            if cache:
                with _cache(pkgname) as (dldir, url2hash):
                    check_free_space_in_dir(dldir, size)

        if show_progress and sys.stdout.isatty():
            progress_stream = sys.stdout
        else:
            progress_stream = io.StringIO()

        if source_url == remote_url:
            dlmsg = f"Downloading {remote_url}"
        else:
            dlmsg = f"Downloading {remote_url} from {source_url}"
        with ProgressBarOrSpinner(size, dlmsg, file=progress_stream) as p:
            with NamedTemporaryFile(prefix=f"astropy-download-{os.getpid()}-",
                                    delete=False) as f:
                try:
                    bytes_read = 0
                    block = remote.read(conf.download_block_size)
                    while block:
                        f.write(block)
                        hasher.update(block)
                        bytes_read += len(block)
                        p.update(bytes_read)
                        block = remote.read(conf.download_block_size)
                        if size is not None and bytes_read > size:
                            raise urllib.error.URLError(
                                f"File was supposed to be {size} bytes but "
                                f"server provides more, at least {bytes_read} "
                                f"bytes. Download failed.")
                    if size is not None and bytes_read < size:
                        raise urllib.error.ContentTooShortError(
                            f"File was supposed to be {size} bytes but we "
                            f"only got {bytes_read} bytes. Download failed.",
                            content=None)
                except BaseException:
                    if os.path.exists(f.name):
                        try:
                            os.remove(f.name)
                        except OSError:
                            pass
                    raise
    return f.name, hasher.hexdigest()


def download_file(remote_url, cache=False, show_progress=True, timeout=None,
                  sources=None, pkgname='astropy', http_headers=None):
    """Downloads a URL and optionally caches the result.

    It returns the filename of a file containing the URL's contents.
    If ``cache=True`` and the file is present in the cache, just
    returns the filename; if the file had to be downloaded, add it
    to the cache.

    The cache is effectively a dictionary mapping URLs to files; by default the
    file contains the contents of the URL that is its key, but in practice
    these can be obtained from a mirror (using ``sources``) or imported from
    the local filesystem (using `~import_file_to_cache` or
    `~import_download_cache`).  Regardless, each file is regarded as
    representing the contents of a particular URL, and this URL should be used
    to look them up or otherwise manipulate them.

    The files in the cache directory are named according to a cryptographic
    hash of their contents (currently MD5, so hackers can cause collisions).
    Thus files with the same content share storage. The modification times on
    these files normally indicate when they were last downloaded from the
    Internet.

    Parameters
    ----------
    remote_url : str
        The URL of the file to download

    cache : bool or "update", optional
        Whether to cache the contents of remote URLs. If "update",
        always download the remote URL in case there is a new version
        and store the result in the cache.

    show_progress : bool, optional
        Whether to display a progress bar during the download (default
        is `True`). Regardless of this setting, the progress bar is only
        displayed when outputting to a terminal.

    timeout : float, optional
        The timeout, in seconds.  Otherwise, use
        `astropy.utils.data.Conf.remote_timeout`.

    sources : list of str, optional
        If provided, a list of URLs to try to obtain the file from. The
        result will be stored under the original URL. The original URL
        will *not* be tried unless it is in this list; this is to prevent
        long waits for a primary server that is known to be inaccessible
        at the moment. If an empty list is passed, then ``download_file``
        will not attempt to connect to the Internet.

    pkgname : `str`, optional
        The package name to use to locate the download cache. i.e. for
        ``pkgname='astropy'`` the default cache location is
        ``~/.astropy/cache``.

    http_headers : dict or None
        HTTP request headers to pass into ``urlopen`` if needed. (These headers
        are ignored if the protocol for the ``name_or_obj``/``sources`` entry
        is not a remote HTTP URL.) In the default case (None), the headers are
        ``User-Agent: some_value`` and ``Accept: */*``, where ``some_value``
        is set by ``astropy.utils.data.conf.default_http_user_agent``.

    Returns
    -------
    local_path : str
        Returns the local path that the file was download to.

    Raises
    ------
    urllib.error.URLError
        Whenever there's a problem getting the remote file.

    Notes
    -----
    Because this returns a filename, another process could run
    clear_download_cache before you actually open the file, leaving
    you with a filename that no longer points to a usable file.
    """
    if timeout is None:
        timeout = conf.remote_timeout
    if sources is None:
        sources = [remote_url]
    if http_headers is None:
        http_headers = {'User-Agent': conf.default_http_user_agent,
                        'Accept': '*/*'}

    missing_cache = ""

    url_key = remote_url

    if cache:
        with _cache(pkgname) as (dldir, url2hash):
            if dldir is None:
                cache = False
                missing_cache = (
                    "Cache directory cannot be read or created, "
                    "providing data in temporary file instead."
                )
            elif cache != "update" and url_key in url2hash:
                return url2hash[url_key]

    errors = {}
    for source_url in sources:
        try:
            f_name, hexdigest = _download_file_from_source(
                    source_url,
                    timeout=timeout,
                    show_progress=show_progress,
                    cache=cache,
                    remote_url=remote_url,
                    pkgname=pkgname,
                    http_headers=http_headers)
            # Success!
            break

        except urllib.error.URLError as e:
            # errno 8 is from SSL "EOF occurred in violation of protocol"
            if (hasattr(e, 'reason')
                    and hasattr(e.reason, 'errno')
                    and e.reason.errno == 8):
                e.reason.strerror = (e.reason.strerror +
                                     '. requested URL: '
                                     + remote_url)
                e.reason.args = (e.reason.errno, e.reason.strerror)
            errors[source_url] = e
        except socket.timeout as e:
            # this isn't supposed to happen, but occasionally a socket.timeout
            # gets through.  It's supposed to be caught in urllib and raised
            # in this way, but for some reason in mysterious circumstances it
            # doesn't (or didn't in python2?). So we'll just re-raise it here
            # instead.
            errors[source_url] = e
    else:   # No success
        if not sources:
            raise KeyError(
                f"No sources listed and file {remote_url} not in cache! "
                f"Please include primary URL in sources if you want it to be "
                f"included as a valid source.")
        elif len(sources) == 1:
            raise errors[sources[0]]
        else:
            raise urllib.error.URLError(
                f"Unable to open any source! Exceptions were {errors}") \
                from errors[sources[0]]

    if cache:
        try:
            return import_file_to_cache(url_key, f_name,
                                        hexdigest=hexdigest,
                                        remove_original=True,
                                        pkgname=pkgname)
        except WrongDBMModule as e:
            missing_cache = (
                f"{e}; Unable to use cache, providing data in temporary file "
                f"{f_name} instead.")
        except PermissionError:
            # Cache is readonly, we can't update it
            missing_cache = (
                f"Cache directory appears to be read-only, unable to import "
                f"downloaded file, providing data in temporary file {f_name} "
                f"instead.")

    if missing_cache:
        warn(CacheMissingWarning(missing_cache, f_name))
    if conf.delete_temporary_downloads_at_exit:
        global _tempfilestodel
        _tempfilestodel.append(f_name)
    return f_name


def is_url_in_cache(url_key, pkgname='astropy'):
    """Check if a download for ``url_key`` is in the cache.

    The provided ``url_key`` will be the name used in the cache. The contents
    may have been downloaded from this URL or from a mirror or they may have
    been provided by the user. See `~download_file` for details.

    Parameters
    ----------
    url_key : str
        The URL retrieved
    pkgname : `str`, optional
        The package name to use to locate the download cache. i.e. for
        ``pkgname='astropy'`` the default cache location is
        ``~/.astropy/cache``.


    Returns
    -------
    in_cache : bool
        `True` if a download for ``url_key`` is in the cache, `False` if not
        or if the cache does not exist at all.

    See Also
    --------
    cache_contents : obtain a dictionary listing everything in the cache
    """
    with _cache(pkgname) as (dldir, url2hash):
        return url_key in url2hash


def cache_total_size(pkgname='astropy'):
    """Return the total size in bytes of all files in the cache."""
    with _cache(pkgname) as (dldir, url2hash):
        return sum(os.path.getsize(os.path.join(dldir, h))
                   for h in url2hash.values())


def _do_download_files_in_parallel(kwargs):
    with astropy.config.paths.set_temp_config(kwargs.pop("temp_config")):
        with astropy.config.paths.set_temp_cache(kwargs.pop("temp_cache")):
            return download_file(**kwargs)


def download_files_in_parallel(urls,
                               cache="update",
                               show_progress=True,
                               timeout=None,
                               sources=None,
                               multiprocessing_start_method=None,
                               pkgname='astropy'):
    """Download multiple files in parallel from the given URLs.

    Blocks until all files have downloaded.  The result is a list of
    local file paths corresponding to the given urls.

    The results will be stored in the cache under the values in ``urls`` even
    if they are obtained from some other location via ``sources``. See
    `~download_file` for details.

    Parameters
    ----------
    urls : list of str
        The URLs to retrieve.

    cache : bool or "update", optional
        Whether to use the cache (default is `True`). If "update",
        always download the remote URLs to see if new data is available
        and store the result in cache.

        .. versionchanged:: 4.0
            The default was changed to ``"update"`` and setting it to
            ``False`` will print a Warning and set it to ``"update"`` again,
            because the function will not work properly without cache. Using
            ``True`` will work as expected.

        .. versionchanged:: 3.0
            The default was changed to ``True`` and setting it to ``False``
            will print a Warning and set it to ``True`` again, because the
            function will not work properly without cache.

    show_progress : bool, optional
        Whether to display a progress bar during the download (default
        is `True`)

    timeout : float, optional
        Timeout for each individual requests in seconds (default is the
        configurable `astropy.utils.data.Conf.remote_timeout`).

    sources : dict, optional
        If provided, for each URL a list of URLs to try to obtain the
        file from. The result will be stored under the original URL.
        For any URL in this dictionary, the original URL will *not* be
        tried unless it is in this list; this is to prevent long waits
        for a primary server that is known to be inaccessible at the
        moment.

    multiprocessing_start_method : str, optional
        Useful primarily for testing; if in doubt leave it as the default.
        When using multiprocessing, certain anomalies occur when starting
        processes with the "spawn" method (the only option on Windows);
        other anomalies occur with the "fork" method (the default on
        Linux).

    pkgname : `str`, optional
        The package name to use to locate the download cache. i.e. for
        ``pkgname='astropy'`` the default cache location is
        ``~/.astropy/cache``.

    Returns
    -------
    paths : list of str
        The local file paths corresponding to the downloaded URLs.

    Notes
    -----
    If a URL is unreachable, the downloading will grind to a halt and the
    exception will propagate upward, but an unpredictable number of
    files will have been successfully downloaded and will remain in
    the cache.
    """
    from .console import ProgressBar

    if timeout is None:
        timeout = conf.remote_timeout
    if sources is None:
        sources = {}

    if not cache:
        # See issue #6662, on windows won't work because the files are removed
        # again before they can be used. On *NIX systems it will behave as if
        # cache was set to True because multiprocessing cannot insert the items
        # in the list of to-be-removed files. This could be fixed, but really,
        # just use the cache, with update_cache if appropriate.
        warn('Disabling the cache does not work because of multiprocessing, '
             'it will be set to ``"update"``. You may need to manually remove '
             'the cached files with clear_download_cache() afterwards.',
             AstropyWarning)
        cache = "update"

    if show_progress:
        progress = sys.stdout
    else:
        progress = io.BytesIO()

    # Combine duplicate URLs
    combined_urls = list(set(urls))
    combined_paths = ProgressBar.map(
        _do_download_files_in_parallel,
        [dict(remote_url=u,
              cache=cache,
              show_progress=False,
              timeout=timeout,
              sources=sources.get(u, None),
              pkgname=pkgname,
              temp_cache=astropy.config.paths.set_temp_cache._temp_path,
              temp_config=astropy.config.paths.set_temp_config._temp_path)
         for u in combined_urls],
        file=progress,
        multiprocess=True,
        multiprocessing_start_method=multiprocessing_start_method,
    )
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
                try:
                    os.remove(fn)
                except OSError:
                    # oh well we tried
                    # could be held open by some process, on Windows
                    pass


def clear_download_cache(hashorurl=None, pkgname='astropy'):
    """Clears the data file cache by deleting the local file(s).

    If a URL is provided, it will be the name used in the cache. The contents
    may have been downloaded from this URL or from a mirror or they may have
    been provided by the user. See `~download_file` for details.

    For the purposes of this function, a file can also be identified by a hash
    of its contents or by the filename under which the data is stored (as
    returned by `~download_file`, for example).

    Parameters
    ----------
    hashorurl : str or None
        If None, the whole cache is cleared.  Otherwise, specify
        a hash for the cached file that is supposed to be deleted,
        the full path to a file in the cache that should be deleted,
        or a URL that should be removed from the cache if present.

    pkgname : `str`, optional
        The package name to use to locate the download cache. i.e. for
        ``pkgname='astropy'`` the default cache location is
        ``~/.astropy/cache``.
    """
    zapped_cache = False
    try:
        with contextlib.ExitStack() as stack:
            try:
                dldir, url2hash = stack.enter_context(
                    _cache(pkgname, write=True))
            except (RuntimeError, WrongDBMModule):  # Couldn't get lock
                if hashorurl is None:
                    # Release lock by blowing away cache
                    # Need to get locations
                    dldir, _ = _get_download_cache_locs(pkgname)
                else:
                    # Can't do specific deletion without the lock
                    raise
            except OSError as e:
                # Problem arose when trying to open the cache
                msg = 'Not clearing data cache - cache inaccessible due to '
                estr = '' if len(e.args) < 1 else (': ' + str(e))
                warn(CacheMissingWarning(msg + e.__class__.__name__ + estr))
                return
            if hashorurl is None:
                if os.path.exists(dldir):
                    shutil.rmtree(dldir)
                    zapped_cache = True
            elif _is_url(hashorurl):
                try:
                    filepath = url2hash.pop(hashorurl)
                    if not any(v == filepath for v in url2hash.values()):
                        try:
                            os.unlink(filepath)
                        except FileNotFoundError:
                            # Maybe someone else got it first
                            pass
                    return
                except KeyError:
                    pass
            else:  # it's a path
                filepath = os.path.join(dldir, hashorurl)
                if not _is_inside(filepath, dldir):
                    # Should this be ValueError? IOError?
                    raise RuntimeError(
                        f"attempted to use clear_download_cache on the path "
                        f"{filepath} outside the data cache directory {dldir}")
                if os.path.exists(filepath):
                    # Convert to list because we'll be modifying it as we go
                    for k, v in list(url2hash.items()):
                        if v == filepath:
                            del url2hash[k]
                    os.unlink(filepath)
                    return
                # Otherwise could not find file or url, but no worries.
                # Clearing download cache just makes sure that the file or url
                # is no longer in the cache regardless of starting condition.
    except OSError as e:
        if zapped_cache and e.errno == errno.ENOENT:
            # We just deleted the directory and, on Windows (?) the "dumb"
            # backend tried to write itself out to a nonexistent directory.
            # It's fine for this to fail.
            return
        else:
            msg = 'Not clearing data from cache - problem arose '
            estr = '' if len(e.args) < 1 else (': ' + str(e))
            warn(CacheMissingWarning(msg + e.__class__.__name__ + estr))
            return


def _get_download_cache_locs(pkgname='astropy'):
    """Finds the path to the cache directory and makes them if they don't exist.

    Parameters
    ----------
    pkgname : `str`, optional
        The package name to use to locate the download cache. i.e. for
        ``pkgname='astropy'`` the default cache location is
        ``~/.astropy/cache``.

    Returns
    -------
    datadir : str
        The path to the data cache directory.
    shelveloc : str
        The path to the shelve object that stores the cache info.
    """
    from astropy.config.paths import get_cache_dir

    # datadir includes both the download files and the shelveloc.  This structure
    # is required since we cannot know a priori the actual file name corresponding
    # to the shelve map named shelveloc.  (The backend can vary and is allowed to
    # do whatever it wants with the filename.  Filename munging can and does happen
    # in practice).
    py_version = 'py' + str(sys.version_info.major)
    datadir = os.path.join(get_cache_dir(pkgname), 'download', py_version)
    shelveloc = os.path.join(datadir, 'urlmap')

    if not os.path.exists(datadir):
        try:
            os.makedirs(datadir)
        except OSError:
            if not os.path.exists(datadir):
                raise
    elif not os.path.isdir(datadir):
        raise OSError(f'Data cache directory {datadir} is not a directory')

    if os.path.isdir(shelveloc):
        raise OSError(
            f'Data cache shelve object location {shelveloc} is a directory')

    return datadir, shelveloc


def _keep_trying(timeout):
    waited = 0.
    yield waited
    while waited < timeout:
        pid_based_random = (hash(str(os.getpid()))
                            / 2.**sys.hash_info.width)
        dt = 0.05*(1+pid_based_random)
        time.sleep(dt)
        waited += dt
        yield waited


# the cache directory must be locked before any writes are performed.  Same for
# the hash shelve, so this should be used for both.
#
# In fact the shelve, if you're using gdbm (the default on Linux), it can't
# be open for reading if someone has it open for writing. So we need to lock
# access to the shelve even for reading.
@contextlib.contextmanager
def _cache_lock(pkgname, need_write=False):
    lockdir = os.path.join(_get_download_cache_locs(pkgname)[0], 'lock')
    pidfn = os.path.join(lockdir, 'pid')
    got_lock = False
    try:
        msg = f"Config file requests {conf.download_cache_lock_attempts} tries"
        for waited in _keep_trying(conf.download_cache_lock_attempts):
            try:
                os.mkdir(lockdir)
            except FileExistsError:
                # It is not safe to open and inspect the pid file here
                # on Windows, because while it is open that prevents its
                # deletion and that of the directory containing it. This
                # gets in the way of useful error messages.
                msg = (
                    f"Cache is locked after {waited:.2f} s. This may indicate "
                    f"an astropy bug or that kill -9 was used. If you want to "
                    f"unlock the cache remove the directory {lockdir}.")
            except OSError as e:
                # PermissionError doesn't cover all read-only-ness, just EACCES
                if e.errno in [errno.EPERM,    # Operation not permitted
                               errno.EACCES,   # Permission denied
                               errno.EROFS]:   # File system is read-only
                    if need_write:
                        raise
                    else:
                        break
                else:
                    raise
            else:
                got_lock = True
                # write the pid of this process for informational purposes
                with open(pidfn, 'w') as f:
                    f.write(str(os.getpid()))
                break
        else:
            # Never did get the lock.
            try:
                # Might as well try to read and report the PID file, at
                # this point it's unlikely we'll block someone else trying
                # to exit cleanly.
                pid = get_file_contents(pidfn)
            except OSError:
                pass
            else:
                msg += f" Lock claims to be held by process {pid}."
            raise RuntimeError(msg)

        yield

    finally:
        # clear_download_cache might have deleted the lockdir
        if got_lock and os.path.exists(lockdir):
            if os.path.isdir(lockdir):
                # if the pid file is present, be sure to remove it
                if os.path.exists(pidfn):
                    os.remove(pidfn)
                os.rmdir(lockdir)
            else:
                raise RuntimeError(
                    f'Error releasing lock. {lockdir} exists but is not '
                    f'a directory.')
        else:
            # Just in case we were called from _clear_download_cache
            # or something went wrong before creating the directory
            # or the cache was readonly; no need to clean it up then.
            pass


class ReadOnlyDict(dict):
    def __setitem__(self, key, value):
        raise TypeError("This object is read-only.")


_NOTHING = ReadOnlyDict()  # might as well share.


class WrongDBMModule(dbm.error[0]):
    pass


class WrongDBMModuleWarning(CacheMissingWarning):
    """
    This warning indicates the standard cache directory is not accessible,
    specifically because it exists but is in a format that the current python
    interpreter cannot understand.
    """


@contextlib.contextmanager
def _cache(pkgname, write=False):
    """Download cache context manager.

    Yields a pair consisting of the download cache directory and a dict-like
    object that maps downloaded URLs to the filenames that contain their
    contents; these files reside in the download cache directory.

    If writing is requested, this holds the lock and yields a modifiable
    dict-like object (actually a shelve object from the shelve module).  If
    there is something wrong with the cache setup, an appropriate exception
    will be raised.

    If reading is requested, the lock will be briefly acquired, the URL map
    will be copied into a read-only dict-like object which is yielded, and the
    lock will be released. If some problem occurs and the cache is inaccessible
    or non-functional, a CacheMissingWarning will be emitted and this context
    manager will yield an empty dict-like object and None as the download
    directory.

    Although this cache lives behind a lock, and files are not normally
    modified, it is possible to break concurrent access - most easily by using
    clear_download_cache on a URL while someone has the filename and wants to
    read the file. Since download_file returns a filename, there is not much we
    can do about this.  get_readable_fileobj doesn't quite avoid the problem,
    though it almost immediately opens the filename, preserving the contents
    from deletion.
    """
    try:
        dldir, urlmapfn = _get_download_cache_locs(pkgname)
    except OSError as e:
        if write:
            raise
        else:
            msg = 'Remote data cache could not be accessed due to '
            estr = '' if len(e.args) < 1 else (': ' + str(e))
            warn(CacheMissingWarning(msg + e.__class__.__name__ + estr))
            yield None, _NOTHING
            return
    wrong_dbm_message = (
        "Existing astropy cache is in an unsupported format, "
        "either install the appropriate package or use "
        "astropy.utils.data.clear_download_cache() to delete the "
        "whole cache; ")
    if write:
        try:
            with _cache_lock(pkgname, need_write=True), \
                    shelve.open(urlmapfn, flag="c") as url2hash:
                yield dldir, url2hash
        except dbm.error as e:
            if "module is not available" in str(e):
                raise WrongDBMModule(wrong_dbm_message + str(e))
            else:
                raise
    else:
        try:
            with _cache_lock(pkgname), shelve.open(urlmapfn, flag="r") as url2hash:
                # Copy so we can release the lock.
                d = ReadOnlyDict(url2hash.items())
        except dbm.error as e:
            # Might be a "file not found" - that is, an un-initialized cache,
            # might be something serious, no way to tell as shelve just gives
            # you a plain dbm.error
            # Also the file doesn't have a platform-independent name, so good
            # luck diagnosing the problem if it is one.
            if "module is not available" in str(e):
                warn(WrongDBMModuleWarning(wrong_dbm_message + str(e)))
            d = _NOTHING
        yield dldir, d


class CacheDamaged(ValueError):
    """Record the URL or file that was a problem.

    Using clear_download_cache on the .bad_file or .bad_url attribute,
    whichever is not None, should resolve this particular problem.
    """
    def __init__(self, *args, bad_urls=None, bad_files=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.bad_urls = bad_urls if bad_urls is not None else []
        self.bad_files = bad_files if bad_files is not None else []


def check_download_cache(check_hashes=False, pkgname='astropy'):
    """Do a consistency check on the cache.

    Because the cache is shared by all versions of astropy in all virtualenvs
    run by your user, possibly concurrently, it could accumulate problems.
    This could lead to hard-to-debug problems or wasted space. This function
    detects a number of incorrect conditions, including nonexistent files that
    are indexed, files that are indexed but in the wrong place, and, if you
    request it, files whose content does not match the hash that is indexed.

    This function also returns a list of non-indexed files. A few will be
    associated with the shelve object; their exact names depend on the backend
    used but will probably be based on ``urlmap``. The presence of other files
    probably indicates that something has gone wrong and inaccessible files
    have accumulated in the cache. These can be removed with
    `clear_download_cache`, either passing the filename returned here, or
    with no arguments to empty the entire cache and return it to a
    reasonable, if empty, state.

    Parameters
    ----------
    check_hashes : bool, optional
        Whether to compute the hashes of the contents of all files in the
        cache and compare them to the names. This can take some time if
        the cache contains large files.
    pkgname : `str`, optional
        The package name to use to locate the download cache. i.e. for
        ``pkgname='astropy'`` the default cache location is
        ``~/.astropy/cache``.

    Returns
    -------
    strays : set of strings
        This is the set of files in the cache directory that do not correspond
        to known URLs. This may include some files associated with the cache
        index.

    Raises
    ------
    CacheDamaged
        To indicate a problem with the cache contents; the exception contains
        a .bad_url or a .bad_file attribute to allow the user to use
        `clear_download_cache` to remove the offending item.
    OSError, RuntimeError
        To indicate some problem with the cache structure. This may need a full
        `clear_download_cache` to resolve, or may indicate some kind of
        misconfiguration.
    """
    with _cache(pkgname) as (dldir, url2hash):
        if dldir is None:
            raise OSError("Cache directory cannot be created or accessed")
        nonexistent_targets = {}
        bad_locations = {}
        bad_hashes = {}
        abandoned_files = set()
        leftover_files = set(os.path.join(dldir, k)
                             for k in os.listdir(dldir))
        leftover_files.discard(os.path.join(dldir, "urlmap"))
        leftover_files.discard(os.path.join(dldir, "lock"))
        for u, h in url2hash.items():
            if not os.path.exists(h):
                nonexistent_targets[u] = h
            leftover_files.discard(h)
            d, hexdigest = os.path.split(h)
            if dldir != d:
                bad_locations[u] = h
            if check_hashes:
                hexdigest_file = compute_hash(h)
                if hexdigest_file != hexdigest:
                    bad_hashes[u] = h
        for h in leftover_files:
            h_base = os.path.basename(h)
            if len(h_base) >= 32 and re.match(r"[0-9a-f]+", h_base):
                abandoned_files.add(h)
        leftover_files -= abandoned_files

        msgs = []
        if nonexistent_targets:
            msgs.append(
                f"URL(s) point(s) to nonexistent file(s): {nonexistent_targets}")
        if bad_locations:
            msgs.append(
                f"URL(s) point(s) to file(s) outside {dldir}: {bad_locations}")
        if bad_hashes:
            msgs.append(
                f"Filename(s) does not match hash(es) of contents: {bad_hashes}")
        if abandoned_files:
            msgs.append(
                f"Apparently abandoned files: {abandoned_files}")
        if msgs:
            raise CacheDamaged(
                "\n".join(msgs),
                bad_urls=(list(nonexistent_targets.keys())
                          + list(bad_locations.keys())),
                bad_files=(list(bad_hashes.values())
                           + list(abandoned_files)))
        else:
            return leftover_files


def import_file_to_cache(url_key, filename,
                         hexdigest=None,
                         remove_original=False,
                         pkgname='astropy'):
    """Import the on-disk file specified by filename to the cache.

    The provided ``url_key`` will be the name used in the cache. The file
    should contain the contents of this URL, at least notionally (the URL may
    be temporarily or permanently unavailable). It is using ``url_key`` that
    users will request these contents from the cache. See `~download_file` for
    details.

    If ``url_key`` already exists in the cache, it will be updated to point to
    these imported contents, and its old contents will be deleted from the
    cache if nothing else points there.

    If a file already exists in the cache with the same contents (no matter the
    ``url_key``), its modification time will be updated but this file will not
    be copied/moved over it.

    Parameters
    ----------
    url_key : str
        The key to index the file under. This should probably be
        the URL where the file was located, though if you obtained
        it from a mirror you should use the URL of the primary
        location.
    filename : str
        The file whose contents you want to import.
    hexdigest : str, optional
        The cryptographic hash of the file, as computed by
        `compute_hash`. If it is not provided, `compute_hash`
        will be called. This argument is available in case
        it is easier to compute the hash progressively as
        the file is downloaded.
    remove_original : bool
        Whether to remove the original file (``filename``) once import is
        complete. If the original is to be removed, and if the cache lives
        on the same filesystem as the file to be imported, a simple rename
        operation will move the file into the cache. Otherwise the file
        will be copied and then the original deleted if appropriate.
    pkgname : `str`, optional
        The package name to use to locate the download cache. i.e. for
        ``pkgname='astropy'`` the default cache location is
        ``~/.astropy/cache``.
    """
    if hexdigest is None:
        hexdigest = compute_hash(filename, pkgname=pkgname)
    with _cache(pkgname, write=True) as (dldir, url2hash):
        # We check now to see if another process has
        # inadvertently written the file underneath us
        # already
        local_path = os.path.join(dldir, hexdigest)
        if os.path.exists(local_path):
            # Same hash, no problem
            if remove_original:
                os.remove(filename)
            # Update file modification time.
            try:
                with open(local_path, "ab"):
                    pass
            except OSError:
                pass
        elif remove_original:
            shutil.move(filename, local_path)
        else:
            shutil.copy(filename, local_path)
        old_hash = url2hash.get(url_key, None)
        url2hash[url_key] = local_path
        if old_hash is not None:
            if old_hash not in url2hash.values():
                try:
                    os.remove(os.path.join(dldir, old_hash))
                except OSError as e:
                    warn(f"Unable to remove no-longer-referenced previous "
                         f"contents of {url_key} because of: {e}",
                         AstropyWarning)
        return url2hash[url_key]


def get_cached_urls(pkgname='astropy'):
    """
    Get the list of URLs in the cache. Especially useful for looking up what
    files are stored in your cache when you don't have internet access.

    The listed URLs are the keys programs should use to access the file
    contents, but those contents may have actually been obtained from a mirror.
    See `~download_file` for details.

    Parameters
    ----------
    pkgname : `str`, optional
        The package name to use to locate the download cache. i.e. for
        ``pkgname='astropy'`` the default cache location is
        ``~/.astropy/cache``.

    Returns
    -------
    cached_urls : list
        List of cached URLs.

    See Also
    --------
    cache_contents : obtain a dictionary listing everything in the cache
    """
    with _cache(pkgname) as (dldir, url2hash):
        return list(url2hash.keys())


def cache_contents(pkgname='astropy'):
    """Obtain a dict mapping cached URLs to filenames.

    This dictionary is a read-only snapshot of the state of the cache when this
    function was called. If other processes are actively working with the
    cache, it is possible for them to delete files that are listed in this
    dictionary. Use with some caution if you are working on a system that is
    busy with many running astropy processes, although the same issues apply to
    most functions in this module.
    """
    with _cache(pkgname) as (dldir, url2hash):
        return ReadOnlyDict((k, os.path.join(dldir, v))
                            for (k, v) in url2hash.items())


def export_download_cache(filename_or_obj, urls=None, overwrite=False, pkgname='astropy'):
    """Exports the cache contents as a ZIP file.

    Parameters
    ----------
    filename_or_obj : str or file-like
        Where to put the created ZIP file. Must be something the zipfile
        module can write to.
    urls : iterable of str or None
        The URLs to include in the exported cache. The default is all
        URLs currently in the cache. If a URL is included in this list
        but is not currently in the cache, a KeyError will be raised.
        To ensure that all are in the cache use `~download_file`
        or `~download_files_in_parallel`.
    overwrite : bool, optional
        If filename_or_obj is a filename that exists, it will only be
        overwritten if this is True.
    pkgname : `str`, optional
        The package name to use to locate the download cache. i.e. for
        ``pkgname='astropy'`` the default cache location is
        ``~/.astropy/cache``.

    See Also
    --------
    import_download_cache : import the contents of such a ZIP file
    import_file_to_cache : import a single file directly
    """
    if urls is None:
        urls = get_cached_urls(pkgname)
    with zipfile.ZipFile(filename_or_obj, 'w' if overwrite else 'x') as z:
        for u in urls:
            fn = download_file(u, cache=True, sources=[], pkgname=pkgname)
            # Do not use os.path.join because ZIP files want
            # "/" on all platforms
            z_fn = urllib.parse.quote(u, safe="")
            z.write(fn, z_fn)


def import_download_cache(filename_or_obj, urls=None, update_cache=False, pkgname='astropy'):
    """Imports the contents of a ZIP file into the cache.

    Each member of the ZIP file should be named by a quoted version of the
    URL whose contents it stores. These names are decoded with
    :func:`~urllib.parse.unquote`.

    Parameters
    ----------
    filename_or_obj : str or file-like
        Where the stored ZIP file is. Must be something the :mod:`~zipfile`
        module can read from.
    urls : set of str or list of str or None
        The URLs to import from the ZIP file. The default is all
        URLs in the file.
    update_cache : bool, optional
        If True, any entry in the ZIP file will overwrite the value in the
        cache; if False, leave untouched any entry already in the cache.
    pkgname : `str`, optional
        The package name to use to locate the download cache. i.e. for
        ``pkgname='astropy'`` the default cache location is
        ``~/.astropy/cache``.

    See Also
    --------
    export_download_cache : export the contents the cache to of such a ZIP file
    import_file_to_cache : import a single file directly
    """
    with zipfile.ZipFile(filename_or_obj, 'r') as z, TemporaryDirectory() as d:
        for i, zf in enumerate(z.infolist()):
            url = urllib.parse.unquote(zf.filename)
            # FIXME(aarchiba): do we want some kind of validation on this URL?
            # urllib.parse might do something sensible...but what URLs might
            # they have?
            # is_url in this file is probably a good check, not just here
            # but throughout this file.
            if urls is not None and url not in urls:
                continue
            if not update_cache and is_url_in_cache(url, pkgname=pkgname):
                continue
            f_temp_name = os.path.join(d, str(i))
            with z.open(zf) as f_zip, open(f_temp_name, "wb") as f_temp:
                hasher = hashlib.md5()
                block = f_zip.read(conf.download_block_size)
                while block:
                    f_temp.write(block)
                    hasher.update(block)
                    block = f_zip.read(conf.download_block_size)
                hexdigest = hasher.hexdigest()
            import_file_to_cache(url, f_temp_name,
                                 hexdigest=hexdigest,
                                 remove_original=True,
                                 pkgname=pkgname)
