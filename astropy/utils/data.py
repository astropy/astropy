# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Functions for accessing, downloading, and caching data files."""

import atexit
import contextlib
import errno
import fnmatch
import ftplib
import functools
import hashlib
import io
import os
import re
import shutil
import sys
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from importlib import import_module
from tempfile import NamedTemporaryFile, TemporaryDirectory, gettempdir
from types import MappingProxyType
from warnings import warn

import astropy_iers_data

import astropy.config.paths
from astropy import config as _config
from astropy.utils.compat.optional_deps import (
    HAS_BZ2,
    HAS_CERTIFI,
    HAS_FSSPEC,
    HAS_LZMA,
    HAS_UNCOMPRESSPY,
)
from astropy.utils.exceptions import AstropyDeprecationWarning, AstropyWarning
from astropy.utils.introspection import find_current_module

# Order here determines order in the autosummary
__all__ = [
    "CacheDamaged",
    "CacheMissingWarning",
    "Conf",
    "cache_contents",
    "cache_total_size",
    "check_download_cache",
    "check_free_space_in_dir",
    "clear_download_cache",
    "compute_hash",
    "conf",
    "download_file",
    "download_files_in_parallel",
    "export_download_cache",
    "get_cached_urls",
    "get_file_contents",
    "get_free_space_in_dir",
    "get_pkg_data_contents",
    "get_pkg_data_filename",
    "get_pkg_data_filenames",
    "get_pkg_data_fileobj",
    "get_pkg_data_fileobjs",
    "get_pkg_data_path",
    "get_readable_fileobj",
    "import_download_cache",
    "import_file_to_cache",
    "is_url",
    "is_url_in_cache",
]

_dataurls_to_alias = {}


_IERS_DATA_REDIRECTS = {
    "Leap_Second.dat": (
        "IERS_LEAP_SECOND_FILE",
        astropy_iers_data.IERS_LEAP_SECOND_FILE,
    ),
    "ReadMe.finals2000A": ("IERS_A_README", astropy_iers_data.IERS_A_README),
    "ReadMe.eopc04": ("IERS_B_README", astropy_iers_data.IERS_B_README),
    "eopc04.1962-now": ("IERS_B_FILE", astropy_iers_data.IERS_B_FILE),
}


class _NonClosingBufferedReader(io.BufferedReader):
    def __del__(self):
        try:
            # NOTE: self.raw will not be closed, but left in the state
            # it was in at detactment
            self.detach()
        except Exception:
            pass


class _NonClosingTextIOWrapper(io.TextIOWrapper):
    def __del__(self):
        try:
            # NOTE: self.stream will not be closed, but left in the state
            # it was in at detactment
            self.detach()
        except Exception:
            pass


class Conf(_config.ConfigNamespace):
    """
    Configuration parameters for `astropy.utils.data`.
    """

    dataurl = _config.ConfigItem(
        "http://data.astropy.org/", "Primary URL for astropy remote data site."
    )
    dataurl_mirror = _config.ConfigItem(
        "http://www.astropy.org/astropy-data/",
        "Mirror URL for astropy remote data site.",
    )
    default_http_user_agent = _config.ConfigItem(
        "astropy",
        "Default User-Agent for HTTP request headers. This can be overwritten "
        "for a particular call via http_headers option, where available. "
        "This only provides the default value when not set by https_headers.",
    )
    remote_timeout = _config.ConfigItem(
        10.0,
        "Time to wait for remote data queries (in seconds).",
        aliases=["astropy.coordinates.name_resolve.name_resolve_timeout"],
    )
    allow_internet = _config.ConfigItem(
        True, "If False, prevents any attempt to download from Internet."
    )
    compute_hash_block_size = _config.ConfigItem(
        2**16,  # 64K
        "Block size for computing file hashes.",
    )
    download_block_size = _config.ConfigItem(
        2**16,  # 64K
        "Number of bytes of remote data to download per step.",
    )
    delete_temporary_downloads_at_exit = _config.ConfigItem(
        True,
        "If True, temporary download files created when the cache is "
        "inaccessible will be deleted at the end of the python session.",
    )


conf = Conf()


class CacheMissingWarning(AstropyWarning):
    """
    This warning indicates the standard cache directory is not accessible, with
    the first argument providing the warning message. If args[1] is present, it
    is a filename indicating the path to a temporary file that was created to
    store a remote data download in the absence of the cache.
    """


def is_url(string):
    """
    Test whether a string is a valid URL for :func:`download_file`.

    Parameters
    ----------
    string : str
        The string to test.

    Returns
    -------
    status : bool
        String is URL or not.

    """
    url = urllib.parse.urlparse(string)
    # we can't just check that url.scheme is not an empty string, because
    # file paths in windows would return a non-empty scheme (e.g. e:\\
    # returns 'e').
    return url.scheme.lower() in ["http", "https", "ftp", "sftp", "ssh", "file"]


# Backward compatibility because some downstream packages allegedly uses it.
_is_url = is_url


def _requires_fsspec(url):
    """Does the `url` require the optional ``fsspec`` dependency to open?"""
    return isinstance(url, str) and url.startswith(("s3://", "gs://"))


def _is_inside(path, parent_path):
    # We have to try realpath too to avoid issues with symlinks, but we leave
    # abspath because some systems like debian have the absolute path (with no
    # symlinks followed) match, but the real directories in different
    # locations, so need to try both cases.
    return os.path.abspath(path).startswith(
        os.path.abspath(parent_path)
    ) or os.path.realpath(path).startswith(os.path.realpath(parent_path))


@contextlib.contextmanager
def get_readable_fileobj(
    name_or_obj,
    encoding=None,
    cache=False,
    show_progress=True,
    remote_timeout=None,
    sources=None,
    http_headers=None,
    *,
    use_fsspec=None,
    fsspec_kwargs=None,
    close_files=True,
):
    """Yield a readable, seekable file-like object from a file or URL.

    This supports passing filenames, URLs, and readable file-like objects,
    any of which can be compressed in gzip, bzip2, lzma (xz) or lzw (Z) if the
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
    name_or_obj : str or file-like
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
        `astropy.utils.data.Conf.remote_timeout`).

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

    use_fsspec : bool, optional
        Use `fsspec.open` to open the file? Defaults to `False` unless
        ``name_or_obj`` starts with the Amazon S3 storage prefix ``s3://``
        or the Google Cloud Storage prefix ``gs://``.  Can also be used for paths
        with other prefixes (e.g. ``http://``) but in this case you must
        explicitly pass ``use_fsspec=True``.
        Use of this feature requires the optional ``fsspec`` package.
        A ``ModuleNotFoundError`` will be raised if the dependency is missing.

        .. versionadded:: 5.2

    fsspec_kwargs : dict, optional
        Keyword arguments passed on to `fsspec.open`. This can be used to
        configure cloud storage credentials and caching behavior.
        For example, pass ``fsspec_kwargs={"anon": True}`` to enable
        anonymous access to Amazon S3 open data buckets.
        See ``fsspec``'s documentation for available parameters.

        .. versionadded:: 5.2

    close_files : bool, optional
        Close the file object when exiting the context manager.
        Default is `True`.

        .. versionadded:: 5.2

    Returns
    -------
    file : :term:`file-like (readable)`
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

    # Have `use_fsspec` default to ``True`` if the user passed an Amazon S3
    # or Google Cloud Storage URI.
    if use_fsspec is None and _requires_fsspec(name_or_obj):
        use_fsspec = True

    if use_fsspec:
        if not isinstance(name_or_obj, str):
            raise TypeError("`name_or_obj` must be a string when `use_fsspec=True`")
        if fsspec_kwargs is None:
            fsspec_kwargs = {}

    # name_or_obj could be an os.PathLike object
    if isinstance(name_or_obj, os.PathLike):
        name_or_obj = os.fspath(name_or_obj)

    # Get a file object to the content
    if isinstance(name_or_obj, str):
        # Use fsspec to open certain cloud-hosted files (e.g., AWS S3, Google Cloud Storage)
        if use_fsspec:
            if not HAS_FSSPEC:
                raise ModuleNotFoundError("please install `fsspec` to open this file")
            import fsspec  # local import because it is a niche dependency

            openfileobj = fsspec.open(name_or_obj, **fsspec_kwargs)
            close_fds.append(openfileobj)
            fileobj = openfileobj.open()
            close_fds.append(fileobj)
        else:
            is_url = _is_url(name_or_obj)
            if is_url:
                name_or_obj = download_file(
                    name_or_obj,
                    cache=cache,
                    show_progress=show_progress,
                    timeout=remote_timeout,
                    sources=sources,
                    http_headers=http_headers,
                )
            fileobj = io.FileIO(name_or_obj, "r")
            if is_url and not cache:
                delete_fds.append(fileobj)
            close_fds.append(fileobj)
    else:
        fileobj = name_or_obj

    # Check if the file object supports random access, and if not,
    # then wrap it in a BytesIO buffer.  It would be nicer to use a
    # BufferedReader to avoid reading loading the whole file first,
    # but that might not be compatible with all possible I/O classes.
    if not hasattr(fileobj, "seek"):
        try:
            # py.path.LocalPath objects have .read() method but it uses
            # text mode, which won't work. .read_binary() does, and
            # surely other ducks would return binary contents when
            # called like this.
            # py.path.LocalPath is what comes from the legacy tmpdir fixture
            # in pytest.
            fileobj = io.BytesIO(fileobj.read_binary())
        except AttributeError:
            fileobj = io.BytesIO(fileobj.read())

    # Now read enough bytes to look at signature
    signature = fileobj.read(6)
    fileobj.seek(0)

    if signature[:3] == b"\x1f\x8b\x08":  # gzip
        import struct

        try:
            import gzip

            fileobj_new = gzip.GzipFile(fileobj=fileobj, mode="rb")
            fileobj_new.read(1)  # need to check that the file is really gzip
        except (OSError, EOFError, struct.error):  # invalid gzip file
            fileobj.seek(0)
            fileobj_new.close()
        else:
            fileobj_new.seek(0)
            fileobj = fileobj_new
    elif signature[:3] == b"BZh":  # bzip2
        if not HAS_BZ2:
            for fd in close_fds:
                fd.close()
            raise ModuleNotFoundError(
                "This Python installation does not provide the bz2 module."
            )
        import bz2

        try:
            # bz2.BZ2File does not support file objects, only filenames, so we
            # need to write the data to a temporary file
            with NamedTemporaryFile("wb", delete=False) as tmp:
                tmp.write(fileobj.read())
                tmp.close()
                fileobj_new = bz2.BZ2File(tmp.name, mode="rb")
            fileobj_new.read(1)  # need to check that the file is really bzip2
        except OSError:  # invalid bzip2 file
            fileobj.seek(0)
            fileobj_new.close()
            # raise
        else:
            fileobj_new.seek(0)
            close_fds.append(fileobj_new)
            fileobj = fileobj_new
    elif signature[:6] == b"\xfd7zXZ\x00":  # xz
        if not HAS_LZMA:
            for fd in close_fds:
                fd.close()
            raise ModuleNotFoundError(
                "This Python installation does not provide the lzma module."
            )
        import lzma

        try:
            fileobj_new = lzma.LZMAFile(fileobj, mode="rb")
            fileobj_new.read(1)  # need to check that the file is really xz
        except lzma.LZMAError:  # invalid xz file
            fileobj.seek(0)
            fileobj_new.close()
            # should we propagate this to the caller to signal bad content?
            # raise ValueError(e)
        else:
            fileobj_new.seek(0)
            fileobj = fileobj_new
    elif signature[:2] == b"\x1f\x9d":  # LZW
        if not HAS_UNCOMPRESSPY:
            for fd in close_fds:
                fd.close()
            raise ModuleNotFoundError(
                "The optional package uncompresspy is necessary for reading LZW"
                " compressed files (.Z extension)."
            )
        import uncompresspy

        try:
            fileobj_new = uncompresspy.LZWFile(fileobj)
            fileobj_new.read(1)
        except ValueError:
            fileobj.seek(0)
            fileobj_new.close()
        else:
            fileobj_new.seek(0)
            close_fds.append(fileobj)
            fileobj = fileobj_new

    # By this point, we have a file, io.FileIO, gzip.GzipFile, bz2.BZ2File,
    # lzma.LZMAFile or uncompresspy.LZWFile instance opened in binary mode (that
    # is, read returns bytes). Now we need to, if requested, wrap it in a
    # io.TextIOWrapper so read will return unicode based on the
    # encoding parameter.

    needs_textio_wrapper = encoding != "binary"

    if needs_textio_wrapper:
        # A bz2.BZ2File can not be wrapped by a TextIOWrapper,
        # so we decompress it to a temporary file and then
        # return a handle to that.
        if HAS_BZ2:
            import bz2

            if isinstance(fileobj, bz2.BZ2File):
                tmp = NamedTemporaryFile("wb", delete=False)
                data = fileobj.read()
                tmp.write(data)
                tmp.close()
                delete_fds.append(tmp)

                fileobj = io.FileIO(tmp.name, "r")
                close_fds.append(fileobj)

        fileobj = _NonClosingBufferedReader(fileobj)
        fileobj = _NonClosingTextIOWrapper(fileobj, encoding=encoding)

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
        if close_files:
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
    object
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
    datafn = get_pkg_data_path(data_name, package=package)
    if os.path.isdir(datafn):
        raise OSError(
            "Tried to access a data file that's actually a package data directory"
        )
    elif os.path.isfile(datafn):  # local file
        with get_readable_fileobj(datafn, encoding=encoding) as fileobj:
            yield fileobj
    else:  # remote file
        with get_readable_fileobj(
            conf.dataurl + data_name,
            encoding=encoding,
            cache=cache,
            sources=[conf.dataurl + data_name, conf.dataurl_mirror + data_name],
        ) as fileobj:
            # We read a byte to trigger any URLErrors
            fileobj.read(1)
            fileobj.seek(0)
            yield fileobj


def get_pkg_data_filename(
    data_name, package=None, show_progress=True, remote_timeout=None
):
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
        configurable `astropy.utils.data.Conf.remote_timeout`).

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

    if data_name.startswith("hash/"):
        # first try looking for a local version if a hash is specified
        hashfn = _find_hash_fn(data_name[5:])

        if hashfn is None:
            return download_file(
                conf.dataurl + data_name,
                cache=True,
                show_progress=show_progress,
                timeout=remote_timeout,
                sources=[conf.dataurl + data_name, conf.dataurl_mirror + data_name],
            )
        else:
            return hashfn
    else:
        fs_path = os.path.normpath(data_name)
        datafn = get_pkg_data_path(fs_path, package=package)
        if os.path.isdir(datafn):
            raise OSError(
                "Tried to access a data file that's actually a package data directory"
            )
        elif os.path.isfile(datafn):  # local file
            return datafn
        else:  # remote file
            return download_file(
                conf.dataurl + data_name,
                cache=True,
                show_progress=show_progress,
                timeout=remote_timeout,
                sources=[conf.dataurl + data_name, conf.dataurl_mirror + data_name],
            )


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
    with get_pkg_data_fileobj(
        data_name, package=package, encoding=encoding, cache=cache
    ) as fd:
        contents = fd.read()
    return contents


def get_pkg_data_filenames(datadir, package=None, pattern="*"):
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
    path = get_pkg_data_path(datadir, package=package)
    if os.path.isfile(path):
        raise OSError(
            "Tried to access a data directory that's actually a package data file"
        )
    elif os.path.isdir(path):
        for filename in os.listdir(path):
            if fnmatch.fnmatch(filename, pattern):
                yield os.path.join(path, filename)
    else:
        raise OSError("Path not found")


def get_pkg_data_fileobjs(datadir, package=None, pattern="*", encoding=None):
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
    fileobjs : iterator of file object
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
    for fn in get_pkg_data_filenames(datadir, package=package, pattern=pattern):
        with get_readable_fileobj(fn, encoding=encoding) as fd:
            yield fd


def compute_hash(localfn):
    """Computes the MD5 hash for a file.

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
    with open(localfn, "rb") as f:
        h = hashlib.md5(usedforsecurity=False)
        block = f.read(conf.compute_hash_block_size)
        while block:
            h.update(block)
            block = f.read(conf.compute_hash_block_size)

    return h.hexdigest()


def get_pkg_data_path(*path, package=None):
    """Get path from source-included data directories.

    Parameters
    ----------
    *path : str
        Name/location of the desired data file/directory.
        May be a tuple of strings for ``os.path`` joining.

    package : str or None, optional, keyword-only
        If specified, look for a file relative to the given package, rather
        than the calling module's package.

    Returns
    -------
    path : str
        Name/location of the desired data file/directory.

    Raises
    ------
    ImportError
        Given package or module is not importable.
    RuntimeError
        If the local data file is outside of the package's tree.

    """
    if package is None:
        module = find_current_module(1, finddiff=["astropy.utils.data", "contextlib"])
        if module is None:
            # not called from inside an astropy package.  So just pass name
            # through
            return os.path.join(*path)

        if not hasattr(module, "__package__") or not module.__package__:
            # The __package__ attribute may be missing or set to None; see
            # PEP-366, also astropy issue #1256
            if "." in module.__name__:
                package = module.__name__.rpartition(".")[0]
            else:
                package = module.__name__
        else:
            package = module.__package__
    else:
        # Backward-compatibility for files that used to exist in astropy.utils.iers
        if package == "astropy.utils.iers":
            filename = os.path.basename(path[-1])
            if filename in _IERS_DATA_REDIRECTS:
                warn(
                    f"Accessing {filename} in this way is deprecated in v6.0, "
                    f"use astropy.utils.iers.{_IERS_DATA_REDIRECTS[filename][0]} "
                    "instead.",
                    AstropyDeprecationWarning,
                )
                return _IERS_DATA_REDIRECTS[filename][1]

        # package errors if it isn't a str
        # so there is no need for checks in the containing if/else
        module = import_module(package)

    # module path within package
    module_path = os.path.dirname(module.__file__)
    full_path = os.path.join(module_path, *path)

    # Check that file is inside tree.
    rootpkgname = package.partition(".")[0]
    root_dir = os.path.dirname(import_module(rootpkgname).__file__)
    if not _is_inside(full_path, root_dir):
        raise RuntimeError(
            f"attempted to get a local data file outside of the {rootpkgname} tree."
        )

    return full_path


def _find_hash_fn(hexdigest, pkgname="astropy"):
    """
    Looks for a local file by hash - returns file name if found and a valid
    file, otherwise returns None.
    """
    for v in cache_contents(pkgname=pkgname).values():
        if compute_hash(v) == hexdigest:
            return v
    return None


def get_free_space_in_dir(path, unit=False):
    """
    Given a path to a directory, returns the amount of free space
    on that filesystem.

    Parameters
    ----------
    path : str
        The path to a directory.

    unit : bool or `~astropy.units.Unit`
        Return the amount of free space as Quantity in the given unit,
        if provided. Default is `False` for backward-compatibility.

    Returns
    -------
    free_space : int or `~astropy.units.Quantity`
        The amount of free space on the partition that the directory is on.
        If ``unit=False``, it is returned as plain integer (in bytes).

    """
    if not os.path.isdir(path):
        raise OSError(
            "Can only determine free space associated with directories, not files."
        )
        # Actually you can on Linux but I want to avoid code that fails
        # on Windows only.
    free_space = shutil.disk_usage(path).free
    if unit:
        from astropy import units as u

        # TODO: Automatically determine best prefix to use.
        if unit is True:
            unit = u.byte
        free_space = u.Quantity(free_space, u.byte).to(unit)
    return free_space


def check_free_space_in_dir(path, size):
    """
    Determines if a given directory has enough space to hold a file of
    a given size.

    Parameters
    ----------
    path : str
        The path to a directory.

    size : int or `~astropy.units.Quantity`
        A proposed filesize. If not a Quantity, assume it is in bytes.

    Raises
    ------
    OSError
        There is not enough room on the filesystem.
    """
    space = get_free_space_in_dir(path, unit=getattr(size, "unit", False))
    if space < size:
        from astropy.utils.console import human_file_size

        raise OSError(
            f"Not enough free space in {path} "
            f"to download a {human_file_size(size)} file, "
            f"only {human_file_size(space)} left"
        )


class _ftptlswrapper(urllib.request.ftpwrapper):
    def init(self):
        self.busy = 0
        self.ftp = ftplib.FTP_TLS()
        self.ftp.connect(self.host, self.port, self.timeout)
        self.ftp.login(self.user, self.passwd)
        self.ftp.prot_p()
        _target = "/".join(self.dirs)
        self.ftp.cwd(_target)


class _FTPTLSHandler(urllib.request.FTPHandler):
    def connect_ftp(self, user, passwd, host, port, dirs, timeout):
        return _ftptlswrapper(user, passwd, host, port, dirs, timeout, persistent=False)


@functools.lru_cache
def _build_urlopener(ftp_tls=False, ssl_context=None, allow_insecure=False):
    """
    Helper for building a `urllib.request.build_opener` which handles TLS/SSL.
    """
    # Import ssl here to avoid import failure when running in pyodide/Emscripten
    import ssl

    ssl_context = dict(it for it in ssl_context) if ssl_context else {}
    cert_chain = {}
    if "certfile" in ssl_context:
        cert_chain.update(
            {
                "certfile": ssl_context.pop("certfile"),
                "keyfile": ssl_context.pop("keyfile", None),
                "password": ssl_context.pop("password", None),
            }
        )
    elif "password" in ssl_context or "keyfile" in ssl_context:
        raise ValueError(
            "passing 'keyfile' or 'password' in the ssl_context argument "
            "requires passing 'certfile' as well"
        )

    if "cafile" not in ssl_context and HAS_CERTIFI:
        import certifi

        ssl_context["cafile"] = certifi.where()

    ssl_context = ssl.create_default_context(**ssl_context)

    if allow_insecure:
        ssl_context.check_hostname = False
        ssl_context.verify_mode = ssl.CERT_NONE

    if cert_chain:
        ssl_context.load_cert_chain(**cert_chain)

    https_handler = urllib.request.HTTPSHandler(context=ssl_context)

    if ftp_tls:
        urlopener = urllib.request.build_opener(_FTPTLSHandler(), https_handler)
    else:
        urlopener = urllib.request.build_opener(https_handler)

    return urlopener


def _try_url_open(
    source_url,
    timeout=None,
    http_headers=None,
    ftp_tls=False,
    ssl_context=None,
    allow_insecure=False,
):
    """Helper for opening a URL while handling TLS/SSL verification issues."""
    # Import ssl here to avoid import failure when running in pyodide/Emscripten
    import ssl

    # Always try first with a secure connection
    # _build_urlopener uses lru_cache, so the ssl_context argument must be
    # converted to a hashshable type (a set of 2-tuples)
    ssl_context = frozenset(ssl_context.items() if ssl_context else [])
    urlopener = _build_urlopener(
        ftp_tls=ftp_tls, ssl_context=ssl_context, allow_insecure=False
    )
    req = urllib.request.Request(source_url, headers=http_headers)

    try:
        return urlopener.open(req, timeout=timeout)
    except urllib.error.URLError as exc:
        reason = exc.reason
        if (
            isinstance(reason, ssl.SSLError)
            and reason.reason == "CERTIFICATE_VERIFY_FAILED"
        ):
            msg = (
                f"Verification of TLS/SSL certificate at {source_url} "
                "failed: this can mean either the server is "
                "misconfigured or your local root CA certificates are "
                "out-of-date; in the latter case this can usually be "
                'addressed by installing the Python package "certifi" '
                "(see the documentation for astropy.utils.data.download_file)"
            )
            if not allow_insecure:
                msg += (
                    " or in both cases you can work around this by "
                    "passing allow_insecure=True, but only if you "
                    "understand the implications; the original error "
                    f"was: {reason}"
                )
                raise urllib.error.URLError(msg)
            else:
                msg += ". Re-trying with allow_insecure=True."
                warn(msg, AstropyWarning)
                # Try again with a new urlopener allowing insecure connections
                urlopener = _build_urlopener(
                    ftp_tls=ftp_tls, ssl_context=ssl_context, allow_insecure=True
                )
                return urlopener.open(req, timeout=timeout)

        raise


def _download_file_from_source(
    source_url,
    show_progress=True,
    timeout=None,
    remote_url=None,
    cache=False,
    pkgname="astropy",
    http_headers=None,
    ftp_tls=None,
    ssl_context=None,
    allow_insecure=False,
):
    from astropy.utils.console import ProgressBarOrSpinner

    if not conf.allow_internet:
        raise urllib.error.URLError(
            f"URL {remote_url} was supposed to be downloaded but "
            f"allow_internet is {conf.allow_internet}; "
            "if this is unexpected check the astropy.cfg file for the option "
            "allow_internet"
        )

    if remote_url is None:
        remote_url = source_url
    if http_headers is None:
        http_headers = {}

    if ftp_tls is None and urllib.parse.urlparse(remote_url).scheme == "ftp":
        try:
            return _download_file_from_source(
                source_url,
                show_progress=show_progress,
                timeout=timeout,
                remote_url=remote_url,
                cache=cache,
                pkgname=pkgname,
                http_headers=http_headers,
                ftp_tls=False,
            )
        except urllib.error.URLError as e:
            # e.reason might not be a string, e.g. socket.gaierror
            # URLError changed to report original exception in Python 3.11 (bpo-43564)
            if (
                str(e.reason)
                .removeprefix("ftp error: ")
                .startswith(("error_perm", "5"))
            ):
                ftp_tls = True
            else:
                raise

    with _try_url_open(
        source_url,
        timeout=timeout,
        http_headers=http_headers,
        ftp_tls=ftp_tls,
        ssl_context=ssl_context,
        allow_insecure=allow_insecure,
    ) as remote:
        info = remote.info()
        try:
            size = int(info["Content-Length"])
        except (KeyError, ValueError, TypeError):
            size = None

        if size is not None:
            check_free_space_in_dir(gettempdir(), size)
            if cache:
                dldir = _get_download_cache_loc(pkgname)
                check_free_space_in_dir(dldir, size)

        # If a user has overridden sys.stdout it might not have the
        # isatty method, in that case assume it's not a tty
        is_tty = hasattr(sys.stdout, "isatty") and sys.stdout.isatty()
        if show_progress and is_tty:
            progress_stream = sys.stdout
        else:
            progress_stream = io.StringIO()

        if source_url == remote_url:
            dlmsg = f"Downloading {remote_url}"
        else:
            dlmsg = f"Downloading {remote_url} from {source_url}"
        with ProgressBarOrSpinner(size, dlmsg, file=progress_stream) as p:
            with NamedTemporaryFile(
                prefix=f"astropy-download-{os.getpid()}-", delete=False
            ) as f:
                try:
                    bytes_read = 0
                    block = remote.read(conf.download_block_size)
                    while block:
                        f.write(block)
                        bytes_read += len(block)
                        p.update(bytes_read)
                        block = remote.read(conf.download_block_size)
                        if size is not None and bytes_read > size:
                            raise urllib.error.URLError(
                                f"File was supposed to be {size} bytes but "
                                f"server provides more, at least {bytes_read} "
                                "bytes. Download failed."
                            )
                    if size is not None and bytes_read < size:
                        raise urllib.error.ContentTooShortError(
                            f"File was supposed to be {size} bytes but we "
                            f"only got {bytes_read} bytes. Download failed.",
                            content=None,
                        )
                except BaseException:
                    if os.path.exists(f.name):
                        try:
                            os.remove(f.name)
                        except OSError:
                            pass
                    raise
    return f.name


def download_file(
    remote_url,
    cache=False,
    show_progress=True,
    timeout=None,
    sources=None,
    pkgname="astropy",
    http_headers=None,
    ssl_context=None,
    allow_insecure=False,
):
    """Downloads a URL and optionally caches the result.

    It returns the filename of a file containing the URL's contents.
    If ``cache=True`` and the file is present in the cache, just
    returns the filename; if the file had to be downloaded, add it
    to the cache. If ``cache="update"`` always download and add it
    to the cache.

    The cache is effectively a dictionary mapping URLs to files; by default the
    file contains the contents of the URL that is its key, but in practice
    these can be obtained from a mirror (using ``sources``) or imported from
    the local filesystem (using `~import_file_to_cache` or
    `~import_download_cache`).  Regardless, each file is regarded as
    representing the contents of a particular URL, and this URL should be used
    to look them up or otherwise manipulate them.

    The files in the cache directory are named according to a cryptographic
    hash of their URLs (currently MD5, so hackers can cause collisions).
    The modification times on these files normally indicate when they were
    last downloaded from the Internet.

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
        Timeout for remote requests in seconds (default is the configurable
        `astropy.utils.data.Conf.remote_timeout`).

    sources : list of str, optional
        If provided, a list of URLs to try to obtain the file from. The
        result will be stored under the original URL. The original URL
        will *not* be tried unless it is in this list; this is to prevent
        long waits for a primary server that is known to be inaccessible
        at the moment. If an empty list is passed, then ``download_file``
        will not attempt to connect to the Internet, that is, if the file
        is not in the cache a KeyError will be raised.

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

    ssl_context : dict, optional
        Keyword arguments to pass to `ssl.create_default_context` when
        downloading from HTTPS or TLS+FTP sources.  This can be used provide
        alternative paths to root CA certificates.  Additionally, if the key
        ``'certfile'`` and optionally ``'keyfile'`` and ``'password'`` are
        included, they are passed to `ssl.SSLContext.load_cert_chain`.  This
        can be used for performing SSL/TLS client certificate authentication
        for servers that require it.

    allow_insecure : bool, optional
        Allow downloading files over a TLS/SSL connection even when the server
        certificate verification failed.  When set to `True` the potentially
        insecure download is allowed to proceed, but an
        `~astropy.utils.exceptions.AstropyWarning` is issued.  If you are
        frequently getting certificate verification warnings, consider
        installing or upgrading `certifi`_ package, which provides frequently
        updated certificates for common root CAs (i.e., a set similar to those
        used by web browsers).  If installed, Astropy will use it
        automatically.

        .. _certifi: https://pypi.org/project/certifi/

    Returns
    -------
    local_path : str
        Returns the local path that the file was download to.

    Raises
    ------
    urllib.error.URLError
        Whenever there's a problem getting the remote file.
    KeyError
        When a file was requested from the cache but is missing and no
        sources were provided to obtain it from the Internet.

    Notes
    -----
    Because this function returns a filename, another process could run
    `clear_download_cache` before you actually open the file, leaving
    you with a filename that no longer points to a usable file.
    """
    if timeout is None:
        timeout = conf.remote_timeout
    if sources is None:
        sources = [remote_url]
    if http_headers is None:
        http_headers = {"User-Agent": conf.default_http_user_agent, "Accept": "*/*"}

    missing_cache = ""

    url_key = remote_url

    if cache:
        try:
            dldir = _get_download_cache_loc(pkgname)
        except OSError as e:
            cache = False
            missing_cache = (
                f"Cache directory cannot be read or created ({e}), "
                "providing data in temporary file instead."
            )
        else:
            if cache == "update":
                pass
            elif isinstance(cache, str):
                raise ValueError(
                    f"Cache value '{cache}' was requested but "
                    "'update' is the only recognized string; "
                    "otherwise use a boolean"
                )
            else:
                filename = os.path.join(dldir, _url_to_dirname(url_key), "contents")
                if os.path.exists(filename):
                    return os.path.abspath(filename)

    errors = {}
    for source_url in sources:
        try:
            f_name = _download_file_from_source(
                source_url,
                timeout=timeout,
                show_progress=show_progress,
                cache=cache,
                remote_url=remote_url,
                pkgname=pkgname,
                http_headers=http_headers,
                ssl_context=ssl_context,
                allow_insecure=allow_insecure,
            )
            # Success!
            break

        except urllib.error.URLError as e:
            # errno 8 is from SSL "EOF occurred in violation of protocol"
            if (
                hasattr(e, "reason")
                and hasattr(e.reason, "errno")
                and e.reason.errno == 8
            ):
                e.reason.strerror = f"{e.reason.strerror}. requested URL: {remote_url}"
                e.reason.args = (e.reason.errno, e.reason.strerror)
            errors[source_url] = e

        except TimeoutError as e:
            errors[source_url] = e

    else:  # No success
        if not sources:
            raise KeyError(
                f"No sources listed and file {remote_url} not in cache! "
                "Please include primary URL in sources if you want it to be "
                "included as a valid source."
            )
        elif len(sources) == 1:
            raise errors[sources[0]]
        else:
            raise urllib.error.URLError(
                f"Unable to open any source! Exceptions were {errors}"
            ) from errors[sources[0]]

    if cache:
        try:
            return import_file_to_cache(
                url_key,
                f_name,
                remove_original=True,
                replace=(cache == "update"),
                pkgname=pkgname,
            )
        except PermissionError as e:
            # Cache is readonly, we can't update it
            missing_cache = (
                f"Cache directory appears to be read-only ({e}), unable to import "
                f"downloaded file, providing data in temporary file {f_name} "
                "instead."
            )
        # FIXME: other kinds of cache problem can occur?

    if missing_cache:
        warn(CacheMissingWarning(missing_cache, f_name))
    if conf.delete_temporary_downloads_at_exit:
        _tempfilestodel.append(f_name)
    return os.path.abspath(f_name)


def is_url_in_cache(url_key, pkgname="astropy"):
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
    try:
        dldir = _get_download_cache_loc(pkgname)
    except OSError:
        return False
    filename = os.path.join(dldir, _url_to_dirname(url_key), "contents")
    return os.path.exists(filename)


def cache_total_size(pkgname="astropy"):
    """Return the total size in bytes of all files in the cache."""
    size = 0
    dldir = _get_download_cache_loc(pkgname=pkgname)
    for root, _, files in os.walk(dldir):
        size += sum(os.path.getsize(os.path.join(root, name)) for name in files)
    return size


def _do_download_files_in_parallel(kwargs):
    with astropy.config.paths.set_temp_config(kwargs.pop("temp_config")):
        with astropy.config.paths.set_temp_cache(kwargs.pop("temp_cache")):
            return download_file(**kwargs)


def download_files_in_parallel(
    urls,
    cache="update",
    show_progress=True,
    timeout=None,
    sources=None,
    multiprocessing_start_method=None,
    pkgname="astropy",
):
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
        warn(
            "Disabling the cache does not work because of multiprocessing, "
            'it will be set to ``"update"``. You may need to manually remove '
            "the cached files with clear_download_cache() afterwards.",
            AstropyWarning,
        )
        cache = "update"

    if show_progress:
        progress = sys.stdout
    else:
        progress = io.BytesIO()

    # Combine duplicate URLs
    combined_urls = list(set(urls))
    combined_paths = ProgressBar.map(
        _do_download_files_in_parallel,
        [
            dict(
                remote_url=u,
                cache=cache,
                show_progress=False,
                timeout=timeout,
                sources=sources.get(u, None),
                pkgname=pkgname,
                temp_cache=astropy.config.paths.set_temp_cache._temp_path,
                temp_config=astropy.config.paths.set_temp_config._temp_path,
            )
            for u in combined_urls
        ],
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
            elif os.path.isdir(fn):
                try:
                    shutil.rmtree(fn)
                except OSError:
                    # couldn't get rid of it, sorry
                    # could be held open by some process, on Windows
                    pass


def clear_download_cache(hashorurl=None, pkgname="astropy"):
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
    try:
        dldir = _get_download_cache_loc(pkgname)
    except OSError as e:
        # Problem arose when trying to open the cache
        # Just a warning, though
        msg = "Not clearing data cache - cache inaccessible due to "
        estr = "" if len(e.args) < 1 else (": " + str(e))
        warn(CacheMissingWarning(msg + e.__class__.__name__ + estr))
        return
    try:
        if hashorurl is None:
            # Optional: delete old incompatible caches too
            _rmtree(dldir)
        elif _is_url(hashorurl):
            filepath = os.path.join(dldir, _url_to_dirname(hashorurl))
            _rmtree(filepath)
        else:
            # Not a URL, it should be either a filename or a hash
            filepath = os.path.join(dldir, hashorurl)
            rp = os.path.relpath(filepath, dldir)
            if rp.startswith(".."):
                raise RuntimeError(
                    "attempted to use clear_download_cache on the path "
                    f"{filepath} outside the data cache directory {dldir}"
                )
            d, f = os.path.split(rp)
            if d and f in ["contents", "url"]:
                # It's a filename not the hash of a URL
                # so we want to zap the directory containing the
                # files "url" and "contents"
                filepath = os.path.join(dldir, d)
            if os.path.exists(filepath):
                _rmtree(filepath)
            elif len(hashorurl) == 2 * hashlib.md5(
                usedforsecurity=False
            ).digest_size and re.match(r"[0-9a-f]+", hashorurl):
                # It's the hash of some file contents, we have to find the right file
                filename = _find_hash_fn(hashorurl)
                if filename is not None:
                    clear_download_cache(filename)
    except OSError as e:
        msg = "Not clearing data from cache - problem arose "
        estr = "" if len(e.args) < 1 else (": " + str(e))
        warn(CacheMissingWarning(msg + e.__class__.__name__ + estr))


def _get_download_cache_loc(pkgname="astropy"):
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
    """
    try:
        datadir = astropy.config.paths.get_cache_dir_path(pkgname) / "download" / "url"

        if not datadir.exists():
            try:
                datadir.mkdir(parents=True)
            except OSError:
                if not datadir.exists():
                    raise
        elif not datadir.is_dir():
            raise OSError(f"Data cache directory {datadir} is not a directory")

        return datadir
    except OSError as e:
        msg = "Remote data cache could not be accessed due to "
        estr = "" if len(e.args) < 1 else (": " + str(e))
        warn(CacheMissingWarning(msg + e.__class__.__name__ + estr))
        raise


def _url_to_dirname(url):
    if not _is_url(url):
        raise ValueError(f"Malformed URL: '{url}'")
    # Make domain names case-insensitive
    # Also makes the http:// case-insensitive
    urlobj = list(urllib.parse.urlsplit(url))
    urlobj[1] = urlobj[1].lower()
    if urlobj[0].lower() in ["http", "https"] and urlobj[1] and urlobj[2] == "":
        urlobj[2] = "/"
    url_c = urllib.parse.urlunsplit(urlobj)
    return hashlib.md5(url_c.encode("utf-8"), usedforsecurity=False).hexdigest()


_NOTHING = MappingProxyType({})


class CacheDamaged(ValueError):
    """Record the URL or file that was a problem.
    Using clear_download_cache on the .bad_file or .bad_url attribute,
    whichever is not None, should resolve this particular problem.
    """

    def __init__(self, *args, bad_urls=None, bad_files=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.bad_urls = bad_urls if bad_urls is not None else []
        self.bad_files = bad_files if bad_files is not None else []


def check_download_cache(pkgname="astropy"):
    """Do a consistency check on the cache.

    .. note::

        Since v5.0, this function no longer returns anything.

    Because the cache is shared by all versions of ``astropy`` in all virtualenvs
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
    :func:`clear_download_cache`, either passing the filename returned here, or
    with no arguments to empty the entire cache and return it to a
    reasonable, if empty, state.

    Parameters
    ----------
    pkgname : str, optional
        The package name to use to locate the download cache, i.e., for
        ``pkgname='astropy'`` the default cache location is
        ``~/.astropy/cache``.

    Raises
    ------
    `~astropy.utils.data.CacheDamaged`
        To indicate a problem with the cache contents; the exception contains
        a ``.bad_files`` attribute containing a set of filenames to allow the
        user to use :func:`clear_download_cache` to remove the offending items.
    OSError, RuntimeError
        To indicate some problem with the cache structure. This may need a full
        :func:`clear_download_cache` to resolve, or may indicate some kind of
        misconfiguration.
    """
    bad_files = set()
    messages = set()
    dldir = _get_download_cache_loc(pkgname=pkgname)
    with os.scandir(dldir) as it:
        for entry in it:
            f = os.path.abspath(os.path.join(dldir, entry.name))
            if entry.name.startswith("rmtree-"):
                if f not in _tempfilestodel:
                    bad_files.add(f)
                    messages.add(f"Cache entry {entry.name} not scheduled for deletion")
            elif entry.is_dir():
                for sf in os.listdir(f):
                    if sf in ["url", "contents"]:
                        continue
                    sf = os.path.join(f, sf)
                    bad_files.add(sf)
                    messages.add(f"Unexpected file f{sf}")
                urlf = os.path.join(f, "url")
                url = None
                if not os.path.isfile(urlf):
                    bad_files.add(urlf)
                    messages.add(f"Problem with URL file f{urlf}")
                else:
                    url = get_file_contents(urlf, encoding="utf-8")
                    if not _is_url(url):
                        bad_files.add(f)
                        messages.add(f"Malformed URL: {url}")
                    else:
                        hashname = _url_to_dirname(url)
                        if entry.name != hashname:
                            bad_files.add(f)
                            messages.add(
                                f"URL hashes to {hashname} but is stored in"
                                f" {entry.name}"
                            )
                if not os.path.isfile(os.path.join(f, "contents")):
                    bad_files.add(f)
                    if url is None:
                        messages.add(f"Hash {entry.name} is missing contents")
                    else:
                        messages.add(
                            f"URL {url} with hash {entry.name} is missing contents"
                        )
            else:
                bad_files.add(f)
                messages.add(f"Left-over non-directory {f} in cache")
    if bad_files:
        raise CacheDamaged("\n".join(messages), bad_files=bad_files)


def _rmtree(path, replace=None):
    """More-atomic rmtree. Ignores missing directory."""
    with TemporaryDirectory(
        prefix="rmtree-", dir=os.path.dirname(os.path.abspath(path))
    ) as d:
        try:
            os.rename(path, os.path.join(d, "to-zap"))
        except FileNotFoundError:
            pass
        except PermissionError:
            warn(
                CacheMissingWarning(
                    f"Unable to remove directory {path} because a file in it "
                    "is in use and you are on Windows",
                    path,
                )
            )
            raise
        except OSError as e:
            if e.errno == errno.EXDEV:
                warn(e.strerror, AstropyWarning)
                shutil.move(path, os.path.join(d, "to-zap"))
            else:
                raise

        if replace is not None:
            try:
                os.rename(replace, path)
            except FileExistsError:
                # already there, fine
                pass
            except OSError as e:
                if e.errno == errno.ENOTEMPTY:
                    # already there, fine
                    pass
                elif e.errno == errno.EXDEV:
                    warn(e.strerror, AstropyWarning)
                    shutil.move(replace, path)
                else:
                    raise


def import_file_to_cache(
    url_key, filename, remove_original=False, pkgname="astropy", *, replace=True
):
    """Import the on-disk file specified by filename to the cache.

    The provided ``url_key`` will be the name used in the cache. The file
    should contain the contents of this URL, at least notionally (the URL may
    be temporarily or permanently unavailable). It is using ``url_key`` that
    users will request these contents from the cache. See :func:`download_file` for
    details.

    If ``url_key`` already exists in the cache, it will be updated to point to
    these imported contents, and its old contents will be deleted from the
    cache.

    Parameters
    ----------
    url_key : str
        The key to index the file under. This should probably be
        the URL where the file was located, though if you obtained
        it from a mirror you should use the URL of the primary
        location.
    filename : str
        The file whose contents you want to import.
    remove_original : bool
        Whether to remove the original file (``filename``) once import is
        complete.
    pkgname : `str`, optional
        The package name to use to locate the download cache. i.e. for
        ``pkgname='astropy'`` the default cache location is
        ``~/.astropy/cache``.
    replace : boolean, optional
        Whether or not to replace an existing object in the cache, if one exists.
        If replacement is not requested but the object exists, silently pass.
    """
    cache_dir = _get_download_cache_loc(pkgname=pkgname)
    cache_dirname = _url_to_dirname(url_key)
    local_dirname = os.path.join(cache_dir, cache_dirname)
    local_filename = os.path.join(local_dirname, "contents")
    with TemporaryDirectory(
        prefix="temp_dir", dir=cache_dir, ignore_cleanup_errors=True
    ) as temp_dir:
        temp_filename = os.path.join(temp_dir, "contents")
        # Make sure we're on the same filesystem
        # This will raise an exception if the url_key doesn't turn into a valid filename
        shutil.copy(filename, temp_filename)
        with open(os.path.join(temp_dir, "url"), "w", encoding="utf-8") as f:
            f.write(url_key)
        if replace:
            _rmtree(local_dirname, replace=temp_dir)
        else:
            try:
                os.rename(temp_dir, local_dirname)
            except FileExistsError:
                # already there, fine
                pass
            except OSError as e:
                if e.errno == errno.ENOTEMPTY:
                    # already there, fine
                    pass
                else:
                    raise
    if remove_original:
        os.remove(filename)
    return os.path.abspath(local_filename)


def get_cached_urls(pkgname="astropy"):
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
    return sorted(cache_contents(pkgname=pkgname).keys())


def cache_contents(pkgname="astropy"):
    """Obtain a dict mapping cached URLs to filenames.

    This dictionary is a read-only snapshot of the state of the cache when this
    function was called. If other processes are actively working with the
    cache, it is possible for them to delete files that are listed in this
    dictionary. Use with some caution if you are working on a system that is
    busy with many running astropy processes, although the same issues apply to
    most functions in this module.
    """
    r = {}
    try:
        dldir = _get_download_cache_loc(pkgname=pkgname)
    except OSError:
        return _NOTHING
    with os.scandir(dldir) as it:
        for entry in it:
            if entry.is_dir:
                url = get_file_contents(
                    os.path.join(dldir, entry.name, "url"), encoding="utf-8"
                )
                r[url] = os.path.abspath(os.path.join(dldir, entry.name, "contents"))
    return MappingProxyType(r)


def export_download_cache(
    filename_or_obj, urls=None, overwrite=False, pkgname="astropy"
):
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
    with zipfile.ZipFile(filename_or_obj, "w" if overwrite else "x") as z:
        for u in urls:
            fn = download_file(u, cache=True, sources=[], pkgname=pkgname)
            # Do not use os.path.join because ZIP files want
            # "/" on all platforms
            z_fn = urllib.parse.quote(u, safe="")
            z.write(fn, z_fn)


def import_download_cache(
    filename_or_obj, urls=None, update_cache=False, pkgname="astropy"
):
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
    with zipfile.ZipFile(filename_or_obj, "r") as z, TemporaryDirectory() as d:
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
                block = f_zip.read(conf.download_block_size)
                while block:
                    f_temp.write(block)
                    block = f_zip.read(conf.download_block_size)
            import_file_to_cache(
                url, f_temp_name, remove_original=True, pkgname=pkgname
            )
