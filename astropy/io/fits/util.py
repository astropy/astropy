# Licensed under a 3-clause BSD style license - see PYFITS.rst

import gzip
import io
import itertools
import mmap
import operator
import os
import platform
import signal
import sys
import tempfile
import textwrap
import threading
import warnings
import weakref
from contextlib import contextmanager, suppress
from functools import wraps

import numpy as np
from packaging.version import Version

from astropy.utils import data
from astropy.utils.exceptions import AstropyUserWarning

path_like = (str, bytes, os.PathLike)

cmp = lambda a, b: (a > b) - (a < b)

all_integer_types = (int, np.integer)


class NotifierMixin:
    """
    Mixin class that provides services by which objects can register
    listeners to changes on that object.

    All methods provided by this class are underscored, since this is intended
    for internal use to communicate between classes in a generic way, and is
    not machinery that should be exposed to users of the classes involved.

    Use the ``_add_listener`` method to register a listener on an instance of
    the notifier.  This registers the listener with a weak reference, so if
    no other references to the listener exist it is automatically dropped from
    the list and does not need to be manually removed.

    Call the ``_notify`` method on the notifier to update all listeners
    upon changes.  ``_notify('change_type', *args, **kwargs)`` results
    in calling ``listener._update_change_type(*args, **kwargs)`` on all
    listeners subscribed to that notifier.

    If a particular listener does not have the appropriate update method
    it is ignored.

    Examples
    --------
    >>> class Widget(NotifierMixin):
    ...     state = 1
    ...     def __init__(self, name):
    ...         self.name = name
    ...     def update_state(self):
    ...         self.state += 1
    ...         self._notify('widget_state_changed', self)
    ...
    >>> class WidgetListener:
    ...     def _update_widget_state_changed(self, widget):
    ...         print('Widget {0} changed state to {1}'.format(
    ...             widget.name, widget.state))
    ...
    >>> widget = Widget('fred')
    >>> listener = WidgetListener()
    >>> widget._add_listener(listener)
    >>> widget.update_state()
    Widget fred changed state to 2
    """

    _listeners = None

    def _add_listener(self, listener):
        """
        Add an object to the list of listeners to notify of changes to this
        object.  This adds a weakref to the list of listeners that is
        removed from the listeners list when the listener has no other
        references to it.
        """
        if self._listeners is None:
            self._listeners = weakref.WeakValueDictionary()

        self._listeners[id(listener)] = listener

    def _remove_listener(self, listener):
        """
        Removes the specified listener from the listeners list.  This relies
        on object identity (i.e. the ``is`` operator).
        """
        if self._listeners is None:
            return

        with suppress(KeyError):
            del self._listeners[id(listener)]

    def _notify(self, notification, *args, **kwargs):
        """
        Notify all listeners of some particular state change by calling their
        ``_update_<notification>`` method with the given ``*args`` and
        ``**kwargs``.

        The notification does not by default include the object that actually
        changed (``self``), but it certainly may if required.
        """
        if self._listeners is None:
            return

        method_name = f"_update_{notification}"
        for listener in self._listeners.valuerefs():
            # Use valuerefs instead of itervaluerefs; see
            # https://github.com/astropy/astropy/issues/4015
            listener = listener()  # dereference weakref
            if listener is None:
                continue

            if hasattr(listener, method_name):
                method = getattr(listener, method_name)
                if callable(method):
                    method(*args, **kwargs)

    def __getstate__(self):
        """
        Exclude listeners when saving the listener's state, since they may be
        ephemeral.
        """
        # TODO: This hasn't come up often, but if anyone needs to pickle HDU
        # objects it will be necessary when HDU objects' states are restored to
        # re-register themselves as listeners on their new column instances.
        try:
            state = super().__getstate__()
        except AttributeError:
            # Chances are the super object doesn't have a getstate
            state = self.__dict__.copy()

        state["_listeners"] = None
        return state


def first(iterable):
    """
    Returns the first item returned by iterating over an iterable object.

    Examples
    --------
    >>> a = [1, 2, 3]
    >>> first(a)
    1
    """
    return next(iter(iterable))


def itersubclasses(cls, _seen=None):
    """
    Generator over all subclasses of a given class, in depth first order.

    >>> class A: pass
    >>> class B(A): pass
    >>> class C(A): pass
    >>> class D(B,C): pass
    >>> class E(D): pass
    >>>
    >>> for cls in itersubclasses(A):
    ...     print(cls.__name__)
    B
    D
    E
    C
    >>> # get ALL classes currently defined
    >>> [cls.__name__ for cls in itersubclasses(object)]
    [...'tuple', ...'type', ...]

    From http://code.activestate.com/recipes/576949/
    """
    if _seen is None:
        _seen = set()
    try:
        subs = cls.__subclasses__()
    except TypeError:  # fails only when cls is type
        subs = cls.__subclasses__(cls)
    for sub in sorted(subs, key=operator.attrgetter("__name__")):
        if sub not in _seen:
            _seen.add(sub)
            yield sub
            for sub in itersubclasses(sub, _seen):
                yield sub


def ignore_sigint(func):
    """
    This decorator registers a custom SIGINT handler to catch and ignore SIGINT
    until the wrapped function is completed.
    """

    @wraps(func)
    def wrapped(*args, **kwargs):
        # Get the name of the current thread and determine if this is a single
        # threaded application
        curr_thread = threading.current_thread()
        single_thread = (
            threading.active_count() == 1 and curr_thread.name == "MainThread"
        )

        class SigintHandler:
            def __init__(self):
                self.sigint_received = False

            def __call__(self, signum, frame):
                warnings.warn(
                    f"KeyboardInterrupt ignored until {func.__name__} is complete!",
                    AstropyUserWarning,
                )
                self.sigint_received = True

        sigint_handler = SigintHandler()

        # Define new signal interput handler
        if single_thread:
            # Install new handler
            old_handler = signal.signal(signal.SIGINT, sigint_handler)

        try:
            func(*args, **kwargs)
        finally:
            if single_thread:
                if old_handler is not None:
                    signal.signal(signal.SIGINT, old_handler)
                else:
                    signal.signal(signal.SIGINT, signal.SIG_DFL)

                if sigint_handler.sigint_received:
                    raise KeyboardInterrupt

    return wrapped


if sys.version_info[:2] >= (3, 10):
    from itertools import pairwise
else:

    def pairwise(iterable):
        """Return the items of an iterable paired with its next item.

        Ex: s -> (s0,s1), (s1,s2), (s2,s3), ....
        """
        a, b = itertools.tee(iterable)
        for _ in b:
            # Just a little trick to advance b without having to catch
            # StopIter if b happens to be empty
            break
        return zip(a, b)


def encode_ascii(s):
    if isinstance(s, str):
        return s.encode("ascii")
    elif isinstance(s, np.ndarray) and issubclass(s.dtype.type, np.str_):
        ns = np.char.encode(s, "ascii").view(type(s))
        if ns.dtype.itemsize != s.dtype.itemsize / 4:
            ns = ns.astype((np.bytes_, s.dtype.itemsize / 4))
        return ns
    elif isinstance(s, np.ndarray) and not issubclass(s.dtype.type, np.bytes_):
        raise TypeError("string operation on non-string array")
    return s


def decode_ascii(s):
    if isinstance(s, bytes):
        try:
            return s.decode("ascii")
        except UnicodeDecodeError:
            warnings.warn(
                "non-ASCII characters are present in the FITS "
                'file header and have been replaced by "?" characters',
                AstropyUserWarning,
            )
            s = s.decode("ascii", errors="replace")
            return s.replace("\ufffd", "?")
    elif isinstance(s, np.ndarray) and issubclass(s.dtype.type, np.bytes_):
        # np.char.encode/decode annoyingly don't preserve the type of the
        # array, hence the view() call
        # It also doesn't necessarily preserve widths of the strings,
        # hence the astype()
        if s.size == 0:
            # Numpy apparently also has a bug that if a string array is
            # empty calling np.char.decode on it returns an empty float64
            # array : https://github.com/numpy/numpy/issues/13156
            dt = s.dtype.str.replace("S", "U")
            ns = np.array([], dtype=dt).view(type(s))
        else:
            ns = np.char.decode(s, "ascii").view(type(s))
        if ns.dtype.itemsize / 4 != s.dtype.itemsize:
            ns = ns.astype((np.str_, s.dtype.itemsize))
        return ns
    elif isinstance(s, np.ndarray) and not issubclass(s.dtype.type, np.str_):
        # Don't silently pass through on non-string arrays; we don't want
        # to hide errors where things that are not stringy are attempting
        # to be decoded
        raise TypeError("string operation on non-string array")
    return s


def isreadable(f):
    """
    Returns True if the file-like object can be read from.  This is a common-
    sense approximation of io.IOBase.readable.
    """
    if hasattr(f, "readable"):
        return f.readable()

    if hasattr(f, "closed") and f.closed:
        # This mimics the behavior of io.IOBase.readable
        raise ValueError("I/O operation on closed file")

    if not hasattr(f, "read"):
        return False

    if hasattr(f, "mode") and not any(c in f.mode for c in "r+"):
        return False

    # Not closed, has a 'read()' method, and either has no known mode or a
    # readable mode--should be good enough to assume 'readable'
    return True


def iswritable(f):
    """
    Returns True if the file-like object can be written to.  This is a common-
    sense approximation of io.IOBase.writable.
    """
    if hasattr(f, "writable"):
        return f.writable()

    if hasattr(f, "closed") and f.closed:
        # This mimics the behavior of io.IOBase.writable
        raise ValueError("I/O operation on closed file")

    if not hasattr(f, "write"):
        return False

    if hasattr(f, "mode") and not any(c in f.mode for c in "wa+"):
        return False

    # Note closed, has a 'write()' method, and either has no known mode or a
    # mode that supports writing--should be good enough to assume 'writable'
    return True


def isfile(f):
    """
    Returns True if the given object represents an OS-level file (that is,
    ``isinstance(f, file)``).

    On Python 3 this also returns True if the given object is higher level
    wrapper on top of a FileIO object, such as a TextIOWrapper.
    """
    if isinstance(f, io.FileIO):
        return True
    elif hasattr(f, "buffer"):
        return isfile(f.buffer)
    elif hasattr(f, "raw"):
        return isfile(f.raw)
    return False


def fileobj_name(f):
    """
    Returns the 'name' of file-like object *f*, if it has anything that could be
    called its name.  Otherwise f's class or type is returned.  If f is a
    string f itself is returned.
    """
    if isinstance(f, (str, bytes)):
        return f
    elif isinstance(f, gzip.GzipFile):
        # The .name attribute on GzipFiles does not always represent the name
        # of the file being read/written--it can also represent the original
        # name of the file being compressed
        # See the documentation at
        # https://docs.python.org/3/library/gzip.html#gzip.GzipFile
        # As such, for gzip files only return the name of the underlying
        # fileobj, if it exists
        return fileobj_name(f.fileobj)
    elif hasattr(f, "name"):
        return f.name
    elif hasattr(f, "filename"):
        return f.filename
    elif hasattr(f, "__class__"):
        return str(f.__class__)
    else:
        return str(type(f))


def fileobj_closed(f):
    """
    Returns True if the given file-like object is closed or if *f* is a string
    (and assumed to be a pathname).

    Returns False for all other types of objects, under the assumption that
    they are file-like objects with no sense of a 'closed' state.
    """
    if isinstance(f, path_like):
        return True

    if hasattr(f, "closed"):
        return f.closed
    elif hasattr(f, "fileobj") and hasattr(f.fileobj, "closed"):
        return f.fileobj.closed
    elif hasattr(f, "fp") and hasattr(f.fp, "closed"):
        return f.fp.closed
    else:
        return False


def fileobj_mode(f):
    """
    Returns the 'mode' string of a file-like object if such a thing exists.
    Otherwise returns None.
    """
    # Go from most to least specific--for example gzip objects have a 'mode'
    # attribute, but it's not analogous to the file.mode attribute

    # gzip.GzipFile -like
    if hasattr(f, "fileobj") and hasattr(f.fileobj, "mode"):
        fileobj = f.fileobj

    # astropy.io.fits._File -like, doesn't need additional checks because it's
    # already validated
    elif hasattr(f, "fileobj_mode"):
        return f.fileobj_mode

    # PIL-Image -like investigate the fp (filebuffer)
    elif hasattr(f, "fp") and hasattr(f.fp, "mode"):
        fileobj = f.fp

    # FILEIO -like (normal open(...)), keep as is.
    elif hasattr(f, "mode"):
        fileobj = f

    # Doesn't look like a file-like object, for example strings, urls or paths.
    else:
        return None

    return _fileobj_normalize_mode(fileobj)


def _fileobj_normalize_mode(f):
    """Takes care of some corner cases in Python where the mode string
    is either oddly formatted or does not truly represent the file mode.
    """
    mode = f.mode

    # Special case: Gzip modes:
    if isinstance(f, gzip.GzipFile):
        # GzipFiles can be either readonly or writeonly
        if mode == gzip.READ:
            return "rb"
        elif mode == gzip.WRITE:
            return "wb"
        else:
            return None  # This shouldn't happen?

    # Sometimes Python can produce modes like 'r+b' which will be normalized
    # here to 'rb+'
    if "+" in mode:
        mode = mode.replace("+", "")
        mode += "+"

    return mode


def fileobj_is_binary(f):
    """
    Returns True if the give file or file-like object has a file open in binary
    mode.  When in doubt, returns True by default.
    """
    # This is kind of a hack for this to work correctly with _File objects,
    # which, for the time being, are *always* binary
    if hasattr(f, "binary"):
        return f.binary

    if isinstance(f, io.TextIOBase):
        return False

    mode = fileobj_mode(f)
    if mode:
        return "b" in mode
    else:
        return True


def translate(s, table, deletechars):
    if deletechars:
        table = table.copy()
        for c in deletechars:
            table[ord(c)] = None
    return s.translate(table)


def fill(text, width, **kwargs):
    """
    Like :func:`textwrap.wrap` but preserves existing paragraphs which
    :func:`textwrap.wrap` does not otherwise handle well.  Also handles section
    headers.
    """
    paragraphs = text.split("\n\n")

    def maybe_fill(t):
        if all(len(line) < width for line in t.splitlines()):
            return t
        else:
            return textwrap.fill(t, width, **kwargs)

    return "\n\n".join(maybe_fill(p) for p in paragraphs)


# On MacOS X 10.8 and earlier, there is a bug that causes numpy.fromfile to
# fail when reading over 2Gb of data. If we detect these versions of MacOS X,
# we can instead read the data in chunks. To avoid performance penalties at
# import time, we defer the setting of this global variable until the first
# time it is needed.
CHUNKED_FROMFILE = None


def _array_from_file(infile, dtype, count):
    """Create a numpy array from a file or a file-like object."""
    if isfile(infile):
        global CHUNKED_FROMFILE
        if CHUNKED_FROMFILE is None:
            if sys.platform == "darwin" and Version(platform.mac_ver()[0]) < Version(
                "10.9"
            ):
                CHUNKED_FROMFILE = True
            else:
                CHUNKED_FROMFILE = False

        if CHUNKED_FROMFILE:
            chunk_size = int(1024**3 / dtype.itemsize)  # 1Gb to be safe
            if count < chunk_size:
                return np.fromfile(infile, dtype=dtype, count=count)
            else:
                array = np.empty(count, dtype=dtype)
                for beg in range(0, count, chunk_size):
                    end = min(count, beg + chunk_size)
                    array[beg:end] = np.fromfile(infile, dtype=dtype, count=end - beg)
                return array
        else:
            return np.fromfile(infile, dtype=dtype, count=count)
    else:
        # treat as file-like object with "read" method; this includes gzip file
        # objects, because numpy.fromfile just reads the compressed bytes from
        # their underlying file object, instead of the decompressed bytes
        read_size = np.dtype(dtype).itemsize * count
        s = infile.read(read_size)
        array = np.ndarray(buffer=s, dtype=dtype, shape=(count,))
        # copy is needed because np.frombuffer returns a read-only view of the
        # underlying buffer
        array = array.copy()
        return array


_OSX_WRITE_LIMIT = (2**32) - 1
_WIN_WRITE_LIMIT = (2**31) - 1


def _array_to_file(arr, outfile):
    """
    Write a numpy array to a file or a file-like object.

    Parameters
    ----------
    arr : ndarray
        The Numpy array to write.
    outfile : file-like
        A file-like object such as a Python file object, an `io.BytesIO`, or
        anything else with a ``write`` method.  The file object must support
        the buffer interface in its ``write``.

    If writing directly to an on-disk file this delegates directly to
    `ndarray.tofile`.  Otherwise a slower Python implementation is used.
    """
    try:
        seekable = outfile.seekable()
    except AttributeError:
        seekable = False

    if isfile(outfile) and seekable:
        write = lambda a, f: a.tofile(f)
    else:
        write = _array_to_file_like

    # Implements a workaround for a bug deep in OSX's stdlib file writing
    # functions; on 64-bit OSX it is not possible to correctly write a number
    # of bytes greater than 2 ** 32 and divisible by 4096 (or possibly 8192--
    # whatever the default blocksize for the filesystem is).
    # This issue should have a workaround in Numpy too, but hasn't been
    # implemented there yet: https://github.com/astropy/astropy/issues/839
    #
    # Apparently Windows has its own fwrite bug:
    # https://github.com/numpy/numpy/issues/2256

    if (
        sys.platform == "darwin"
        and arr.nbytes >= _OSX_WRITE_LIMIT + 1
        and arr.nbytes % 4096 == 0
    ):
        # chunksize is a count of elements in the array, not bytes
        chunksize = _OSX_WRITE_LIMIT // arr.itemsize
    elif sys.platform.startswith("win"):
        chunksize = _WIN_WRITE_LIMIT // arr.itemsize
    else:
        # Just pass the whole array to the write routine
        return write(arr, outfile)

    # Write one chunk at a time for systems whose fwrite chokes on large
    # writes.
    idx = 0
    arr = arr.view(np.ndarray).flatten()
    while idx < arr.nbytes:
        write(arr[idx : idx + chunksize], outfile)
        idx += chunksize


def _array_to_file_like(arr, fileobj):
    """
    Write a `~numpy.ndarray` to a file-like object (which is not supported by
    `numpy.ndarray.tofile`).
    """
    # If the array is empty, we can simply take a shortcut and return since
    # there is nothing to write.
    if len(arr) == 0:
        return

    if arr.flags.contiguous:
        # It suffices to just pass the underlying buffer directly to the
        # fileobj's write (assuming it supports the buffer interface). If
        # it does not have the buffer interface, a TypeError should be returned
        # in which case we can fall back to the other methods.

        try:
            fileobj.write(arr.data)
        except TypeError:
            pass
        else:
            return

    if hasattr(np, "nditer"):
        # nditer version for non-contiguous arrays
        for item in np.nditer(arr, order="C"):
            fileobj.write(item.tobytes())
    else:
        # Slower version for Numpy versions without nditer;
        # The problem with flatiter is it doesn't preserve the original
        # byteorder
        byteorder = arr.dtype.byteorder
        if (sys.byteorder == "little" and byteorder == ">") or (
            sys.byteorder == "big" and byteorder == "<"
        ):
            for item in arr.flat:
                fileobj.write(item.byteswap().tobytes())
        else:
            for item in arr.flat:
                fileobj.write(item.tobytes())


def _write_string(f, s):
    """
    Write a string to a file, encoding to ASCII if the file is open in binary
    mode, or decoding if the file is open in text mode.
    """
    # Assume if the file object doesn't have a specific mode, that the mode is
    # binary
    binmode = fileobj_is_binary(f)

    if binmode and isinstance(s, str):
        s = encode_ascii(s)
    elif not binmode and not isinstance(f, str):
        s = decode_ascii(s)

    f.write(s)


def _convert_array(array, dtype):
    """
    Converts an array to a new dtype--if the itemsize of the new dtype is
    the same as the old dtype and both types are not numeric, a view is
    returned.  Otherwise a new array must be created.
    """
    if array.dtype == dtype:
        return array
    elif array.dtype.itemsize == dtype.itemsize and not (
        np.issubdtype(array.dtype, np.number) and np.issubdtype(dtype, np.number)
    ):
        # Includes a special case when both dtypes are at least numeric to
        # account for old Trac ticket 218 (now inaccessible).
        return array.view(dtype)
    else:
        return array.astype(dtype)


def _pseudo_zero(dtype):
    """
    Given a numpy dtype, finds its "zero" point, which is exactly in the
    middle of its range.
    """
    # special case for int8
    if dtype.kind == "i" and dtype.itemsize == 1:
        return -128

    assert dtype.kind == "u"
    return 1 << (dtype.itemsize * 8 - 1)


def _is_pseudo_integer(dtype):
    return (dtype.kind == "u" and dtype.itemsize >= 2) or (
        dtype.kind == "i" and dtype.itemsize == 1
    )


def _is_int(val):
    return isinstance(val, all_integer_types)


def _str_to_num(val):
    """Converts a given string to either an int or a float if necessary."""
    try:
        num = int(val)
    except ValueError:
        # If this fails then an exception should be raised anyways
        num = float(val)
    return num


def _words_group(s, width):
    """
    Split a long string into parts where each part is no longer than ``strlen``
    and no word is cut into two pieces.  But if there are any single words
    which are longer than ``strlen``, then they will be split in the middle of
    the word.
    """
    words = []
    slen = len(s)

    # appending one blank at the end always ensures that the "last" blank
    # is beyond the end of the string
    arr = np.frombuffer(s.encode("utf8") + b" ", dtype="S1")

    # locations of the blanks
    blank_loc = np.nonzero(arr == b" ")[0]
    offset = 0
    xoffset = 0

    while True:
        try:
            loc = np.nonzero(blank_loc >= width + offset)[0][0]
        except IndexError:
            loc = len(blank_loc)

        if loc > 0:
            offset = blank_loc[loc - 1] + 1
        else:
            offset = -1

        # check for one word longer than strlen, break in the middle
        if offset <= xoffset:
            offset = min(xoffset + width, slen)

        # collect the pieces in a list
        words.append(s[xoffset:offset])
        if offset >= slen:
            break
        xoffset = offset

    return words


def _tmp_name(input):
    """
    Create a temporary file name which should not already exist.  Use the
    directory of the input file as the base name of the mkstemp() output.
    """
    if input is not None:
        input = os.path.dirname(input)
    f, fn = tempfile.mkstemp(dir=input)
    os.close(f)
    return fn


def _get_array_mmap(array):
    """
    If the array has an mmap.mmap at base of its base chain, return the mmap
    object; otherwise return None.
    """
    if isinstance(array, mmap.mmap):
        return array

    base = array
    while hasattr(base, "base") and base.base is not None:
        if isinstance(base.base, mmap.mmap):
            return base.base
        base = base.base


@contextmanager
def _free_space_check(hdulist, dirname=None):
    try:
        yield
    except OSError as exc:
        error_message = ""
        if not isinstance(hdulist, list):
            hdulist = [hdulist]
        if dirname is None:
            dirname = os.path.dirname(hdulist._file.name)
        if os.path.isdir(dirname):
            free_space = data.get_free_space_in_dir(dirname)
            hdulist_size = sum(hdu.size for hdu in hdulist)
            if free_space < hdulist_size:
                error_message = (
                    "Not enough space on disk: requested {}, available {}. ".format(
                        hdulist_size, free_space
                    )
                )

        for hdu in hdulist:
            hdu._close()

        raise OSError(error_message + str(exc))


def _extract_number(value, default):
    """
    Attempts to extract an integer number from the given value. If the
    extraction fails, the value of the 'default' argument is returned.
    """
    try:
        # The _str_to_num method converts the value to string/float
        # so we need to perform one additional conversion to int on top
        return int(_str_to_num(value))
    except (TypeError, ValueError):
        return default


def get_testdata_filepath(filename):
    """
    Return a string representing the path to the file requested from the
    io.fits test data set.

    .. versionadded:: 2.0.3

    Parameters
    ----------
    filename : str
        The filename of the test data file.

    Returns
    -------
    filepath : str
        The path to the requested file.
    """
    return data.get_pkg_data_filename(f"io/fits/tests/data/{filename}", "astropy")


def _rstrip_inplace(array):
    """
    Performs an in-place rstrip operation on string arrays. This is necessary
    since the built-in `np.char.rstrip` in Numpy does not perform an in-place
    calculation.
    """
    # The following implementation convert the string to unsigned integers of
    # the right length. Trailing spaces (which are represented as 32) are then
    # converted to null characters (represented as zeros). To avoid creating
    # large temporary mask arrays, we loop over chunks (attempting to do that
    # on a 1-D version of the array; large memory may still be needed in the
    # unlikely case that a string array has small first dimension and cannot
    # be represented as a contiguous 1-D array in memory).

    dt = array.dtype

    if dt.kind not in "SU":
        raise TypeError("This function can only be used on string arrays")
    # View the array as appropriate integers. The last dimension will
    # equal the number of characters in each string.
    bpc = 1 if dt.kind == "S" else 4
    dt_int = f"({dt.itemsize // bpc},){dt.byteorder}u{bpc}"
    b = array.view(dt_int, np.ndarray)
    # For optimal speed, work in chunks of the internal ufunc buffer size.
    bufsize = np.getbufsize()
    # Attempt to have the strings as a 1-D array to give the chunk known size.
    # Note: the code will work if this fails; the chunks will just be larger.
    if b.ndim > 2:
        try:
            b.shape = -1, b.shape[-1]
        except AttributeError:  # can occur for non-contiguous arrays
            pass
    for j in range(0, b.shape[0], bufsize):
        c = b[j : j + bufsize]
        # Mask which will tell whether we're in a sequence of trailing spaces.
        mask = np.ones(c.shape[:-1], dtype=bool)
        # Loop over the characters in the strings, in reverse order. We process
        # the i-th character of all strings in the chunk at the same time. If
        # the character is 32, this corresponds to a space, and we then change
        # this to 0. We then construct a new mask to find rows where the
        # i-th character is 0 (null) and the i-1-th is 32 (space) and repeat.
        for i in range(-1, -c.shape[-1], -1):
            mask &= c[..., i] == 32
            c[..., i][mask] = 0
            mask = c[..., i] == 0

    return array


def _is_dask_array(data):
    """Check whether data is a dask array.

    We avoid importing dask unless it is likely it is a dask array,
    so that non-dask code is not slowed down.
    """
    if not hasattr(data, "compute"):
        return False

    try:
        from dask.array import Array
    except ImportError:
        # If we cannot import dask, surely this cannot be a
        # dask array!
        return False
    else:
        return isinstance(data, Array)
