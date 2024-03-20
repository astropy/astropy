# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
A "grab bag" of relatively small general-purpose utilities that don't have
a clear module/package to live in.
"""

import contextlib
import difflib
import inspect
import json
import locale
import os
import re
import signal
import sys
import threading
import traceback
import unicodedata
from contextlib import contextmanager

from astropy.utils import deprecated

__all__ = [
    "isiterable",
    "silence",
    "format_exception",
    "NumpyRNGContext",
    "find_api_page",
    "is_path_hidden",
    "walk_skip_hidden",
    "JsonCustomEncoder",
    "indent",
    "dtype_bytes_or_chars",
]

NOT_OVERWRITING_MSG = (
    "File {} already exists. If you mean to replace it "
    'then use the argument "overwrite=True".'
)
# A useful regex for tests.
_NOT_OVERWRITING_MSG_MATCH = (
    r"File .* already exists\. If you mean to "
    r"replace it then use the argument "
    r'"overwrite=True"\.'
)


def isiterable(obj):
    """Returns `True` if the given object is iterable."""
    try:
        iter(obj)
        return True
    except TypeError:
        return False


@deprecated(since="6.1", alternative="textwrap.indent()")
def indent(s, shift=1, width=4):
    """Indent a block of text.  The indentation is applied to each line."""
    indented = "\n".join(" " * (width * shift) + l if l else "" for l in s.splitlines())
    if s[-1] == "\n":
        indented += "\n"

    return indented


class _DummyFile:
    """A noop writeable object."""

    def write(self, s):
        pass


@contextlib.contextmanager
def silence():
    """A context manager that silences sys.stdout and sys.stderr."""
    old_stdout = sys.stdout
    old_stderr = sys.stderr
    sys.stdout = _DummyFile()
    sys.stderr = _DummyFile()
    yield
    sys.stdout = old_stdout
    sys.stderr = old_stderr


def format_exception(msg, *args, **kwargs):
    """Fill in information about the exception that occurred.

    Given an exception message string, uses new-style formatting arguments
    ``{filename}``, ``{lineno}``, ``{func}`` and/or ``{text}`` to fill in
    information about the exception that occurred.  For example:

        try:
            1/0
        except:
            raise ZeroDivisionError(
                format_except('A divide by zero occurred in {filename} at '
                              'line {lineno} of function {func}.'))

    Any additional positional or keyword arguments passed to this function are
    also used to format the message.

    .. note::
        This uses `sys.exc_info` to gather up the information needed to fill
        in the formatting arguments. Since `sys.exc_info` is not carried
        outside a handled exception, it's not wise to use this
        outside of an ``except`` clause - if it is, this will substitute
        '<unknown>' for the 4 formatting arguments.
    """
    tb = traceback.extract_tb(sys.exc_info()[2], limit=1)
    if len(tb) > 0:
        filename, lineno, func, text = tb[0]
    else:
        filename = lineno = func = text = "<unknown>"

    return msg.format(
        *args, filename=filename, lineno=lineno, func=func, text=text, **kwargs
    )


class NumpyRNGContext:
    """
    A context manager (for use with the ``with`` statement) that will seed the
    numpy random number generator (RNG) to a specific value, and then restore
    the RNG state back to whatever it was before.

    This is primarily intended for use in the astropy testing suit, but it
    may be useful in ensuring reproducibility of Monte Carlo simulations in a
    science context.

    Parameters
    ----------
    seed : int
        The value to use to seed the numpy RNG

    Examples
    --------
    A typical use case might be::

        with NumpyRNGContext(<some seed value you pick>):
            from numpy import random

            randarr = random.randn(100)
            ... run your test using `randarr` ...

        #Any code using numpy.random at this indent level will act just as it
        #would have if it had been before the with statement - e.g. whatever
        #the default seed is.


    """

    def __init__(self, seed):
        self.seed = seed

    def __enter__(self):
        from numpy import random

        self.startstate = random.get_state()
        random.seed(self.seed)

    def __exit__(self, exc_type, exc_value, traceback):
        from numpy import random

        random.set_state(self.startstate)


def find_api_page(obj, version=None, openinbrowser=True, timeout=None):
    """
    Determines the URL of the API page for the specified object, and
    optionally open that page in a web browser.

    .. note::
        You must be connected to the internet for this to function even if
        ``openinbrowser`` is `False`, unless you provide a local version of
        the documentation to ``version`` (e.g., ``file:///path/to/docs``).

    Parameters
    ----------
    obj
        The object to open the docs for or its fully-qualified name
        (as a str).
    version : str
        The doc version - either a version number like '0.1', 'dev' for
        the development/latest docs, or a URL to point to a specific
        location that should be the *base* of the documentation. Defaults to
        latest if you are on aren't on a release, otherwise, the version you
        are on.
    openinbrowser : bool
        If `True`, the `webbrowser` package will be used to open the doc
        page in a new web browser window.
    timeout : number, optional
        The number of seconds to wait before timing-out the query to
        the astropy documentation.  If not given, the default python
        stdlib timeout will be used.

    Returns
    -------
    url : str
        The loaded URL

    Raises
    ------
    ValueError
        If the documentation can't be found

    """
    import webbrowser
    from zlib import decompress

    from astropy.utils.data import get_readable_fileobj

    if (
        not isinstance(obj, str)
        and hasattr(obj, "__module__")
        and hasattr(obj, "__name__")
    ):
        obj = obj.__module__ + "." + obj.__name__
    elif inspect.ismodule(obj):
        obj = obj.__name__

    if version is None:
        from astropy import version

        if version.release:
            version = "v" + version.version
        else:
            version = "dev"

    if "://" in version:
        if version.endswith("index.html"):
            baseurl = version[:-10]
        elif version.endswith("/"):
            baseurl = version
        else:
            baseurl = version + "/"
    elif version == "dev" or version == "latest":
        baseurl = "http://devdocs.astropy.org/"
    else:
        baseurl = f"https://docs.astropy.org/en/{version}/"

    # Custom request headers; see
    # https://github.com/astropy/astropy/issues/8990
    url = baseurl + "objects.inv"
    headers = {"User-Agent": f"Astropy/{version}"}
    with get_readable_fileobj(
        url, encoding="binary", remote_timeout=timeout, http_headers=headers
    ) as uf:
        oiread = uf.read()

        # need to first read/remove the first four lines, which have info before
        # the compressed section with the actual object inventory
        idx = -1
        headerlines = []
        for _ in range(4):
            oldidx = idx
            idx = oiread.index(b"\n", oldidx + 1)
            headerlines.append(oiread[(oldidx + 1) : idx].decode("utf-8"))

        # intersphinx version line, project name, and project version
        ivers, proj, vers, compr = headerlines
        if "The remainder of this file is compressed using zlib" not in compr:
            raise ValueError(
                f"The file downloaded from {baseurl}objects.inv does not seem to be"
                "the usual Sphinx objects.inv format.  Maybe it "
                "has changed?"
            )

        compressed = oiread[(idx + 1) :]

    decompressed = decompress(compressed).decode("utf-8")

    resurl = None

    for l in decompressed.strip().splitlines():
        ls = l.split()
        name = ls[0]
        loc = ls[3]
        if loc.endswith("$"):
            loc = loc[:-1] + name

        if name == obj:
            resurl = baseurl + loc
            break

    if resurl is None:
        raise ValueError(f"Could not find the docs for the object {obj}")
    elif openinbrowser:
        webbrowser.open(resurl)

    return resurl


def signal_number_to_name(signum):
    """
    Given an OS signal number, returns a signal name.  If the signal
    number is unknown, returns ``'UNKNOWN'``.
    """
    # Since these numbers and names are platform specific, we use the
    # builtin signal module and build a reverse mapping.

    signal_to_name_map = {
        k: v for v, k in signal.__dict__.items() if v.startswith("SIG")
    }

    return signal_to_name_map.get(signum, "UNKNOWN")


# _has_hidden_attribute() can be deleted together with deprecated is_path_hidden() and
# walk_skip_hidden().
if sys.platform == "win32":
    import ctypes

    def _has_hidden_attribute(filepath):
        """
        Returns True if the given filepath has the hidden attribute on
        MS-Windows.  Based on a post here:
        https://stackoverflow.com/questions/284115/cross-platform-hidden-file-detection.
        """
        if isinstance(filepath, bytes):
            filepath = filepath.decode(sys.getfilesystemencoding())
        try:
            attrs = ctypes.windll.kernel32.GetFileAttributesW(filepath)
            result = bool(attrs & 2) and attrs != -1
        except AttributeError:
            result = False
        return result

else:

    def _has_hidden_attribute(filepath):
        return False


@deprecated(since="6.0")
def is_path_hidden(filepath):
    """
    Determines if a given file or directory is hidden.

    Parameters
    ----------
    filepath : str
        The path to a file or directory

    Returns
    -------
    hidden : bool
        Returns `True` if the file is hidden
    """
    name = os.path.basename(os.path.abspath(filepath))
    if isinstance(name, bytes):
        is_dotted = name.startswith(b".")
    else:
        is_dotted = name.startswith(".")
    return is_dotted or _has_hidden_attribute(filepath)


@deprecated(since="6.0")
def walk_skip_hidden(top, onerror=None, followlinks=False):
    """
    A wrapper for `os.walk` that skips hidden files and directories.

    This function does not have the parameter ``topdown`` from
    `os.walk`: the directories must always be recursed top-down when
    using this function.

    See Also
    --------
    os.walk : For a description of the parameters
    """
    for root, dirs, files in os.walk(
        top, topdown=True, onerror=onerror, followlinks=followlinks
    ):
        # These lists must be updated in-place so os.walk will skip
        # hidden directories
        dirs[:] = [d for d in dirs if not is_path_hidden(d)]
        files[:] = [f for f in files if not is_path_hidden(f)]
        yield root, dirs, files


class JsonCustomEncoder(json.JSONEncoder):
    """Support for data types that JSON default encoder
    does not do.

    This includes:

        * Numpy array or number
        * Complex number
        * Set
        * Bytes
        * astropy.UnitBase
        * astropy.Quantity

    Examples
    --------
    >>> import json
    >>> import numpy as np
    >>> from astropy.utils.misc import JsonCustomEncoder
    >>> json.dumps(np.arange(3), cls=JsonCustomEncoder)
    '[0, 1, 2]'

    """

    def default(self, obj):
        import numpy as np

        from astropy import units as u

        if isinstance(obj, u.Quantity):
            return dict(value=obj.value, unit=obj.unit.to_string())
        if isinstance(obj, (np.number, np.ndarray)):
            return obj.tolist()
        elif isinstance(obj, complex):
            return [obj.real, obj.imag]
        elif isinstance(obj, set):
            return list(obj)
        elif isinstance(obj, bytes):  # pragma: py3
            return obj.decode()
        elif isinstance(obj, (u.UnitBase, u.FunctionUnitBase)):
            if obj == u.dimensionless_unscaled:
                obj = "dimensionless_unit"
            else:
                return obj.to_string()

        return json.JSONEncoder.default(self, obj)


def strip_accents(s):
    """
    Remove accents from a Unicode string.

    This helps with matching "ångström" to "angstrom", for example.
    """
    return "".join(
        c for c in unicodedata.normalize("NFD", s) if unicodedata.category(c) != "Mn"
    )


def did_you_mean(s, candidates, n=3, cutoff=0.8, fix=None):
    """
    When a string isn't found in a set of candidates, we can be nice
    to provide a list of alternatives in the exception.  This
    convenience function helps to format that part of the exception.

    Parameters
    ----------
    s : str

    candidates : sequence of str or dict of str keys

    n : int
        The maximum number of results to include.  See
        `difflib.get_close_matches`.

    cutoff : float
        In the range [0, 1]. Possibilities that don't score at least
        that similar to word are ignored.  See
        `difflib.get_close_matches`.

    fix : callable
        A callable to modify the results after matching.  It should
        take a single string and return a sequence of strings
        containing the fixed matches.

    Returns
    -------
    message : str
        Returns the string "Did you mean X, Y, or Z?", or the empty
        string if no alternatives were found.
    """
    if isinstance(s, str):
        s = strip_accents(s)
    s_lower = s.lower()

    # Create a mapping from the lower case name to all capitalization
    # variants of that name.
    candidates_lower = {}
    for candidate in candidates:
        candidate_lower = candidate.lower()
        candidates_lower.setdefault(candidate_lower, [])
        candidates_lower[candidate_lower].append(candidate)

    # The heuristic here is to first try "singularizing" the word.  If
    # that doesn't match anything use difflib to find close matches in
    # original, lower and upper case.
    if s_lower.endswith("s") and s_lower[:-1] in candidates_lower:
        matches = [s_lower[:-1]]
    else:
        matches = difflib.get_close_matches(
            s_lower, candidates_lower, n=n, cutoff=cutoff
        )

    if len(matches):
        capitalized_matches = set()
        for match in matches:
            capitalized_matches.update(candidates_lower[match])
        matches = capitalized_matches

        if fix is not None:
            mapped_matches = []
            for match in matches:
                mapped_matches.extend(fix(match))
            matches = mapped_matches

        matches = list(set(matches))
        matches = sorted(matches)

        if len(matches) == 1:
            matches = matches[0]
        else:
            matches = ", ".join(matches[:-1]) + " or " + matches[-1]
        return f"Did you mean {matches}?"

    return ""


LOCALE_LOCK = threading.Lock()


@contextmanager
def _set_locale(name):
    """
    Context manager to temporarily set the locale to ``name``.

    An example is setting locale to "C" so that the C strtod()
    function will use "." as the decimal point to enable consistent
    numerical string parsing.

    Note that one cannot nest multiple _set_locale() context manager
    statements as this causes a threading lock.

    This code taken from https://stackoverflow.com/questions/18593661/how-do-i-strftime-a-date-object-in-a-different-locale.

    Parameters
    ----------
    name : str
        Locale name, e.g. "C" or "fr_FR".
    """
    name = str(name)

    with LOCALE_LOCK:
        saved = locale.setlocale(locale.LC_ALL)
        if saved == name:
            # Don't do anything if locale is already the requested locale
            yield
        else:
            try:
                locale.setlocale(locale.LC_ALL, name)
                yield
            finally:
                locale.setlocale(locale.LC_ALL, saved)


def dtype_bytes_or_chars(dtype):
    """
    Parse the number out of a dtype.str value like '<U5' or '<f8'.

    See #5819 for discussion on the need for this function for getting
    the number of characters corresponding to a string dtype.

    Parameters
    ----------
    dtype : numpy dtype object
        Input dtype

    Returns
    -------
    bytes_or_chars : int or None
        Bits (for numeric types) or characters (for string types)
    """
    match = re.search(r"(\d+)$", dtype.str)
    out = int(match.group(1)) if match else None
    return out


def _hungry_for(option):  # pragma: no cover
    """
    Open browser loaded with ``option`` options near you.

    *Disclaimers: Payments not included. Astropy is not
    responsible for any liability from using this function.*

    .. note:: Accuracy depends on your browser settings.

    """
    import webbrowser

    webbrowser.open(f"https://www.google.com/search?q={option}+near+me")


def pizza():  # pragma: no cover
    """``/pizza``."""
    _hungry_for("pizza")


def coffee(is_adam=False, is_brigitta=False):  # pragma: no cover
    """``/coffee``."""
    if is_adam and is_brigitta:
        raise ValueError("There can be only one!")
    if is_adam:
        option = "fresh+third+wave+coffee"
    elif is_brigitta:
        option = "decent+espresso"
    else:
        option = "coffee"
    _hungry_for(option)
