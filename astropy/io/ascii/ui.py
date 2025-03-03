# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""An extensible ASCII table reader and writer.

ui.py:
  Provides the main user functions for reading and writing tables.

:Copyright: Smithsonian Astrophysical Observatory (2010)
:Author: Tom Aldcroft (aldcroft@head.cfa.harvard.edu)
"""

import collections
import contextlib
import copy
import os
import pydoc
import re
import sys
import time
import typing
import warnings
from io import StringIO
from pathlib import Path

import numpy as np

from astropy.table import Table
from astropy.utils.data import get_readable_fileobj
from astropy.utils.decorators import format_doc
from astropy.utils.exceptions import AstropyWarning
from astropy.utils.misc import NOT_OVERWRITING_MSG

from . import (
    basic,
    cds,
    core,
    cparser,
    daophot,
    ecsv,
    fastbasic,
    fixedwidth,
    html,
    ipac,
    latex,
    mrt,
    rst,
    sextractor,
)
from .docs import READ_KWARG_TYPES, WRITE_KWARG_TYPES

_read_trace = []

# Default setting for guess parameter in read()
_GUESS = True


def _read_write_help(
    read_write: str, format: str | None = None, out: typing.IO | None = None
) -> None:
    """Helper function to output help documentation for read() or write().

    This uses the ``Table.read/write.help()`` functionality and modifies the output
    to look like the ``ascii.read/write()`` syntax.
    """
    help_func = getattr(Table, read_write).help
    format_parts = ["ascii"] + ([format] if format else [])
    help_str_io = StringIO()
    help_func(format=".".join(format_parts), out=help_str_io)
    help_str = help_str_io.getvalue()
    # Replace e.g. Table.read(format='ascii.ecsv') with ascii.read(format='ecsv').
    # Special case for format='ascii' which goes to ascii.read() or ascii.write().
    help_str = re.sub(
        r"Table\.(read|write)\(format='ascii\.(\w+)'\)",
        r"ascii.\1(format='\2')",
        help_str,
    )
    help_str = re.sub(
        r"Table\.(read|write)\(format='ascii'\)",
        r"ascii.\1()",
        help_str,
    )

    if out is None:
        pydoc.pager(help_str)
    else:
        out.write(help_str)


READ_WRITE_HELP = """Output help documentation for ``ascii.{read_write}()`` for the specified ``format``.

    By default the help output is printed to the console via ``pydoc.pager``.
    Instead one can supplied a file handle object as ``out`` and the output
    will be written to that handle.

    Parameters
    ----------
    format : str, None
        Format name, e.g. 'basic', 'ecsv' or 'html'
    out : None or file-like
        Output destination (default is stdout via a pager)
    """


@format_doc(READ_WRITE_HELP, read_write="read")
def read_help(format: str | None = None, out: typing.IO | None = None) -> None:
    _read_write_help("read", format=format, out=out)


@format_doc(READ_WRITE_HELP, read_write="write")
def write_help(format: str | None = None, out: typing.IO | None = None) -> None:
    _read_write_help("write", format=format, out=out)


def _probably_html(table, maxchars=100000):
    """
    Determine if ``table`` probably contains HTML content.  See PR #3693 and issue
    #3691 for context.
    """
    if not isinstance(table, str):
        try:
            # If table is an iterable (list of strings) then take the first
            # maxchars of these.  Make sure this is something with random
            # access to exclude a file-like object
            table[0]
            table[:1]
            size = 0
            for i, line in enumerate(table):
                size += len(line)
                if size > maxchars:
                    table = table[: i + 1]
                    break
            table = os.linesep.join(table)
        except Exception:
            pass

    if isinstance(table, str):
        # Look for signs of an HTML table in the first maxchars characters
        table = table[:maxchars]

        # URL ending in .htm or .html
        if re.match(
            r"( http[s]? | ftp | file ) :// .+ \.htm[l]?$",
            table,
            re.IGNORECASE | re.VERBOSE,
        ):
            return True

        # Filename ending in .htm or .html which exists
        if (
            re.search(r"\.htm[l]?$", table[-5:], re.IGNORECASE)
            and Path(table).expanduser().is_file()
        ):
            return True

        # Table starts with HTML document type declaration
        if re.match(r"\s* <! \s* DOCTYPE \s* HTML", table, re.IGNORECASE | re.VERBOSE):
            return True

        # Look for <TABLE .. >, <TR .. >, <TD .. > tag openers.
        if all(
            re.search(rf"< \s* {element} [^>]* >", table, re.IGNORECASE | re.VERBOSE)
            for element in ("table", "tr", "td")
        ):
            return True

    return False


def set_guess(guess):
    """
    Set the default value of the ``guess`` parameter for read().

    Parameters
    ----------
    guess : bool
        New default ``guess`` value (e.g., True or False)

    """
    global _GUESS
    _GUESS = guess


def get_reader(reader_cls=None, inputter_cls=None, outputter_cls=None, **kwargs):
    """
    Initialize a table reader allowing for common customizations.

    Most of the default behavior for various parameters is determined by the Reader
    class specified by ``reader_cls``.

    Parameters
    ----------
    reader_cls : `~astropy.io.ascii.BaseReader`
        Reader class. Default is :class:`Basic`.
    inputter_cls : `~astropy.io.ascii.BaseInputter`
        Inputter class
    outputter_cls : `~astropy.io.ascii.BaseOutputter`
        Outputter class
    delimiter : str
        Column delimiter string
    comment : str
        Regular expression defining a comment line in table
    quotechar : str
        One-character string to quote fields containing special characters
    header_start : int
        Line index for the header line not counting comment or blank lines. A line with
        only whitespace is considered blank.
    data_start : int
        Line index for the start of data not counting comment or blank lines. A line
        with only whitespace is considered blank.
    data_end : int
        Line index for the end of data not counting comment or blank lines. This value
        can be negative to count from the end.
    converters : dict
        Dict of converters.
    data_splitter_cls : `~astropy.io.ascii.BaseSplitter`
        Splitter class to split data columns.
    header_splitter_cls : `~astropy.io.ascii.BaseSplitter`
        Splitter class to split header columns.
    names : list
        List of names corresponding to each data column.
    include_names : list, optional
        List of names to include in output.
    exclude_names : list
        List of names to exclude from output (applied after ``include_names``).
    fill_values : tuple, list of tuple
        Specification of fill values for bad or missing table values.
    fill_include_names : list
        List of names to include in fill_values.
    fill_exclude_names : list
        List of names to exclude from fill_values (applied after
        ``fill_include_names``).

    Returns
    -------
    reader : `~astropy.io.ascii.BaseReader` subclass
        ASCII format reader instance
    """
    # This function is a light wrapper around core._get_reader to provide a
    # public interface with a default Reader.
    if reader_cls is None:
        # Default reader is Basic unless fast reader is forced
        fast_reader = _get_fast_reader_dict(kwargs)
        if fast_reader["enable"] == "force":
            reader_cls = fastbasic.FastBasic
        else:
            reader_cls = basic.Basic

    reader = core._get_reader(
        reader_cls, inputter_cls=inputter_cls, outputter_cls=outputter_cls, **kwargs
    )
    return reader


def _get_format_class(format):
    if format is None:
        return None
    if format in core.FORMAT_CLASSES:
        return core.FORMAT_CLASSES[format]
    raise ValueError(
        f"ASCII format {format!r} not in allowed list {sorted(core.FORMAT_CLASSES)}"
    )


def _get_fast_reader_dict(kwargs):
    """Convert 'fast_reader' key in kwargs into a dict if not already and make sure
    'enable' key is available.
    """
    fast_reader = copy.deepcopy(kwargs.get("fast_reader", True))
    if isinstance(fast_reader, dict):
        fast_reader.setdefault("enable", "force")
    else:
        fast_reader = {"enable": fast_reader}
    return fast_reader


def _validate_read_write_kwargs(read_write, **kwargs):
    """Validate types of keyword arg inputs to read() or write()."""

    def is_ducktype(val, cls):
        """Check if ``val`` is an instance of ``cls`` or "seems" like one:
        ``cls(val) == val`` does not raise and exception and is `True`. In
        this way you can pass in ``np.int16(2)`` and have that count as `int`.

        This has a special-case of ``cls`` being 'list-like', meaning it is
        an iterable but not a string.
        """
        if cls == "list-like":
            ok = not isinstance(val, str) and isinstance(val, collections.abc.Iterable)
        else:
            ok = isinstance(val, cls)
            if not ok:
                # See if ``val`` walks and quacks like a ``cls```.
                try:
                    new_val = cls(val)
                except Exception:
                    ok = False
                else:
                    ok = new_val == val
        return ok

    kwarg_types = READ_KWARG_TYPES if read_write == "read" else WRITE_KWARG_TYPES

    for arg, val in kwargs.items():
        # Kwarg type checking is opt-in, so kwargs not in the list are considered OK.
        # This reflects that some readers allow additional arguments that may not
        # be well-specified, e.g. ```__init__(self, **kwargs)`` is an option.
        if arg not in kwarg_types or val is None:
            continue

        # Single type or tuple of types for this arg (like isinstance())
        types = kwarg_types[arg]
        err_msg = (
            f"{read_write}() argument '{arg}' must be a "
            f"{types} object, got {type(val)} instead"
        )

        # Force `types` to be a tuple for the any() check below
        if not isinstance(types, tuple):
            types = (types,)

        if not any(is_ducktype(val, cls) for cls in types):
            raise TypeError(err_msg)


def _expand_user_if_path(argument):
    if isinstance(argument, (str, bytes, os.PathLike)):
        # For the `read()` method, a `str` input can be either a file path or
        # the table data itself. File names for io.ascii cannot have newlines
        # in them and io.ascii does not accept table data as `bytes`, so we can
        # attempt to detect data strings like this.
        is_str_data = isinstance(argument, str) and (
            "\n" in argument or "\r" in argument
        )
        if not is_str_data:
            # Remain conservative in expanding the presumed-path
            if (ex_user := Path(argument).expanduser()).exists():
                argument = str(ex_user)
    return argument


def read(table, guess=None, **kwargs):
    # This the final output from reading. Static analysis indicates the reading
    # logic (which is indeed complex) might not define `dat`, thus do so here.
    dat = None

    # Specifically block `reader_cls` kwarg, which will otherwise cause a confusing
    # exception later in the call to get_reader().
    if "reader_cls" in kwargs:
        raise TypeError("read() got an unexpected keyword argument 'reader_cls'")

    # Docstring defined below
    del _read_trace[:]

    # Downstream readers might munge kwargs
    kwargs = copy.deepcopy(kwargs)

    _validate_read_write_kwargs("read", **kwargs)

    # Convert 'fast_reader' key in kwargs into a dict if not already and make sure
    # 'enable' key is available.
    fast_reader = _get_fast_reader_dict(kwargs)
    kwargs["fast_reader"] = fast_reader

    if fast_reader["enable"] and fast_reader.get("chunk_size"):
        return _read_in_chunks(table, **kwargs)

    if "fill_values" not in kwargs:
        kwargs["fill_values"] = [("", "0")]

    # If an outputter_cls is supplied in kwargs that will take precedence.
    if (
        "outputter_cls" in kwargs
    ):  # user specified Outputter, not supported for fast reading
        fast_reader["enable"] = False

    format = kwargs.get("format")
    # Dictionary arguments are passed by reference per default and thus need
    # special protection:
    new_kwargs = copy.deepcopy(kwargs)
    kwargs["fast_reader"] = copy.deepcopy(fast_reader)

    reader_cls = _get_format_class(format)
    if reader_cls is not None:
        new_kwargs["reader_cls"] = reader_cls
        format = reader_cls._format_name

    # Remove format keyword if there, this is only allowed in read() not get_reader()
    if "format" in new_kwargs:
        del new_kwargs["format"]

    if guess is None:
        guess = _GUESS

    if guess:
        # If ``table`` is probably an HTML file then tell guess function to add
        # the HTML reader at the top of the guess list.  This is in response to
        # issue #3691 (and others) where libxml can segfault on a long non-HTML
        # file, thus prompting removal of the HTML reader from the default
        # guess list.
        new_kwargs["guess_html"] = _probably_html(table)

        # If `table` is a filename or readable file object then read in the
        # file now.  This prevents problems in Python 3 with the file object
        # getting closed or left at the file end.  See #3132, #3013, #3109,
        # #2001.  If a `readme` arg was passed that implies CDS format, in
        # which case the original `table` as the data filename must be left
        # intact.
        if "readme" not in new_kwargs:
            encoding = kwargs.get("encoding")
            try:
                table = _expand_user_if_path(table)
                with get_readable_fileobj(table, encoding=encoding) as fileobj:
                    table = fileobj.read()
            except ValueError:  # unreadable or invalid binary file
                raise
            except Exception:
                pass
            else:
                # Ensure that `table` has at least one \r or \n in it
                # so that the core.BaseInputter test of
                # ('\n' not in table and '\r' not in table)
                # will fail and so `table` cannot be interpreted there
                # as a filename.  See #4160.
                if not re.search(r"[\r\n]", table):
                    table = table + os.linesep

                # If the table got successfully read then look at the content
                # to see if is probably HTML, but only if it wasn't already
                # identified as HTML based on the filename.
                if not new_kwargs["guess_html"]:
                    new_kwargs["guess_html"] = _probably_html(table)

        # Get the table from guess in ``dat``.  If ``dat`` comes back as None
        # then there was just one set of kwargs in the guess list so fall
        # through below to the non-guess way so that any problems result in a
        # more useful traceback.
        dat = _guess(table, new_kwargs, format, fast_reader)
        if dat is None:
            guess = False

    if not guess:
        if format is None:
            reader = get_reader(**new_kwargs)
            format = reader._format_name

        table = _expand_user_if_path(table)

        # Try the fast reader version of `format` first if applicable.  Note that
        # if user specified a fast format (e.g. format='fast_basic') this test
        # will fail and the else-clause below will be used.
        if fast_reader["enable"] and f"fast_{format}" in core.FAST_CLASSES:
            fast_kwargs = copy.deepcopy(new_kwargs)
            fast_kwargs["reader_cls"] = core.FAST_CLASSES[f"fast_{format}"]
            fast_reader_rdr = get_reader(**fast_kwargs)
            try:
                dat = fast_reader_rdr.read(table)
                _read_trace.append(
                    {
                        "kwargs": copy.deepcopy(fast_kwargs),
                        "reader_cls": fast_reader_rdr.__class__,
                        "status": "Success with fast reader (no guessing)",
                    }
                )
            except (
                core.ParameterError,
                cparser.CParserError,
                UnicodeEncodeError,
            ) as err:
                # special testing value to avoid falling back on the slow reader
                if fast_reader["enable"] == "force":
                    raise core.InconsistentTableError(
                        f"fast reader {fast_reader_rdr.__class__} exception: {err}"
                    )
                # If the fast reader doesn't work, try the slow version
                reader = get_reader(**new_kwargs)
                dat = reader.read(table)
                _read_trace.append(
                    {
                        "kwargs": copy.deepcopy(new_kwargs),
                        "reader_cls": reader.__class__,
                        "status": (
                            "Success with slow reader after failing"
                            " with fast (no guessing)"
                        ),
                    }
                )
        else:
            reader = get_reader(**new_kwargs)
            dat = reader.read(table)
            _read_trace.append(
                {
                    "kwargs": copy.deepcopy(new_kwargs),
                    "reader_cls": reader.__class__,
                    "status": "Success with specified Reader class (no guessing)",
                }
            )

    # Static analysis (pyright) indicates `dat` might be left undefined, so just
    # to be sure define it at the beginning and check here.
    if dat is None:
        raise RuntimeError(
            "read() function failed due to code logic error, "
            "please report this bug on github"
        )

    return dat


read.__doc__ = core.READ_DOCSTRING
read.help = read_help


def _guess(table, read_kwargs, format, fast_reader):
    """
    Try to read the table using various sets of keyword args.  Start with the
    standard guess list and filter to make it unique and consistent with
    user-supplied read keyword args.  Finally, if none of those work then
    try the original user-supplied keyword args.

    Parameters
    ----------
    table : str, file-like, list
        Input table as a file name, file-like object, list of strings, or
        single newline-separated string.
    read_kwargs : dict
        Keyword arguments from user to be supplied to reader
    format : str
        Table format
    fast_reader : dict
        Options for the C engine fast reader.  See read() function for details.

    Returns
    -------
    dat : `~astropy.table.Table` or None
        Output table or None if only one guess format was available
    """
    # Keep a trace of all failed guesses kwarg
    failed_kwargs = []

    # Get an ordered list of read() keyword arg dicts that will be cycled
    # through in order to guess the format.
    full_list_guess = _get_guess_kwargs_list(read_kwargs)

    # If a fast version of the reader is available, try that before the slow version
    if (
        fast_reader["enable"]
        and format is not None
        and f"fast_{format}" in core.FAST_CLASSES
    ):
        fast_kwargs = copy.deepcopy(read_kwargs)
        fast_kwargs["reader_cls"] = core.FAST_CLASSES[f"fast_{format}"]
        full_list_guess = [fast_kwargs] + full_list_guess
    else:
        fast_kwargs = None

    # Filter the full guess list so that each entry is consistent with user kwarg inputs.
    # This also removes any duplicates from the list.
    filtered_guess_kwargs = []
    fast_reader = read_kwargs.get("fast_reader")

    for guess_kwargs in full_list_guess:
        # If user specified slow reader then skip all fast readers
        if (
            fast_reader["enable"] is False
            and guess_kwargs["reader_cls"] in core.FAST_CLASSES.values()
        ):
            _read_trace.append(
                {
                    "kwargs": copy.deepcopy(guess_kwargs),
                    "reader_cls": guess_kwargs["reader_cls"].__class__,
                    "status": "Disabled: reader only available in fast version",
                    "dt": f"{0.0:.3f} ms",
                }
            )
            continue

        # If user required a fast reader then skip all non-fast readers
        if (
            fast_reader["enable"] == "force"
            and guess_kwargs["reader_cls"] not in core.FAST_CLASSES.values()
        ):
            _read_trace.append(
                {
                    "kwargs": copy.deepcopy(guess_kwargs),
                    "reader_cls": guess_kwargs["reader_cls"].__class__,
                    "status": "Disabled: no fast version of reader available",
                    "dt": f"{0.0:.3f} ms",
                }
            )
            continue

        guess_kwargs_ok = True  # guess_kwargs are consistent with user_kwargs?
        for key, val in read_kwargs.items():
            # Do guess_kwargs.update(read_kwargs) except that if guess_args has
            # a conflicting key/val pair then skip this guess entirely.
            if key not in guess_kwargs:
                guess_kwargs[key] = copy.deepcopy(val)
            elif val != guess_kwargs[key] and guess_kwargs != fast_kwargs:
                guess_kwargs_ok = False
                break

        if not guess_kwargs_ok:
            # User-supplied kwarg is inconsistent with the guess-supplied kwarg, e.g.
            # user supplies delimiter="|" but the guess wants to try delimiter=" ",
            # so skip the guess entirely.
            continue

        # Add the guess_kwargs to filtered list only if it is not already there.
        if guess_kwargs not in filtered_guess_kwargs:
            filtered_guess_kwargs.append(guess_kwargs)

    # If there are not at least two formats to guess then return no table
    # (None) to indicate that guessing did not occur.  In that case the
    # non-guess read() will occur and any problems will result in a more useful
    # traceback.
    if len(filtered_guess_kwargs) <= 1:
        return None

    # Define whitelist of exceptions that are expected from readers when
    # processing invalid inputs.  Note that OSError must fall through here
    # so one cannot simply catch any exception.
    guess_exception_classes = (
        core.InconsistentTableError,
        ValueError,
        TypeError,
        AttributeError,
        core.OptionalTableImportError,
        core.ParameterError,
        cparser.CParserError,
    )

    # Determine whether we should limit the number of lines used in the guessing.
    # Note that this does not currently work for file objects, so we set this to
    # False if a file object was passed in.
    from astropy.io.ascii import conf  # avoid circular imports

    limit_lines = conf.guess_limit_lines if not hasattr(table, "read") else False

    # Don't limit the number of lines if there are fewer than this number of
    # lines in the table. In fact, we also don't limit the number of lines if
    # there are just above the number of lines compared to the limit, up to a
    # factor of 2, since it is fast to just go straight to the full table read.
    table_guess_subset = None

    if limit_lines:
        if isinstance(table, list):
            if len(table) > 2 * limit_lines:
                table_guess_subset = table[:limit_lines]
        else:
            # Now search for the position of the Nth line ending
            pos = -1
            for idx in range(limit_lines * 2):
                pos = table.find("\n", pos + 1)
                if pos == -1:
                    # Fewer than 2 * limit_lines line endings found so no guess subset.
                    break
                if idx == limit_lines - 1:
                    pos_limit = pos
            else:
                table_guess_subset = table[:pos_limit]

    # Now cycle through each possible reader and associated keyword arguments.
    # Try to read the table using those args, and if an exception occurs then
    # keep track of the failed guess and move on.
    for guess_kwargs in filtered_guess_kwargs:
        t0 = time.time()
        try:
            # If guessing will try all Readers then use strict req'ts on column names
            if "reader_cls" not in read_kwargs:
                guess_kwargs["strict_names"] = True

            reader = get_reader(**guess_kwargs)

            reader.guessing = True

            if table_guess_subset:
                # First try with subset of lines - if this fails we can skip this
                # format early. If it works, we still proceed to check with the
                # full table since we need to then return the read data.
                reader.read(table_guess_subset)

            dat = reader.read(table)
            _read_trace.append(
                {
                    "kwargs": copy.deepcopy(guess_kwargs),
                    "reader_cls": reader.__class__,
                    "status": "Success (guessing)",
                    "dt": f"{(time.time() - t0) * 1000:.3f} ms",
                }
            )
            return dat

        except guess_exception_classes as err:
            _read_trace.append(
                {
                    "kwargs": copy.deepcopy(guess_kwargs),
                    "status": f"{err.__class__.__name__}: {str(err)}",
                    "dt": f"{(time.time() - t0) * 1000:.3f} ms",
                }
            )
            failed_kwargs.append(guess_kwargs)

    # Failed all guesses, try the original read_kwargs without column requirements
    try:
        reader = get_reader(**read_kwargs)
        dat = reader.read(table)
        _read_trace.append(
            {
                "kwargs": copy.deepcopy(read_kwargs),
                "reader_cls": reader.__class__,
                "status": (
                    "Success with original kwargs without strict_names (guessing)"
                ),
            }
        )
        return dat

    except guess_exception_classes as err:
        _read_trace.append(
            {
                "kwargs": copy.deepcopy(read_kwargs),
                "status": f"{err.__class__.__name__}: {str(err)}",
            }
        )
        failed_kwargs.append(read_kwargs)
        lines = ["\nERROR: Unable to guess table format with the guesses listed below:"]
        for kwargs in failed_kwargs:
            sorted_keys = sorted(
                x for x in sorted(kwargs) if x not in ("reader_cls", "outputter_cls")
            )
            reader_repr = repr(kwargs.get("reader_cls", basic.Basic))
            keys_vals = ["reader_cls:" + re.search(r"\.(\w+)'>", reader_repr).group(1)]
            kwargs_sorted = ((key, kwargs[key]) for key in sorted_keys)
            keys_vals.extend([f"{key}: {val!r}" for key, val in kwargs_sorted])
            lines.append(" ".join(keys_vals))

    msg = [
        "",
        "************************************************************************",
        "** ERROR: Unable to guess table format with the guesses listed above. **",
        "**                                                                    **",
        "** To figure out why the table did not read, use guess=False and      **",
        "** fast_reader=False, along with any appropriate arguments to read(). **",
        "** In particular specify the format and any known attributes like the **",
        "** delimiter.                                                         **",
        "************************************************************************",
    ]
    lines.extend(msg)
    raise core.InconsistentTableError("\n".join(lines)) from None


def _get_guess_kwargs_list(read_kwargs):
    """Get the full list of reader keyword argument dicts.

    These are the basis for the format guessing process.
    The returned full list will then be:

    - Filtered to be consistent with user-supplied kwargs
    - Cleaned to have only unique entries
    - Used one by one to try reading the input table

    Note that the order of the guess list has been tuned over years of usage.
    Maintainers need to be very careful about any adjustments as the
    reasoning may not be immediately evident in all cases.

    This list can (and usually does) include duplicates.  This is a result
    of the order tuning, but these duplicates get removed later.

    Parameters
    ----------
    read_kwargs : dict
        User-supplied read keyword args

    Returns
    -------
    guess_kwargs_list : list
        List of read format keyword arg dicts
    """
    guess_kwargs_list = []

    # If the table is probably HTML based on some heuristics then start with the
    # HTML reader.
    if read_kwargs.pop("guess_html", None):
        guess_kwargs_list.append({"reader_cls": html.HTML})

    # Start with ECSV because an ECSV file will be read by Basic.  This format
    # has very specific header requirements and fails out quickly.
    guess_kwargs_list.append({"reader_cls": ecsv.Ecsv})

    # Now try readers that accept the user-supplied keyword arguments
    # (actually include all here - check for compatibility of arguments later).
    # FixedWidthTwoLine would also be read by Basic, so it needs to come first;
    # same for RST.
    for reader in (
        fixedwidth.FixedWidthTwoLine,
        rst.RST,
        fastbasic.FastBasic,
        basic.Basic,
        fastbasic.FastRdb,
        basic.Rdb,
        fastbasic.FastTab,
        basic.Tab,
        cds.Cds,
        mrt.Mrt,
        daophot.Daophot,
        sextractor.SExtractor,
        ipac.Ipac,
        latex.Latex,
        latex.AASTex,
    ):
        guess_kwargs_list.append({"reader_cls": reader})

    # Cycle through the basic-style readers using all combinations of delimiter
    # and quotechar.
    for reader_cls in (
        fastbasic.FastCommentedHeader,
        basic.CommentedHeader,
        fastbasic.FastBasic,
        basic.Basic,
        fastbasic.FastNoHeader,
        basic.NoHeader,
    ):
        for delimiter in ("|", ",", " ", r"\s"):
            for quotechar in ('"', "'"):
                guess_kwargs_list.append(
                    {
                        "reader_cls": reader_cls,
                        "delimiter": delimiter,
                        "quotechar": quotechar,
                    }
                )

    return guess_kwargs_list


def _read_in_chunks(table, **kwargs):
    """
    For fast_reader read the ``table`` in chunks and vstack to create
    a single table, OR return a generator of chunk tables.
    """
    fast_reader = kwargs["fast_reader"]
    chunk_size = fast_reader.pop("chunk_size")
    chunk_generator = fast_reader.pop("chunk_generator", False)

    tbl_chunks = _read_in_chunks_generator(table, chunk_size, **kwargs)
    if chunk_generator:
        return tbl_chunks

    tbl0 = next(tbl_chunks)
    masked = tbl0.masked

    # Numpy won't allow resizing the original so make a copy here.
    out_cols = {col.name: col.data.copy() for col in tbl0.itercols()}

    str_kinds = ("S", "U")
    for tbl in tbl_chunks:
        masked |= tbl.masked
        for name, col in tbl.columns.items():
            # Concatenate current column data and new column data

            # If one of the inputs is string-like and the other is not, then
            # convert the non-string to a string.  In a perfect world this would
            # be handled by numpy, but as of numpy 1.13 this results in a string
            # dtype that is too long (https://github.com/numpy/numpy/issues/10062).

            col1, col2 = out_cols[name], col.data
            if col1.dtype.kind in str_kinds and col2.dtype.kind not in str_kinds:
                col2 = np.array(col2.tolist(), dtype=col1.dtype.kind)
            elif col2.dtype.kind in str_kinds and col1.dtype.kind not in str_kinds:
                col1 = np.array(col1.tolist(), dtype=col2.dtype.kind)

            # Choose either masked or normal concatenation
            concatenate = np.ma.concatenate if masked else np.concatenate

            out_cols[name] = concatenate([col1, col2])

    # Make final table from numpy arrays, converting dict to list
    out_cols = [out_cols[name] for name in tbl0.colnames]
    out = tbl0.__class__(out_cols, names=tbl0.colnames, meta=tbl0.meta, copy=False)

    return out


def _read_in_chunks_generator(table, chunk_size, **kwargs):
    """
    For fast_reader read the ``table`` in chunks and return a generator
    of tables for each chunk.
    """

    @contextlib.contextmanager
    def passthrough_fileobj(fileobj, encoding=None):
        """Stub for get_readable_fileobj, which does not seem to work in Py3
        for input file-like object, see #6460.
        """
        yield fileobj

    # Set up to coerce `table` input into a readable file object by selecting
    # an appropriate function.

    # Convert table-as-string to a File object.  Finding a newline implies
    # that the string is not a filename.
    if isinstance(table, str) and ("\n" in table or "\r" in table):
        table = StringIO(table)
        fileobj_context = passthrough_fileobj
    elif hasattr(table, "read") and hasattr(table, "seek"):
        fileobj_context = passthrough_fileobj
    else:
        # string filename or pathlib
        fileobj_context = get_readable_fileobj

    # Set up for iterating over chunks
    kwargs["fast_reader"]["return_header_chars"] = True
    header = ""  # Table header (up to start of data)
    prev_chunk_chars = ""  # Chars from previous chunk after last newline
    first_chunk = True  # True for the first chunk, False afterward

    with fileobj_context(table, encoding=kwargs.get("encoding")) as fh:
        while True:
            chunk = fh.read(chunk_size)
            # Got fewer chars than requested, must be end of file
            final_chunk = len(chunk) < chunk_size

            # If this is the last chunk and there is only whitespace then break
            if final_chunk and not re.search(r"\S", chunk):
                break

            # Step backwards from last character in chunk and find first newline
            for idx in range(len(chunk) - 1, -1, -1):
                if final_chunk or chunk[idx] == "\n":
                    break
            else:
                raise ValueError("no newline found in chunk (chunk_size too small?)")

            # Stick on the header to the chunk part up to (and including) the
            # last newline.  Make sure the small strings are concatenated first.
            complete_chunk = (header + prev_chunk_chars) + chunk[: idx + 1]
            prev_chunk_chars = chunk[idx + 1 :]

            # Now read the chunk as a complete table
            tbl = read(complete_chunk, guess=False, **kwargs)

            # For the first chunk pop the meta key which contains the header
            # characters (everything up to the start of data) then fix kwargs
            # so it doesn't return that in meta any more.
            if first_chunk:
                header = tbl.meta.pop("__ascii_fast_reader_header_chars__")
                first_chunk = False

            yield tbl

            if final_chunk:
                break


extra_writer_pars = (
    "delimiter",
    "comment",
    "quotechar",
    "formats",
    "names",
    "include_names",
    "exclude_names",
    "strip_whitespace",
)


def get_writer(writer_cls=None, fast_writer=True, **kwargs):
    """
    Initialize a table writer allowing for common customizations.

    Most of the default behavior for various parameters is determined by the Writer
    class.

    Parameters
    ----------
    writer_cls : ``writer_cls``
        Writer class. Defaults to :class:`Basic`.
    delimiter : str
        Column delimiter string
    comment : str
        String defining a comment line in table
    quotechar : str
        One-character string to quote fields containing special characters
    formats : dict
        Dictionary of format specifiers or formatting functions
    strip_whitespace : bool
        Strip surrounding whitespace from column values.
    names : list
        List of names corresponding to each data column
    include_names : list
        List of names to include in output.
    exclude_names : list
        List of names to exclude from output (applied after ``include_names``)
    fast_writer : bool
        Whether to use the fast Cython writer.

    Returns
    -------
    writer : `~astropy.io.ascii.BaseReader` subclass
        ASCII format writer instance
    """
    if writer_cls is None:
        writer_cls = basic.Basic
    if "strip_whitespace" not in kwargs:
        kwargs["strip_whitespace"] = True
    writer = core._get_writer(writer_cls, fast_writer, **kwargs)

    # Handle the corner case of wanting to disable writing table comments for the
    # commented_header format.  This format *requires* a string for `write_comment`
    # because that is used for the header column row, so it is not possible to
    # set the input `comment` to None.  Without adding a new keyword or assuming
    # a default comment character, there is no other option but to tell user to
    # simply remove the meta['comments'].
    if isinstance(
        writer, (basic.CommentedHeader, fastbasic.FastCommentedHeader)
    ) and not isinstance(kwargs.get("comment", ""), str):
        raise ValueError(
            "for the commented_header writer you must supply a string\n"
            "value for the `comment` keyword.  In order to disable writing\n"
            "table comments use `del t.meta['comments']` prior to writing."
        )

    return writer


def write(
    table,
    output=None,
    format=None,
    fast_writer=True,
    *,
    overwrite=False,
    **kwargs,
):
    # Docstring inserted below

    # Specifically block the legacy `writer_cls` kwarg, which will otherwise cause a confusing
    # exception later in the call to get_writer().
    if "writer_cls" in kwargs:
        raise TypeError("write() got an unexpected keyword argument 'writer_cls'")

    _validate_read_write_kwargs(
        "write", format=format, fast_writer=fast_writer, overwrite=overwrite, **kwargs
    )

    if isinstance(output, (str, bytes, os.PathLike)):
        output = os.path.expanduser(output)  # noqa: PTH111
        if not overwrite and os.path.lexists(output):
            raise OSError(NOT_OVERWRITING_MSG.format(output))

    if output is None:
        output = sys.stdout

    # Ensure that `table` is a Table subclass.
    names = kwargs.get("names")
    if isinstance(table, Table):
        # While we are only going to read data from columns, we may need to
        # to adjust info attributes such as format, so we make a shallow copy.
        table = table.__class__(table, names=names, copy=False)
    else:
        # Otherwise, create a table from the input.
        table = Table(table, names=names, copy=False)

    table0 = table[:0].copy()
    core._apply_include_exclude_names(
        table0,
        kwargs.get("names"),
        kwargs.get("include_names"),
        kwargs.get("exclude_names"),
    )
    diff_format_with_names = set(kwargs.get("formats", [])) - set(table0.colnames)

    if diff_format_with_names:
        warnings.warn(
            (
                f"The key(s) {diff_format_with_names} specified in the formats "
                "argument do not match a column name."
            ),
            AstropyWarning,
        )

    if table.has_mixin_columns:
        fast_writer = False

    writer_cls = _get_format_class(format)
    writer = get_writer(writer_cls=writer_cls, fast_writer=fast_writer, **kwargs)
    if writer._format_name in core.FAST_CLASSES:
        writer.write(table, output)
        return

    lines = writer.write(table)

    # Write the lines to output
    outstr = os.linesep.join(lines)
    if not hasattr(output, "write"):
        # NOTE: we need to specify newline='', otherwise the default
        # behavior is for Python to translate \r\n (which we write because
        # of os.linesep) into \r\r\n. Specifying newline='' disables any
        # auto-translation.
        with open(output, "w", newline="") as output:
            output.write(outstr)
            output.write(os.linesep)
    else:
        output.write(outstr)
        output.write(os.linesep)


write.__doc__ = core.WRITE_DOCSTRING
write.help = write_help


def get_read_trace():
    """
    Return a traceback of the attempted read formats for the last call to
    `~astropy.io.ascii.read` where guessing was enabled.  This is primarily for
    debugging.

    The return value is a list of dicts, where each dict includes the keyword
    args ``kwargs`` used in the read call and the returned ``status``.

    Returns
    -------
    trace : list of dict
        Ordered list of format guesses and status
    """
    return copy.deepcopy(_read_trace)
