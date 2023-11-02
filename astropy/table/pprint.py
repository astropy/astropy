# Licensed under a 3-clause BSD style license - see LICENSE.rst

import fnmatch
import os
import re
import sys

import numpy as np

from astropy import log
from astropy.utils.console import Getch, color_print, conf, terminal_size
from astropy.utils.data_info import dtype_info_name

__all__ = []


def default_format_func(format_, val):
    if isinstance(val, bytes):
        return val.decode("utf-8", errors="replace")
    else:
        return str(val)


# The first three functions are helpers for _auto_format_func


def _use_str_for_masked_values(format_func):
    """Wrap format function to trap masked values.

    String format functions and most user functions will not be able to deal
    with masked values, so we wrap them to ensure they are passed to str().
    """
    return lambda format_, val: (
        str(val) if val is np.ma.masked else format_func(format_, val)
    )


def _possible_string_format_functions(format_):
    """Iterate through possible string-derived format functions.

    A string can either be a format specifier for the format built-in,
    a new-style format string, or an old-style format string.
    """
    yield lambda format_, val: format(val, format_)
    yield lambda format_, val: format_.format(val)
    yield lambda format_, val: format_ % val
    yield lambda format_, val: format_.format(**{k: val[k] for k in val.dtype.names})


def get_auto_format_func(
    col=None, possible_string_format_functions=_possible_string_format_functions
):
    """
    Return a wrapped ``auto_format_func`` function which is used in
    formatting table columns.  This is primarily an internal function but
    gets used directly in other parts of astropy, e.g. `astropy.io.ascii`.

    Parameters
    ----------
    col_name : object, optional
        Hashable object to identify column like id or name. Default is None.

    possible_string_format_functions : func, optional
        Function that yields possible string formatting functions
        (defaults to internal function to do this).

    Returns
    -------
    Wrapped ``auto_format_func`` function
    """

    def _auto_format_func(format_, val):
        """Format ``val`` according to ``format_`` for a plain format specifier,
        old- or new-style format strings, or using a user supplied function.
        More importantly, determine and cache (in _format_funcs) a function
        that will do this subsequently.  In this way this complicated logic is
        only done for the first value.

        Returns the formatted value.
        """
        if format_ is None:
            return default_format_func(format_, val)

        if format_ in col.info._format_funcs:
            return col.info._format_funcs[format_](format_, val)

        if callable(format_):
            format_func = lambda format_, val: format_(val)
            try:
                out = format_func(format_, val)
                if not isinstance(out, str):
                    raise ValueError(
                        f"Format function for value {val} returned {type(val)} "
                        "instead of string type"
                    )
            except Exception as err:
                # For a masked element, the format function call likely failed
                # to handle it.  Just return the string representation for now,
                # and retry when a non-masked value comes along.
                if val is np.ma.masked:
                    return str(val)

                raise ValueError(f"Format function for value {val} failed.") from err
            # If the user-supplied function handles formatting masked elements, use
            # it directly.  Otherwise, wrap it in a function that traps them.
            try:
                format_func(format_, np.ma.masked)
            except Exception:
                format_func = _use_str_for_masked_values(format_func)
        else:
            # For a masked element, we cannot set string-based format functions yet,
            # as all tests below will fail.  Just return the string representation
            # of masked for now, and retry when a non-masked value comes along.
            if val is np.ma.masked:
                return str(val)

            for format_func in possible_string_format_functions(format_):
                try:
                    # Does this string format method work?
                    out = format_func(format_, val)
                    # Require that the format statement actually did something.
                    if out == format_:
                        raise ValueError("the format passed in did nothing.")
                except Exception:
                    continue
                else:
                    break
            else:
                # None of the possible string functions passed muster.
                raise ValueError(
                    f"unable to parse format string {format_} for its column."
                )

            # String-based format functions will fail on masked elements;
            # wrap them in a function that traps them.
            format_func = _use_str_for_masked_values(format_func)

        col.info._format_funcs[format_] = format_func
        return out

    return _auto_format_func


def _get_pprint_include_names(table):
    """Get the set of names to show in pprint from the table pprint_include_names
    and pprint_exclude_names attributes.

    These may be fnmatch unix-style globs.
    """

    def get_matches(name_globs, default):
        match_names = set()
        if name_globs:  # For None or () use the default
            for name in table.colnames:
                for name_glob in name_globs:
                    if fnmatch.fnmatch(name, name_glob):
                        match_names.add(name)
                        break
        else:
            match_names.update(default)
        return match_names

    include_names = get_matches(table.pprint_include_names(), table.colnames)
    exclude_names = get_matches(table.pprint_exclude_names(), [])

    return include_names - exclude_names


class TableFormatter:
    @staticmethod
    def _get_pprint_size(max_lines=None, max_width=None):
        """Get the output size (number of lines and character width) for Column and
        Table pformat/pprint methods.

        If no value of ``max_lines`` is supplied then the height of the
        screen terminal is used to set ``max_lines``.  If the terminal
        height cannot be determined then the default will be determined
        using the ``astropy.table.conf.max_lines`` configuration item. If a
        negative value of ``max_lines`` is supplied then there is no line
        limit applied.

        The same applies for max_width except the configuration item is
        ``astropy.table.conf.max_width``.

        Parameters
        ----------
        max_lines : int or None
            Maximum lines of output (header + data rows)

        max_width : int or None
            Maximum width (characters) output

        Returns
        -------
        max_lines, max_width : int

        """
        # Declare to keep static type checker happy.
        lines = None
        width = None

        if max_lines is None:
            max_lines = conf.max_lines

        if max_width is None:
            max_width = conf.max_width

        if max_lines is None or max_width is None:
            lines, width = terminal_size()

        if max_lines is None:
            max_lines = lines
        elif max_lines < 0:
            max_lines = sys.maxsize
        if max_lines < 8:
            max_lines = 8

        if max_width is None:
            max_width = width
        elif max_width < 0:
            max_width = sys.maxsize
        if max_width < 10:
            max_width = 10

        return max_lines, max_width

    def _pformat_col(
        self,
        col,
        max_lines=None,
        show_name=True,
        show_unit=None,
        show_dtype=False,
        show_length=None,
        html=False,
        align=None,
    ):
        """Return a list of formatted string representation of column values.

        Parameters
        ----------
        max_lines : int
            Maximum lines of output (header + data rows)

        show_name : bool
            Include column name. Default is True.

        show_unit : bool
            Include a header row for unit.  Default is to show a row
            for units only if one or more columns has a defined value
            for the unit.

        show_dtype : bool
            Include column dtype. Default is False.

        show_length : bool
            Include column length at end.  Default is to show this only
            if the column is not shown completely.

        html : bool
            Output column as HTML

        align : str
            Left/right alignment of columns. Default is '>' (right) for all
            columns. Other allowed values are '<', '^', and '0=' for left,
            centered, and 0-padded, respectively.

        Returns
        -------
        lines : list
            List of lines with formatted column values

        outs : dict
            Dict which is used to pass back additional values
            defined within the iterator.

        """
        if show_unit is None:
            show_unit = col.info.unit is not None

        outs = {}  # Some values from _pformat_col_iter iterator that are needed here
        col_strs_iter = self._pformat_col_iter(
            col,
            max_lines,
            show_name=show_name,
            show_unit=show_unit,
            show_dtype=show_dtype,
            show_length=show_length,
            outs=outs,
        )

        # Replace tab and newline with text representations so they display nicely.
        # Newline in particular is a problem in a multicolumn table.
        col_strs = [
            val.replace("\t", "\\t").replace("\n", "\\n") for val in col_strs_iter
        ]
        if len(col_strs) > 0:
            col_width = max(len(x) for x in col_strs)

        if html:
            from astropy.utils.xml.writer import xml_escape

            n_header = outs["n_header"]
            for i, col_str in enumerate(col_strs):
                # _pformat_col output has a header line '----' which is not needed here
                if i == n_header - 1:
                    continue
                td = "th" if i < n_header else "td"
                val = f"<{td}>{xml_escape(col_str.strip())}</{td}>"
                row = "<tr>" + val + "</tr>"
                if i < n_header:
                    row = "<thead>" + row + "</thead>"
                col_strs[i] = row

            if n_header > 0:
                # Get rid of '---' header line
                col_strs.pop(n_header - 1)
            col_strs.insert(0, "<table>")
            col_strs.append("</table>")

        # Now bring all the column string values to the same fixed width
        else:
            col_width = max(len(x) for x in col_strs) if col_strs else 1

            # Center line header content and generate dashed headerline
            for i in outs["i_centers"]:
                col_strs[i] = col_strs[i].center(col_width)
            if outs["i_dashes"] is not None:
                col_strs[outs["i_dashes"]] = "-" * col_width

            # Format columns according to alignment.  `align` arg has precedent, otherwise
            # use `col.format` if it starts as a legal alignment string.  If neither applies
            # then right justify.
            re_fill_align = re.compile(r"(?P<fill>.?)(?P<align>[<^>=])")
            match = None
            if align:
                # If there is an align specified then it must match
                match = re_fill_align.match(align)
                if not match:
                    raise ValueError(
                        "column align must be one of '<', '^', '>', or '='"
                    )
            elif isinstance(col.info.format, str):
                # col.info.format need not match, in which case rjust gets used
                match = re_fill_align.match(col.info.format)

            if match:
                fill_char = match.group("fill")
                align_char = match.group("align")
                if align_char == "=":
                    if fill_char != "0":
                        raise ValueError("fill character must be '0' for '=' align")
                    # str.zfill gets used which does not take fill char arg
                    fill_char = ""
            else:
                fill_char = ""
                align_char = ">"

            justify_methods = {"<": "ljust", "^": "center", ">": "rjust", "=": "zfill"}
            justify_method = justify_methods[align_char]
            justify_args = (col_width, fill_char) if fill_char else (col_width,)

            for i, col_str in enumerate(col_strs):
                col_strs[i] = getattr(col_str, justify_method)(*justify_args)

        if outs["show_length"]:
            col_strs.append(f"Length = {len(col)} rows")

        return col_strs, outs

    def _name_and_structure(self, name, dtype, sep=" "):
        """Format a column name, including a possible structure.

        Normally, just returns the name, but if it has a structured dtype,
        will add the parts in between square brackets.  E.g.,
        "name [f0, f1]" or "name [f0[sf0, sf1], f1]".
        """
        if dtype is None or dtype.names is None:
            return name

        structure = ", ".join(
            [
                self._name_and_structure(name, dt, sep="")
                for name, (dt, _) in dtype.fields.items()
            ]
        )
        return f"{name}{sep}[{structure}]"

    def _pformat_col_iter(
        self,
        col,
        max_lines,
        show_name,
        show_unit,
        outs,
        show_dtype=False,
        show_length=None,
    ):
        """Iterator which yields formatted string representation of column values.

        Parameters
        ----------
        max_lines : int
            Maximum lines of output (header + data rows)

        show_name : bool
            Include column name. Default is True.

        show_unit : bool
            Include a header row for unit.  Default is to show a row
            for units only if one or more columns has a defined value
            for the unit.

        outs : dict
            Must be a dict which is used to pass back additional values
            defined within the iterator.

        show_dtype : bool
            Include column dtype. Default is False.

        show_length : bool
            Include column length at end.  Default is to show this only
            if the column is not shown completely.
        """
        max_lines, _ = self._get_pprint_size(max_lines, -1)
        dtype = getattr(col, "dtype", None)
        multidims = getattr(col, "shape", [0])[1:]
        if multidims:
            multidim0 = tuple(0 for n in multidims)
            multidim1 = tuple(n - 1 for n in multidims)
            multidims_all_ones = np.prod(multidims) == 1
            multidims_has_zero = 0 in multidims

        i_dashes = None
        i_centers = []  # Line indexes where content should be centered
        n_header = 0
        if show_name:
            i_centers.append(n_header)
            # Get column name (or 'None' if not set)
            col_name = str(col.info.name)
            n_header += 1
            yield self._name_and_structure(col_name, dtype)
        if show_unit:
            i_centers.append(n_header)
            n_header += 1
            yield str(col.info.unit or "")
        if show_dtype:
            i_centers.append(n_header)
            n_header += 1
            if dtype is not None:
                col_dtype = dtype_info_name((dtype, multidims))
            else:
                col_dtype = col.__class__.__qualname__ or "object"
            yield col_dtype
        if show_unit or show_name or show_dtype:
            i_dashes = n_header
            n_header += 1
            yield "---"

        max_lines -= n_header
        n_print2 = max_lines // 2
        n_rows = len(col)

        # This block of code is responsible for producing the function that
        # will format values for this column.  The ``format_func`` function
        # takes two args (col_format, val) and returns the string-formatted
        # version.  Some points to understand:
        #
        # - col_format could itself be the formatting function, so it will
        #    actually end up being called with itself as the first arg.  In
        #    this case the function is expected to ignore its first arg.
        #
        # - auto_format_func is a function that gets called on the first
        #    column value that is being formatted.  It then determines an
        #    appropriate formatting function given the actual value to be
        #    formatted.  This might be deterministic or it might involve
        #    try/except.  The latter allows for different string formatting
        #    options like %f or {:5.3f}.  When auto_format_func is called it:

        #    1. Caches the function in the _format_funcs dict so for subsequent
        #       values the right function is called right away.
        #    2. Returns the formatted value.
        #
        # - possible_string_format_functions is a function that yields a
        #    succession of functions that might successfully format the
        #    value.  There is a default, but Mixin methods can override this.
        #    See Quantity for an example.
        #
        # - get_auto_format_func() returns a wrapped version of auto_format_func
        #    with the column id and possible_string_format_functions as
        #    enclosed variables.
        col_format = col.info.format or getattr(col.info, "default_format", None)
        pssf = (
            getattr(col.info, "possible_string_format_functions", None)
            or _possible_string_format_functions
        )
        auto_format_func = get_auto_format_func(col, pssf)
        format_func = col.info._format_funcs.get(col_format, auto_format_func)

        if len(col) > max_lines:
            if show_length is None:
                show_length = True
            i0 = n_print2 - (1 if show_length else 0)
            i1 = n_rows - n_print2 - max_lines % 2
            indices = np.concatenate(
                [np.arange(0, i0 + 1), np.arange(i1 + 1, len(col))]
            )
        else:
            i0 = -1
            indices = np.arange(len(col))

        def format_col_str(idx):
            if multidims:
                # Prevents columns like Column(data=[[(1,)],[(2,)]], name='a')
                # with shape (n,1,...,1) from being printed as if there was
                # more than one element in a row
                if multidims_all_ones:
                    return format_func(col_format, col[(idx,) + multidim0])
                elif multidims_has_zero:
                    # Any zero dimension means there is no data to print
                    return ""
                else:
                    left = format_func(col_format, col[(idx,) + multidim0])
                    right = format_func(col_format, col[(idx,) + multidim1])
                    return f"{left} .. {right}"
            else:
                return format_func(col_format, col[idx])

        # Add formatted values if within bounds allowed by max_lines
        for idx in indices:
            if idx == i0:
                yield "..."
            else:
                try:
                    yield format_col_str(idx)
                except ValueError:
                    raise ValueError(
                        'Unable to parse format string "{}" for entry "{}" '
                        'in column "{}"'.format(col_format, col[idx], col.info.name)
                    )

        outs["show_length"] = show_length
        outs["n_header"] = n_header
        outs["i_centers"] = i_centers
        outs["i_dashes"] = i_dashes

    def _pformat_table(
        self,
        table,
        max_lines=None,
        max_width=None,
        show_name=True,
        show_unit=None,
        show_dtype=False,
        html=False,
        tableid=None,
        tableclass=None,
        align=None,
    ):
        """Return a list of lines for the formatted string representation of
        the table.

        Parameters
        ----------
        max_lines : int or None
            Maximum number of rows to output

        max_width : int or None
            Maximum character width of output

        show_name : bool
            Include a header row for column names. Default is True.

        show_unit : bool
            Include a header row for unit.  Default is to show a row
            for units only if one or more columns has a defined value
            for the unit.

        show_dtype : bool
            Include a header row for column dtypes. Default is to False.

        html : bool
            Format the output as an HTML table. Default is False.

        tableid : str or None
            An ID tag for the table; only used if html is set.  Default is
            "table{id}", where id is the unique integer id of the table object,
            id(table)

        tableclass : str or list of str or None
            CSS classes for the table; only used if html is set.  Default is
            none

        align : str or list or tuple
            Left/right alignment of columns. Default is '>' (right) for all
            columns. Other allowed values are '<', '^', and '0=' for left,
            centered, and 0-padded, respectively. A list of strings can be
            provided for alignment of tables with multiple columns.

        Returns
        -------
        rows : list
            Formatted table as a list of strings

        outs : dict
            Dict which is used to pass back additional values
            defined within the iterator.

        """
        # "Print" all the values into temporary lists by column for subsequent
        # use and to determine the width
        max_lines, max_width = self._get_pprint_size(max_lines, max_width)

        if show_unit is None:
            show_unit = any(col.info.unit for col in table.columns.values())

        # Coerce align into a correctly-sized list of alignments (if possible)
        n_cols = len(table.columns)
        if align is None or isinstance(align, str):
            align = [align] * n_cols

        elif isinstance(align, (list, tuple)):
            if len(align) != n_cols:
                raise ValueError(
                    f"got {len(align)} alignment values instead of "
                    f"the number of columns ({n_cols})"
                )
        else:
            raise TypeError(
                f"align keyword must be str or list or tuple (got {type(align)})"
            )

        # Process column visibility from table pprint_include_names and
        # pprint_exclude_names attributes and get the set of columns to show.
        pprint_include_names = _get_pprint_include_names(table)

        cols = []
        outs = None  # Initialize so static type checker is happy
        for align_, col in zip(align, table.columns.values()):
            if col.info.name not in pprint_include_names:
                continue

            lines, outs = self._pformat_col(
                col,
                max_lines,
                show_name=show_name,
                show_unit=show_unit,
                show_dtype=show_dtype,
                align=align_,
            )
            if outs["show_length"]:
                lines = lines[:-1]
            cols.append(lines)

        if not cols:
            return ["<No columns>"], {"show_length": False}

        # Use the values for the last column since they are all the same
        n_header = outs["n_header"]

        n_rows = len(cols[0])

        def outwidth(cols):
            return sum(len(c[0]) for c in cols) + len(cols) - 1

        dots_col = ["..."] * n_rows
        middle = len(cols) // 2
        while outwidth(cols) > max_width:
            if len(cols) == 1:
                break
            if len(cols) == 2:
                cols[1] = dots_col
                break
            if cols[middle] is dots_col:
                cols.pop(middle)
                middle = len(cols) // 2
            cols[middle] = dots_col

        # Now "print" the (already-stringified) column values into a
        # row-oriented list.
        rows = []
        if html:
            from astropy.utils.xml.writer import xml_escape

            if tableid is None:
                tableid = f"table{id(table)}"

            if tableclass is not None:
                if isinstance(tableclass, list):
                    tableclass = " ".join(tableclass)
                rows.append(f'<table id="{tableid}" class="{tableclass}">')
            else:
                rows.append(f'<table id="{tableid}">')

            for i in range(n_rows):
                # _pformat_col output has a header line '----' which is not needed here
                if i == n_header - 1:
                    continue
                td = "th" if i < n_header else "td"
                vals = (f"<{td}>{xml_escape(col[i].strip())}</{td}>" for col in cols)
                row = "<tr>" + "".join(vals) + "</tr>"
                if i < n_header:
                    row = "<thead>" + row + "</thead>"
                rows.append(row)
            rows.append("</table>")
        else:
            for i in range(n_rows):
                row = " ".join(col[i] for col in cols)
                rows.append(row)

        return rows, outs

    def _more_tabcol(
        self,
        tabcol,
        max_lines=None,
        max_width=None,
        show_name=True,
        show_unit=None,
        show_dtype=False,
    ):
        """Interactive "more" of a table or column.

        Parameters
        ----------
        max_lines : int or None
            Maximum number of rows to output

        max_width : int or None
            Maximum character width of output

        show_name : bool
            Include a header row for column names. Default is True.

        show_unit : bool
            Include a header row for unit.  Default is to show a row
            for units only if one or more columns has a defined value
            for the unit.

        show_dtype : bool
            Include a header row for column dtypes. Default is False.
        """
        allowed_keys = "f br<>qhpn"

        # Count the header lines
        n_header = 0
        if show_name:
            n_header += 1
        if show_unit:
            n_header += 1
        if show_dtype:
            n_header += 1
        if show_name or show_unit or show_dtype:
            n_header += 1

        # Set up kwargs for pformat call.  Only Table gets max_width.
        kwargs = dict(
            max_lines=-1,
            show_name=show_name,
            show_unit=show_unit,
            show_dtype=show_dtype,
        )
        if hasattr(tabcol, "columns"):  # tabcol is a table
            kwargs["max_width"] = max_width

        # If max_lines is None (=> query screen size) then increase by 2.
        # This is because get_pprint_size leaves 6 extra lines so that in
        # ipython you normally see the last input line.
        max_lines1, max_width = self._get_pprint_size(max_lines, max_width)
        if max_lines is None:
            max_lines1 += 2
        delta_lines = max_lines1 - n_header

        # Set up a function to get a single character on any platform
        inkey = Getch()

        i0 = 0  # First table/column row to show
        showlines = True
        while True:
            i1 = i0 + delta_lines  # Last table/col row to show
            if showlines:  # Don't always show the table (e.g. after help)
                try:
                    os.system("cls" if os.name == "nt" else "clear")
                except Exception:
                    pass  # No worries if clear screen call fails
                lines = tabcol[i0:i1].pformat(**kwargs)
                colors = (
                    "red" if i < n_header else "default" for i in range(len(lines))
                )
                for color, line in zip(colors, lines):
                    color_print(line, color)
            showlines = True
            print()
            print("-- f, <space>, b, r, p, n, <, >, q h (help) --", end=" ")
            # Get a valid key
            while True:
                try:
                    key = inkey().lower()
                except Exception:
                    print("\n")
                    log.error(
                        "Console does not support getting a character"
                        " as required by more().  Use pprint() instead."
                    )
                    return
                if key in allowed_keys:
                    break
            print(key)

            if key.lower() == "q":
                break

            if key == " " or key == "f":
                i0 += delta_lines
            elif key == "b":
                i0 = i0 - delta_lines
            elif key == "r":
                pass
            elif key == "<":
                i0 = 0
            elif key == ">":
                i0 = len(tabcol)
            elif key == "p":
                i0 -= 1
            elif key == "n":
                i0 += 1
            elif key == "h":
                showlines = False
                print(
                    """
    Browsing keys:
       f, <space> : forward one page
       b : back one page
       r : refresh same page
       n : next row
       p : previous row
       < : go to beginning
       > : go to end
       q : quit browsing
       h : print this help""",
                    end=" ",
                )
            if i0 < 0:
                i0 = 0
            if i0 >= len(tabcol) - delta_lines:
                i0 = len(tabcol) - delta_lines
            print("\n")
