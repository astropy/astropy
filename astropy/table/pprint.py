# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from ..extern import six
from ..extern.six import text_type
from ..extern.six.moves import zip as izip
from ..extern.six.moves import xrange

import os
import sys

import numpy as np

from .. import log
# Note, in numpy <= 1.6, some classes do not properly represent themselves.
from ..utils.compat import NUMPY_LT_1_6_1
from ..utils.console import Getch, color_print, terminal_size, conf
from ..utils.data_info import dtype_info_name

if six.PY3:
    def default_format_func(format_, val):
        if isinstance(val, bytes):
            return val.decode('utf-8')
        else:
            return str(val)
    _format_funcs = {None: default_format_func}
elif six.PY2:
    _format_funcs = {None: lambda format_, val: text_type(val)}


### The first three functions are helpers for _auto_format_func


def _use_val_tolist(format_func):
    """Wrap format function to work with values converted to python equivalents.

    In numpy <= 1.6, classes such as np.float32 do not properly represent
    themselves as floats, and hence cannot easily be formatted; see
    https://github.com/astropy/astropy/issues/148#issuecomment-3930809
    Hence, we force the value to a python type using tolist()
    (except for np.ma.masked, since np.ma.masked.tolist() is None).
    """
    return lambda format_, val: format_func(format_,
                                            val if val is np.ma.masked
                                            else val.tolist())

def _use_str_for_masked_values(format_func):
    """Wrap format function to trap masked values.

    String format functions and most user functions will not be able to deal
    with masked values, so we wrap them to ensure they are passed to str().
    """
    return lambda format_, val: (str(val) if val is np.ma.masked
                                 else format_func(format_, val))

def _possible_string_format_functions(format_):
    """Iterate through possible string-derived format functions.

    A string can either be a format specifier for the format built-in,
    a new-style format string, or an old-style format string.
    """
    yield lambda format_, val: format(val, format_)
    yield lambda format_, val: format_.format(val)
    yield lambda format_, val: format_ % val

def _auto_format_func(format_, val):
    """Format ``val`` according to ``format_`` for a plain format specifier,
    old- or new-style format strings, or using a user supplied function.
    More importantly, determine and cache (in _format_funcs) a function
    that will do this subsequently.  In this way this complicated logic is
    only done for the first value.

    Returns the formatted value.
    """
    if format_ in _format_funcs:
        return _format_funcs[format_](format_, val)

    if six.callable(format_):
        format_func = lambda format_, val: format_(val)
        if NUMPY_LT_1_6_1:
            format_func = _use_val_tolist(format_func)
        try:
            out = format_func(format_, val)
            if not isinstance(out, six.string_types):
                raise ValueError('Format function for value {0} returned {1} '
                                 'instead of string type'
                                 .format(val, type(val)))
        except Exception as err:
            # For a masked element, the format function call likely failed
            # to handle it.  Just return the string representation for now,
            # and retry when a non-masked value comes along.
            if val is np.ma.masked:
                return str(val)

            raise ValueError('Format function for value {0} failed: {1}'
                             .format(val, err))
        # If the user-supplied function handles formatting masked elements, use
        # it directly.  Otherwise, wrap it in a function that traps them.
        try:
            format_func(format_, np.ma.masked)
        except:
            format_func = _use_str_for_masked_values(format_func)
    else:
        # For a masked element, we cannot set string-based format functions yet,
        # as all tests below will fail.  Just return the string representation
        # of masked for now, and retry when a non-masked value comes along.
        if val is np.ma.masked:
            return str(val)

        for format_func in _possible_string_format_functions(format_):
            if NUMPY_LT_1_6_1:
                format_func = _use_val_tolist(format_func)

            try:
                # Does this string format method work?
                out = format_func(format_, val)
                # Require that the format statement actually did something.
                assert out != format_
            except:
                continue
            else:
                break
        else:
            # None of the possible string functions passed muster.
            raise ValueError('Unable to parse format string {0}'
                             .format(format_))

        # String-based format functions will fail on masked elements;
        # wrap them in a function that traps them.
        format_func = _use_str_for_masked_values(format_func)

    _format_funcs[format_] = format_func
    return out


class TableFormatter(object):
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

    def _pformat_col(self, col, max_lines=None, show_name=True, show_unit=None,
                     show_dtype=False, show_length=None, html=False, align=None):
        """Return a list of formatted string representation of column values.

        Parameters
        ----------
        max_lines : int
            Maximum lines of output (header + data rows)

        show_name : bool
            Include column name (default=True)

        show_unit : bool
            Include a header row for unit.  Default is to show a row
            for units only if one or more columns has a defined value
            for the unit.

        show_dtype : bool
            Include column dtype (default=False)

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
        col_strs_iter = self._pformat_col_iter(col, max_lines, show_name=show_name,
                                               show_unit=show_unit,
                                               show_dtype=show_dtype,
                                               show_length=show_length,
                                               outs=outs)
        col_strs = list(col_strs_iter)
        if len(col_strs) > 0:
            col_width = max(len(x) for x in col_strs)

        if html:
            from ..utils.xml.writer import xml_escape
            n_header = outs['n_header']
            for i, col_str in enumerate(col_strs):
                # _pformat_col output has a header line '----' which is not needed here
                if i == n_header - 1:
                    continue
                td = 'th' if i < n_header else 'td'
                val = '<{0}>{1}</{2}>'.format(td, xml_escape(col_str.strip()), td)
                row = ('<tr>' + val + '</tr>')
                if i < n_header:
                    row = ('<thead>' + row + '</thead>')
                col_strs[i] = row

            if n_header > 0:
                # Get rid of '---' header line
                col_strs.pop(n_header - 1)
            col_strs.insert(0, '<table>')
            col_strs.append('</table>')

        # Now bring all the column string values to the same fixed width
        else:
            col_width = max(len(x) for x in col_strs) if col_strs else 1

            # Center line header content and generate dashed headerline
            for i in outs['i_centers']:
                col_strs[i] = col_strs[i].center(col_width)
            if outs['i_dashes'] is not None:
                col_strs[outs['i_dashes']] = '-' * col_width

            # Format columns according to alignment.  `align` arg has precedent, otherwise
            # use `col.format` if it is a legal alignment string.  If neither applies
            # then right justify.
            justify_methods = {'<': 'ljust', '^': 'center', '>': 'rjust', '0=': 'zfill'}
            align = align or (col.info.format
                              if col.info.format in justify_methods
                              else '>')

            # This can only occur if `align` was explicitly provided.
            if align not in justify_methods:
                raise ValueError("column align must be one of '<', '^', '>', or '0='")

            justify_method = justify_methods[align]
            for i, col_str in enumerate(col_strs):
                col_strs[i] = getattr(col_str, justify_method)(col_width)

        if outs['show_length']:
            col_strs.append('Length = {0} rows'.format(len(col)))

        return col_strs, outs

    def _pformat_col_iter(self, col, max_lines, show_name, show_unit, outs,
                          show_dtype=False, show_length=None):
        """Iterator which yields formatted string representation of column values.

        Parameters
        ----------
        max_lines : int
            Maximum lines of output (header + data rows)

        show_name : bool
            Include column name (default=True)

        show_unit : bool
            Include a header row for unit.  Default is to show a row
            for units only if one or more columns has a defined value
            for the unit.

        outs : dict
            Must be a dict which is used to pass back additional values
            defined within the iterator.

        show_dtype : bool
            Include column dtype (default=False)

        show_length : bool
            Include column length at end.  Default is to show this only
            if the column is not shown completely.
        """
        max_lines, _ = self._get_pprint_size(max_lines, -1)

        multidims = getattr(col, 'shape', [0])[1:]
        if multidims:
            multidim0 = tuple(0 for n in multidims)
            multidim1 = tuple(n - 1 for n in multidims)
            trivial_multidims = np.prod(multidims) == 1

        i_dashes = None
        i_centers = []  # Line indexes where content should be centered
        n_header = 0
        if show_name:
            i_centers.append(n_header)
            # Get column name (or 'None' if not set)
            col_name = six.text_type(col.info.name)
            if multidims:
                col_name += ' [{0}]'.format(
                    ','.join(six.text_type(n) for n in multidims))
            n_header += 1
            yield col_name
        if show_unit:
            i_centers.append(n_header)
            n_header += 1
            yield six.text_type(col.info.unit or '')
        if show_dtype:
            i_centers.append(n_header)
            n_header += 1
            try:
                dtype = dtype_info_name(col.dtype)
            except AttributeError:
                dtype = 'object'
            yield six.text_type(dtype)
        if show_unit or show_name or show_dtype:
            i_dashes = n_header
            n_header += 1
            yield '---'

        max_lines -= n_header
        n_print2 = max_lines // 2
        n_rows = len(col)

        col_format = col.info.format or getattr(col.info, 'default_format', None)
        format_func = _format_funcs.get(col_format, _auto_format_func)
        if len(col) > max_lines:
            if show_length is None:
                show_length = True
            i0 = n_print2 - (1 if show_length else 0)
            i1 = n_rows - n_print2 - max_lines % 2
            ii = np.concatenate([np.arange(0, i0 + 1), np.arange(i1 + 1, len(col))])
        else:
            i0 = -1
            ii = np.arange(len(col))

        # Add formatted values if within bounds allowed by max_lines
        for i in ii:
            if i == i0:
                yield '...'
            else:
                if multidims:
                    # Prevents columns like Column(data=[[(1,)],[(2,)]], name='a')
                    # with shape (n,1,...,1) from being printed as if there was
                    # more than one element in a row
                    if trivial_multidims:
                        col_str = format_func(col_format, col[(i,) + multidim0])
                    else:
                        col_str = (format_func(col_format, col[(i,) + multidim0]) +
                                  ' .. ' +
                                  format_func(col_format, col[(i,) + multidim1]))
                else:
                    col_str = format_func(col_format, col[i])
                yield col_str

        outs['show_length'] = show_length
        outs['n_header'] = n_header
        outs['i_centers'] = i_centers
        outs['i_dashes'] = i_dashes

    def _pformat_table(self, table, max_lines=None, max_width=None,
                       show_name=True, show_unit=None, show_dtype=False,
                       html=False, tableid=None, tableclass=None, align=None):
        """Return a list of lines for the formatted string representation of
        the table.

        Parameters
        ----------
        max_lines : int or None
            Maximum number of rows to output

        max_width : int or None
            Maximum character width of output

        show_name : bool
            Include a header row for column names (default=True)

        show_unit : bool
            Include a header row for unit.  Default is to show a row
            for units only if one or more columns has a defined value
            for the unit.

        show_dtype : bool
            Include a header row for column dtypes (default=False)

        html : bool
            Format the output as an HTML table (default=False)

        tableid : str or None
            An ID tag for the table; only used if html is set.  Default is
            "table{id}", where id is the unique integer id of the table object,
            id(table)

        tableclass : str or list of str or `None`
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
        cols = []

        if show_unit is None:
            show_unit = any([col.info.unit for col in six.itervalues(table.columns)])

        # Coerce align into a correctly-sized list of alignments (if possible)
        n_cols = len(table.columns)
        if align is None:
            align = [None] * n_cols

        if isinstance(align, six.string_types):
            align = [align]

        if isinstance(align, (list, tuple)):
            if len(align) == 1:
                align = align * n_cols
            elif len(align) != n_cols:
                raise ValueError('got {0} alignment values instead of 1 or '
                                 'the number of columns ({1})'
                                 .format(len(align), n_cols))
        else:
            raise TypeError('align keyword must be str or list or tuple (got {0})'
                            .format(type(align)))

        for align_, col in izip(align, table.columns.values()):
            lines, outs = self._pformat_col(col, max_lines, show_name=show_name,
                                            show_unit=show_unit, show_dtype=show_dtype,
                                            align=align_)
            if outs['show_length']:
                lines = lines[:-1]
            cols.append(lines)

        if not cols:
            return ['<No columns>'], {'show_length': False}

        # Use the values for the last column since they are all the same
        n_header = outs['n_header']

        n_rows = len(cols[0])
        outwidth = lambda cols: sum(len(c[0]) for c in cols) + len(cols) - 1
        dots_col = ['...'] * n_rows
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
            from ..utils.xml.writer import xml_escape

            if tableid is None:
                tableid = 'table{id}'.format(id=id(table))

            if tableclass is not None:
                if isinstance(tableclass, list):
                    tableclass = ' '.join(tableclass)
                rows.append('<table id="{tid}" class="{tcls}">'.format(
                    tid=tableid, tcls=tableclass))
            else:
                rows.append('<table id="{tid}">'.format(tid=tableid))

            for i in range(n_rows):
                # _pformat_col output has a header line '----' which is not needed here
                if i == n_header - 1:
                    continue
                td = 'th' if i < n_header else 'td'
                vals = ('<{0}>{1}</{2}>'.format(td, xml_escape(col[i].strip()), td)
                        for col in cols)
                row = ('<tr>' + ''.join(vals) + '</tr>')
                if i < n_header:
                    row = ('<thead>' + row + '</thead>')
                rows.append(row)
            rows.append('</table>')
        else:
            for i in range(n_rows):
                row = ' '.join(col[i] for col in cols)
                rows.append(row)

        return rows, outs

    def _more_tabcol(self, tabcol, max_lines=None, max_width=None,
                     show_name=True, show_unit=None, show_dtype=False):
        """Interactive "more" of a table or column.

        Parameters
        ----------
        max_lines : int or None
            Maximum number of rows to output

        max_width : int or None
            Maximum character width of output

        show_name : bool
            Include a header row for column names (default=True)

        show_unit : bool
            Include a header row for unit.  Default is to show a row
            for units only if one or more columns has a defined value
            for the unit.

        show_dtype : bool
            Include a header row for column dtypes (default=False)
        """
        allowed_keys = 'f br<>qhpn'

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
        kwargs = dict(max_lines=-1, show_name=show_name, show_unit=show_unit,
                      show_dtype=show_dtype)
        if hasattr(tabcol, 'columns'):  # tabcol is a table
            kwargs['max_width'] = max_width

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
                    os.system('cls' if os.name == 'nt' else 'clear')
                except:
                    pass  # No worries if clear screen call fails
                lines = tabcol[i0:i1].pformat(**kwargs)
                colors = ('red' if i < n_header else 'default'
                          for i in xrange(len(lines)))
                for color, line in izip(colors, lines):
                    color_print(line, color)
            showlines = True
            print()
            print("-- f, <space>, b, r, p, n, <, >, q h (help) --", end=' ')
            # Get a valid key
            while True:
                try:
                    key = inkey().lower()
                except:
                    print("\n")
                    log.error('Console does not support getting a character'
                              ' as required by more().  Use pprint() instead.')
                    return
                if key in allowed_keys:
                    break
            print(key)

            if key.lower() == 'q':
                break
            elif key == ' ' or key == 'f':
                i0 += delta_lines
            elif key == 'b':
                i0 = i0 - delta_lines
            elif key == 'r':
                pass
            elif key == '<':
                i0 = 0
            elif key == '>':
                i0 = len(tabcol)
            elif key == 'p':
                i0 -= 1
            elif key == 'n':
                i0 += 1
            elif key == 'h':
                showlines = False
                print("""
    Browsing keys:
       f, <space> : forward one page
       b : back one page
       r : refresh same page
       n : next row
       p : previous row
       < : go to beginning
       > : go to end
       q : quit browsing
       h : print this help""", end=' ')
            if i0 < 0:
                i0 = 0
            if i0 >= len(tabcol) - delta_lines:
                i0 = len(tabcol) - delta_lines
            print("\n")
