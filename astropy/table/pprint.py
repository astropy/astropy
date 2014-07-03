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
from ..utils.console import Getch, color_print, terminal_size, conf

if six.PY3:
    def default_format_func(format_, val):
        if isinstance(val, bytes):
            return val.decode('utf-8')
        else:
            return str(val)
    _format_funcs = {None: default_format_func}
elif six.PY2:
    _format_funcs = {None: lambda format_, val: text_type(val)}


def _auto_format_func(format_, val):
    """Format ``val`` according to ``format_`` for both old- and new-
    style format specifications or using a user supplied function.
    More importantly, determine and cache (in _format_funcs) a function
    that will do this subsequently.  In this way this complicated logic is
    only done for the first value.

    Returns the formatted value.
    """
    if six.callable(format_):
        format_func = lambda format_, val: format_(val.tolist())
        try:
            out = format_func(format_, val)
            if not isinstance(out, six.string_types):
                raise ValueError('Format function for value {0} returned {1} instead of string type'
                                 .format(val, type(val)))
        except Exception as err:
            raise ValueError('Format function for value {0} failed: {1}'
                             .format(val, err))
    else:
        try:
            # Convert val to Python object with tolist().  See
            # https://github.com/astropy/astropy/issues/148#issuecomment-3930809
            out = format_.format(val.tolist())
            # Require that the format statement actually did something
            if out == format_:
                raise ValueError
            format_func = lambda format_, val: format_.format(val.tolist())
        except:  # Not sure what exceptions might be raised
            try:
                out = format_ % val
                if out == format_:
                    raise ValueError
                format_func = lambda format_, val: format_ % val
            except:
                raise ValueError('Unable to parse format string {0}'
                                 .format(format_))
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
        if max_lines < 6:
            max_lines = 6

        if max_width is None:
            max_width = width
        elif max_width < 0:
            max_width = sys.maxsize
        if max_width < 10:
            max_width = 10

        return max_lines, max_width


    def _pformat_col(self, col, max_lines=None, show_name=True, show_unit=None):
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

        Returns
        -------
        lines : list
            List of lines with formatted column values

        n_header : int
            Number of lines in the header

        """
        outs = {}  # Some values from _pformat_col_iter iterator that are needed here
        col_strs = list(self._pformat_col_iter(col, max_lines, show_name, show_unit, outs))
        col_width = max(len(x) for x in col_strs)

        # Center line content and generate dashed headerline
        for i in outs['i_centers']:
            col_strs[i] = col_strs[i].center(col_width)
        if outs['i_dashes'] is not None:
            col_strs[outs['i_dashes']] = '-' * col_width

        # Now bring all the column string values to the same fixed width
        for i, col_str in enumerate(col_strs):
            col_strs[i] = col_str.rjust(col_width)

        return col_strs, outs['n_header']


    def _pformat_col_iter(self, col, max_lines, show_name, show_unit, outs):
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

        out : dict
            Must be a dict which is used to pass back additional values
            defined within the iterator.
        """
        max_lines, _ = self._get_pprint_size(max_lines, -1)

        multidims = col.shape[1:]
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
            col_name = six.text_type(col.name)
            if multidims:
                col_name += ' [{0}]'.format(
                    ','.join(six.text_type(n) for n in multidims))
            n_header += 1
            yield col_name
        if show_unit:
            i_centers.append(n_header)
            n_header += 1
            yield six.text_type(col.unit or '')
        if show_unit or show_name:
            i_dashes = n_header
            n_header += 1
            yield '---'

        max_lines -= n_header
        n_print2 = max_lines // 2
        n_rows = len(col)

        format_func = _format_funcs.get(col.format, _auto_format_func)
        if len(col) > max_lines:
            i0 = n_print2
            i1 = n_rows - n_print2 - max_lines % 2
        else:
            i0 = len(col)
            i1 = 0

        # Add formatted values if within bounds allowed by max_lines
        for i in xrange(n_rows):
            if i < i0 or i > i1:
                if multidims:
                    # Prevents colums like Column(data=[[(1,)],[(2,)]], name='a')
                    # with shape (n,1,...,1) from being printed as if there was
                    # more than one element in a row
                    if trivial_multidims:
                        col_str = format_func(col.format, col[(i,) + multidim0])
                    else:
                        col_str = (format_func(col.format, col[(i,) + multidim0]) +
                                  ' .. ' +
                                  format_func(col.format, col[(i,) + multidim1]))
                else:
                    col_str = format_func(col.format, col[i])
                yield col_str
            elif i == i0:
                yield '...'

        outs['n_header'] = n_header
        outs['i_centers'] = i_centers
        outs['i_dashes'] = i_dashes


    def _pformat_table(self, table, max_lines=None, max_width=None, show_name=True,
                       show_unit=None, html=False, tableid=None):
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

        html : bool
            Format the output as an HTML table (default=False)

        tableid : str or None
            An ID tag for the table; only used if html is set.  Default is
            "table{id}", where id is the unique integer id of the table object,
            id(table)

        Returns
        -------
        out : str
            Formatted table as a single string

        n_header : int
            Number of lines in the header
        """
        # "Print" all the values into temporary lists by column for subsequent
        # use and to determine the width
        max_lines, max_width = self._get_pprint_size(max_lines, max_width)
        cols = []

        if show_unit is None:
            show_unit = any([col.unit for col in six.itervalues(table.columns)])

        for col in six.itervalues(table.columns):
            lines, n_header = self._pformat_col(col, max_lines, show_name,
                                                show_unit)
            cols.append(lines)

        if not cols:
            return []

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

        return rows, n_header


    def _more_tabcol(self, tabcol, max_lines=None, max_width=None, show_name=True,
                     show_unit=None):
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
        """
        allowed_keys = 'f br<>qhpn'

        # Count the header lines
        n_header = 0
        if show_name:
            n_header += 1
        if show_unit:
            n_header += 1
        if show_name or show_unit:
            n_header += 1

        # Set up kwargs for pformat call.  Only Table gets max_width.
        kwargs = dict(max_lines=-1, show_name=show_name, show_unit=show_unit)
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
