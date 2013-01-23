# Licensed under a 3-clause BSD style license - see LICENSE.rst
import os
import sys
from itertools import izip

from .. import log
from ..utils.console import Getch, color_print
from ..config import ConfigurationItem

_format_funcs = {None: lambda format_, val: str(val)}

MAX_LINES = ConfigurationItem('max_lines', 25, 'Maximum number of lines for '
    'the pretty-printer to use if it cannot determine the terminal size. '
    'Negative numbers mean no limit.')
MAX_WIDTH = ConfigurationItem('max_width', 80, 'Maximum number of characters '
    'for the pretty-printer to use per line if it cannot determine the '
    'terminal size.  Negative numbers mean no limit.')

def _get_pprint_size(max_lines=None, max_width=None):
    """Get the output size (number of lines and character width) for Column and
    Table pformat/pprint methods.

    If no value of `max_lines` is supplied then the height of the screen
    terminal is used to set `max_lines`.  If the terminal height cannot be
    determined then the default will be determined using the
    `astropy.table.pprint.MAX_LINES` configuration item. If a negative value
    of `max_lines` is supplied then there is no line limit applied.

    The same applies for max_width except the configuration item is
    `astropy.table.pprint.MAX_WIDTH`.

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
    if max_lines is None or max_width is None:
        try:  # Will likely fail on Windows
            import termios
            import fcntl
            import struct
            s = struct.pack("HHHH", 0, 0, 0, 0)
            fd_stdout = sys.stdout.fileno()
            x = fcntl.ioctl(fd_stdout, termios.TIOCGWINSZ, s)
            (lines, width, xpixels, ypixels) = struct.unpack("HHHH", x)
            if lines > 12:
                lines -= 6
            if width > 10:
                width -= 1
        except:
            lines, width = MAX_LINES(), MAX_WIDTH()

    if max_lines is None:
        max_lines = lines
    elif max_lines < 0:
        max_lines = sys.maxint
    if max_lines < 6:
        max_lines = 6

    if max_width is None:
        max_width = width
    elif max_width < 0:
        max_width = sys.maxint
    if max_width < 10:
        max_width = 10

    return max_lines, max_width


def _auto_format_func(format_, val):
    """Format ``val`` according to ``format_`` for both old- and new-
    style format specifications.  More importantly, determine and cache
    (in _format_funcs) a function that will do this subsequently.  In
    this way this complicated logic is only done for the first value.

    Returns the formatted value.
    """
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


def _pformat_col(col, max_lines=None, show_name=True, show_units=False):
    """Return a list of formatted string representation of column values.

    Parameters
    ----------
    max_lines : int
        Maximum lines of output (header + data rows)

    show_name : bool
        Include column name (default=True)

    show_units : bool
        Include a header row for units (default=False)

    Returns
    -------
    lines : list
        List of lines with formatted column values

    n_header : int
        Number of lines in the header

    """
    outs = {}  # Some values from _pformat_col_iter iterator that are needed here
    col_strs = list(_pformat_col_iter(col, max_lines, show_name, show_units, outs))
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


def _pformat_col_iter(col, max_lines, show_name, show_units, outs):
    """Iterator which yields formatted string representation of column values.

    Parameters
    ----------
    max_lines : int
        Maximum lines of output (header + data rows)

    show_name : bool
        Include column name (default=True)

    show_units : bool
        Include a header row for units (default=False)

    out : dict
        Must be a dict which is used to pass back additional values
        defined within the iterator.
    """
    max_lines, _ = _get_pprint_size(max_lines, -1)

    multidims = col.shape[1:]
    if multidims:
        multidim0 = tuple(0 for n in multidims)
        multidim1 = tuple(n - 1 for n in multidims)

    col_strs = []  # List of formatted column values
    i_dashes = None
    i_centers = []  # Line indexes where content should be centered
    n_header = 0
    if show_name:
        i_centers.append(n_header)
        if multidims:
            col_name = col.name + ' [{0}]'.format(
                ','.join(str(n) for n in multidims))
        else:
            col_name = col.name
        n_header += 1
        yield col_name
    if show_units:
        i_centers.append(n_header)
        n_header += 1
        yield str(col.units or '')
    if show_units or show_name:
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

def _pformat_table(table, max_lines=None, max_width=None, show_name=True,
                   show_units=False, html=False):
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

    show_units : bool
        Include a header row for units (default=False)

    html : bool
        Format the output as an HTML table (default=False)

    Returns
    -------
    out : str
        Formatted table as a single string

    n_header : int
        Number of lines in the header
    """
    # "Print" all the values into temporary lists by column for subsequent
    # use and to determine the width
    max_lines, max_width = _get_pprint_size(max_lines, max_width)
    cols = []
    for col in table.columns.values():
        lines, n_header = _pformat_col(col, max_lines, show_name,
                                       show_units)
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
        rows.append('<table>')
        for i in range(n_rows):
            # _pformat_col output has a header line '----' which is not needed here
            if i == n_header - 1:
                continue
            td = 'th' if i < n_header else 'td'
            vals = ('<{0}>{1}</{2}>'.format(td, col[i].strip(), td) for col in cols)
            row = ('<tr>' + ''.join(vals) + '</tr>')
            rows.append(row)
        rows.append('</table>')
    else:
        for i in range(n_rows):
            row = ' '.join(col[i] for col in cols)
            rows.append(row)

    return rows, n_header


def _more_tabcol(tabcol, max_lines=None, max_width=None, show_name=True,
                show_units=False):
    """Interactive "more" of a table or column.

    Parameters
    ----------
    max_lines : int or None
        Maximum number of rows to output

    max_width : int or None
        Maximum character width of output

    show_name : bool
        Include a header row for column names (default=True)

    show_units : bool
        Include a header row for units (default=False)
    """
    allowed_keys = 'f br<>qhpn'

    # Count the header lines
    n_header = 0
    if show_name:
        n_header += 1
    if show_units:
        n_header += 1
    if show_name or show_units:
        n_header += 1

    # Set up kwargs for pformat call.  Only Table gets max_width.
    kwargs = dict(max_lines=-1, show_name=show_name, show_units=show_units)
    if hasattr(tabcol, 'columns'):  # tabcol is a table
        kwargs['max_width'] = max_width

    # If max_lines is None (=> query screen size) then increase by 2.
    # This is because get_pprint_size leaves 6 extra lines so that in
    # ipython you normally see the last input line.
    max_lines1, max_width = _get_pprint_size(max_lines, max_width)
    if max_lines == None:
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
        print
        print "-- f, <space>, b, r, p, n, <, >, q h (help) --",
        # Get a valid key
        while True:
            try:
                key = inkey().lower()
            except:
                print "\n"
                log.error('Console does not support getting a character'
                          ' as required by more().  Use pprint() instead.')
                return
            if key in allowed_keys:
                break
        print key

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
            print """
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
        if i0 < 0:
            i0 = 0
        if i0 >= len(tabcol) - delta_lines:
            i0 = len(tabcol) - delta_lines
        print "\n"
