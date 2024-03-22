# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""An extensible ASCII table reader and writer.

fixedwidth.py:
  Read or write a table with fixed width columns.

:Copyright: Smithsonian Astrophysical Observatory (2011)
:Author: Tom Aldcroft (aldcroft@head.cfa.harvard.edu)
"""

from . import basic, core
from .core import DefaultSplitter, InconsistentTableError


class FixedWidthSplitter(core.BaseSplitter):
    """
    Split line based on fixed start and end positions for each ``col`` in
    ``self.cols``.

    This class requires that the Header class will have defined ``col.start``
    and ``col.end`` for each column.  The reference to the ``header.cols`` gets
    put in the splitter object by the base Reader.read() function just in time
    for splitting data lines by a ``data`` object.

    Note that the ``start`` and ``end`` positions are defined in the pythonic
    style so line[start:end] is the desired substring for a column.  This splitter
    class does not have a hook for ``process_lines`` since that is generally not
    useful for fixed-width input.

    """

    delimiter_pad = ""
    bookend = False
    delimiter = "|"

    def __call__(self, lines):
        for line in lines:
            vals = [line[x.start : x.end] for x in self.cols]
            if self.process_val:
                yield [self.process_val(x) for x in vals]
            else:
                yield vals

    def join(self, vals, widths):
        pad = self.delimiter_pad or ""
        delimiter = self.delimiter or ""
        padded_delim = pad + delimiter + pad
        if self.bookend:
            bookend_left = delimiter + pad
            bookend_right = pad + delimiter
        else:
            bookend_left = ""
            bookend_right = ""
        vals = [" " * (width - len(val)) + val for val, width in zip(vals, widths)]
        return bookend_left + padded_delim.join(vals) + bookend_right


class FixedWidthHeaderSplitter(DefaultSplitter):
    """Splitter class that splits on ``|``."""

    delimiter = "|"


class FixedWidthHeader(basic.BasicHeader):
    """
    Fixed width table header reader.
    """

    splitter_class = FixedWidthHeaderSplitter
    """ Splitter class for splitting data lines into columns """
    position_line = None  # secondary header line position
    """ row index of line that specifies position (default = 1) """
    set_of_position_line_characters = set(r"""`~!#$%^&*-_+=\|":'""")

    def get_line(self, lines, index):
        for i, line in enumerate(self.process_lines(lines)):
            if i == index:
                break
        else:  # No header line matching
            raise InconsistentTableError("No header line found in table")
        return line

    def get_cols(self, lines):
        """
        Initialize the header Column objects from the table ``lines``.

        Based on the previously set Header attributes find or create the column names.
        Sets ``self.cols`` with the list of Columns.

        Parameters
        ----------
        lines : list
            List of table lines

        """
        header_rows = getattr(self, "header_rows", ["name"])

        # See "else" clause below for explanation of start_line and position_line
        start_line = core._get_line_index(self.start_line, self.process_lines(lines))
        position_line = core._get_line_index(
            self.position_line, self.process_lines(lines)
        )

        # If start_line is none then there is no header line.  Column positions are
        # determined from first data line and column names are either supplied by user
        # or auto-generated.
        if start_line is None:
            if position_line is not None:
                raise ValueError(
                    "Cannot set position_line without also setting header_start"
                )

            # data.data_lines attribute already set via self.data.get_data_lines(lines)
            # in BaseReader.read().  This includes slicing for data_start / data_end.
            data_lines = self.data.data_lines

            if not data_lines:
                raise InconsistentTableError(
                    "No data lines found so cannot autogenerate column names"
                )
            vals, starts, ends = self.get_fixedwidth_params(data_lines[0])

            self.names = [self.auto_format.format(i) for i in range(1, len(vals) + 1)]

        else:
            # This bit of code handles two cases:
            # start_line = <index> and position_line = None
            #    Single header line where that line is used to determine both the
            #    column positions and names.
            # start_line = <index> and position_line = <index2>
            #    Two header lines where the first line defines the column names and
            #    the second line defines the column positions

            if position_line is not None:
                # Define self.col_starts and self.col_ends so that the call to
                # get_fixedwidth_params below will use those to find the header
                # column names.  Note that get_fixedwidth_params returns Python
                # slice col_ends but expects inclusive col_ends on input (for
                # more intuitive user interface).
                line = self.get_line(lines, position_line)
                if len(set(line) - {self.splitter.delimiter, " "}) != 1:
                    raise InconsistentTableError(
                        "Position line should only contain delimiters and "
                        'one other character, e.g. "--- ------- ---".'
                    )
                    # The line above lies. It accepts white space as well.
                    # We don't want to encourage using three different
                    # characters, because that can cause ambiguities, but white
                    # spaces are so common everywhere that practicality beats
                    # purity here.
                charset = self.set_of_position_line_characters.union(
                    {self.splitter.delimiter, " "}
                )
                if not set(line).issubset(charset):
                    raise InconsistentTableError(
                        f"Characters in position line must be part of {charset}"
                    )
                vals, self.col_starts, col_ends = self.get_fixedwidth_params(line)
                self.col_ends = [x - 1 if x is not None else None for x in col_ends]

            # Get the column names from the header line
            line = self.get_line(lines, start_line + header_rows.index("name"))
            self.names, starts, ends = self.get_fixedwidth_params(line)

        self._set_cols_from_names()

        for ii, attr in enumerate(header_rows):
            if attr != "name":
                line = self.get_line(lines, start_line + ii)
                vals = self.get_fixedwidth_params(line)[0]
                for col, val in zip(self.cols, vals):
                    if val:
                        setattr(col, attr, val)

        # Set column start and end positions.
        for i, col in enumerate(self.cols):
            col.start = starts[i]
            col.end = ends[i]

    def get_fixedwidth_params(self, line):
        """
        Split ``line`` on the delimiter and determine column values and
        column start and end positions.  This might include null columns with
        zero length (e.g. for ``header row = "| col1 || col2 | col3 |"`` or
        ``header2_row = "----- ------- -----"``).  The null columns are
        stripped out.  Returns the values between delimiters and the
        corresponding start and end positions.

        Parameters
        ----------
        line : str
            Input line

        Returns
        -------
        vals : list
            List of values.
        starts : list
            List of starting indices.
        ends : list
            List of ending indices.

        """
        # If column positions are already specified then just use those.
        # If neither column starts or ends are given, figure out positions
        # between delimiters. Otherwise, either the starts or the ends have
        # been given, so figure out whichever wasn't given.
        if self.col_starts is not None and self.col_ends is not None:
            starts = list(self.col_starts)  # could be any iterable, e.g. np.array
            # user supplies inclusive endpoint
            ends = [x + 1 if x is not None else None for x in self.col_ends]
            if len(starts) != len(ends):
                raise ValueError(
                    "Fixed width col_starts and col_ends must have the same length"
                )
            vals = [line[start:end].strip() for start, end in zip(starts, ends)]
        elif self.col_starts is None and self.col_ends is None:
            # There might be a cleaner way to do this but it works...
            vals = line.split(self.splitter.delimiter)
            starts = [0]
            ends = []
            for val in vals:
                if val:
                    ends.append(starts[-1] + len(val))
                    starts.append(ends[-1] + 1)
                else:
                    starts[-1] += 1
            starts = starts[:-1]
            vals = [x.strip() for x in vals if x]
            if len(vals) != len(starts) or len(vals) != len(ends):
                raise InconsistentTableError("Error parsing fixed width header")
        else:
            # exactly one of col_starts or col_ends is given...
            if self.col_starts is not None:
                starts = list(self.col_starts)
                ends = starts[1:] + [None]  # Assume each col ends where the next starts
            else:  # self.col_ends is not None
                ends = [x + 1 for x in self.col_ends]
                starts = [0] + ends[:-1]  # Assume each col starts where the last ended
            vals = [line[start:end].strip() for start, end in zip(starts, ends)]

        return vals, starts, ends

    def write(self, lines):
        # Header line not written until data are formatted.  Until then it is
        # not known how wide each column will be for fixed width.
        pass


class FixedWidthData(basic.BasicData):
    """
    Base table data reader.
    """

    splitter_class = FixedWidthSplitter
    """ Splitter class for splitting data lines into columns """
    start_line = None

    def write(self, lines):
        default_header_rows = [] if self.header.start_line is None else ["name"]
        header_rows = getattr(self, "header_rows", default_header_rows)
        # First part is getting the widths of each column.
        # List (rows) of list (column values) for data lines
        vals_list = list(zip(*self.str_vals()))

        # List (rows) of list (columns values) for header lines.
        hdrs_list = []
        for col_attr in header_rows:
            vals = [
                "" if (val := getattr(col.info, col_attr)) is None else str(val)
                for col in self.cols
            ]
            hdrs_list.append(vals)

        # Widths for data columns
        widths = [
            max(len(vals[i_col]) for vals in vals_list)
            for i_col in range(len(self.cols))
        ]
        # Incorporate widths for header columns (if there are any)
        if hdrs_list:
            for i_col in range(len(self.cols)):
                widths[i_col] = max(
                    widths[i_col], *(len(vals[i_col]) for vals in hdrs_list)
                )

        # Now collect formatted header and data lines into the output lines
        for vals in hdrs_list:
            lines.append(self.splitter.join(vals, widths))

        if self.header.position_line is not None:
            vals = [self.header.position_char * width for width in widths]
            lines.append(self.splitter.join(vals, widths))

        for vals in vals_list:
            lines.append(self.splitter.join(vals, widths))

        return lines


class FixedWidth(basic.Basic):
    """Fixed width table with single header line defining column names and positions.

    Examples::

      # Bar delimiter in header and data

      |  Col1 |   Col2      |  Col3 |
      |  1.2  | hello there |     3 |
      |  2.4  | many words  |     7 |

      # Bar delimiter in header only

      Col1 |   Col2      | Col3
      1.2    hello there    3
      2.4    many words     7

      # No delimiter with column positions specified as input

      Col1       Col2Col3
       1.2hello there   3
       2.4many words    7

    See the :ref:`astropy:fixed_width_gallery` for specific usage examples.

    """

    _format_name = "fixed_width"
    _description = "Fixed width"

    header_class = FixedWidthHeader
    data_class = FixedWidthData

    def __init__(
        self,
        col_starts=None,
        col_ends=None,
        delimiter_pad=" ",
        bookend=True,
        header_rows=None,
    ):
        if header_rows is None:
            header_rows = ["name"]
        super().__init__()
        self.data.splitter.delimiter_pad = delimiter_pad
        self.data.splitter.bookend = bookend
        self.header.col_starts = col_starts
        self.header.col_ends = col_ends
        self.header.header_rows = header_rows
        self.data.header_rows = header_rows
        if self.data.start_line is None:
            self.data.start_line = len(header_rows)


class FixedWidthNoHeaderHeader(FixedWidthHeader):
    """Header reader for fixed with tables with no header line."""

    start_line = None


class FixedWidthNoHeaderData(FixedWidthData):
    """Data reader for fixed width tables with no header line."""

    start_line = 0


class FixedWidthNoHeader(FixedWidth):
    """Fixed width table which has no header line.

    When reading, column names are either input (``names`` keyword) or
    auto-generated.  Column positions are determined either by input
    (``col_starts`` and ``col_stops`` keywords) or by splitting the first data
    line.  In the latter case a ``delimiter`` is required to split the data
    line.

    Examples::

      # Bar delimiter in header and data

      |  1.2  | hello there |     3 |
      |  2.4  | many words  |     7 |

      # Compact table having no delimiter and column positions specified as input

      1.2hello there3
      2.4many words 7

    This class is just a convenience wrapper around the ``FixedWidth`` reader
    but with ``header_start=None`` and ``data_start=0``.

    See the :ref:`astropy:fixed_width_gallery` for specific usage examples.

    """

    _format_name = "fixed_width_no_header"
    _description = "Fixed width with no header"
    header_class = FixedWidthNoHeaderHeader
    data_class = FixedWidthNoHeaderData

    def __init__(self, col_starts=None, col_ends=None, delimiter_pad=" ", bookend=True):
        super().__init__(
            col_starts,
            col_ends,
            delimiter_pad=delimiter_pad,
            bookend=bookend,
            header_rows=[],
        )


class FixedWidthTwoLineHeader(FixedWidthHeader):
    """Header reader for fixed width tables splitting on whitespace.

    For fixed width tables with several header lines, there is typically
    a white-space delimited format line, so splitting on white space is
    needed.
    """

    splitter_class = DefaultSplitter


class FixedWidthTwoLineDataSplitter(FixedWidthSplitter):
    """Splitter for fixed width tables splitting on ``' '``."""

    delimiter = " "


class FixedWidthTwoLineData(FixedWidthData):
    """Data reader for fixed with tables with two header lines."""

    splitter_class = FixedWidthTwoLineDataSplitter


class FixedWidthTwoLine(FixedWidth):
    """Fixed width table which has two header lines.

    The first header line defines the column names and the second implicitly
    defines the column positions.

    Examples::

      # Typical case with column extent defined by ---- under column names.

       col1    col2         <== header_start = 0
      -----  ------------   <== position_line = 1, position_char = "-"
        1     bee flies     <== data_start = 2
        2     fish swims

      # Pretty-printed table

      +------+------------+
      | Col1 |   Col2     |
      +------+------------+
      |  1.2 | "hello"    |
      |  2.4 | there world|
      +------+------------+

    See the :ref:`astropy:fixed_width_gallery` for specific usage examples.

    """

    _format_name = "fixed_width_two_line"
    _description = "Fixed width with second header line"
    data_class = FixedWidthTwoLineData
    header_class = FixedWidthTwoLineHeader

    def __init__(
        self,
        position_line=None,
        position_char="-",
        delimiter_pad=None,
        bookend=False,
        header_rows=None,
    ):
        if len(position_char) != 1:
            raise ValueError(
                f'Position_char="{position_char}" must be a single character'
            )
        super().__init__(
            delimiter_pad=delimiter_pad, bookend=bookend, header_rows=header_rows
        )
        if position_line is None:
            position_line = len(self.header.header_rows)
        self.header.position_line = position_line
        self.header.position_char = position_char
        self.data.start_line = position_line + 1
