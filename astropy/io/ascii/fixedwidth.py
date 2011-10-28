"""Asciitable: an extensible ASCII table reader and writer.

fixedwidth.py:
  Read or write a table with fixed width columns.

:Copyright: Smithsonian Astrophysical Observatory (2011)
:Author: Tom Aldcroft (aldcroft@head.cfa.harvard.edu)
"""

## 
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are met:
##     * Redistributions of source code must retain the above copyright
##       notice, this list of conditions and the following disclaimer.
##     * Redistributions in binary form must reproduce the above copyright
##       notice, this list of conditions and the following disclaimer in the
##       documentation and/or other materials provided with the distribution.
##     * Neither the name of the Smithsonian Astrophysical Observatory nor the
##       names of its contributors may be used to endorse or promote products
##       derived from this software without specific prior written permission.
## 
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
## ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
## WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
## DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
## DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
## (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
## LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
## ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS  
## SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import re
import itertools
import asciitable.core as core
from asciitable.core import io, next, izip, any

class FixedWidthSplitter(core.BaseSplitter):
    """Split line based on fixed start and end positions for each ``col`` in
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
    delimiter_pad = ''
    bookend = False

    def __call__(self, lines):
        for line in lines:
            vals = [line[x.start:x.end] for x in self.cols]
            if self.process_val:
                yield [self.process_val(x) for x in vals]
            else:
                yield vals

    def join(self, vals, widths):
        pad = self.delimiter_pad or ''
        delimiter = self.delimiter or ''
        padded_delim = pad + delimiter + pad
        if self.bookend:
            bookend_left = delimiter + pad
            bookend_right = pad + delimiter
        else:
            bookend_left = ''
            bookend_right = ''
        vals = [' ' * (width - len(val)) + val for val, width in zip(vals, widths)]
        return bookend_left + padded_delim.join(vals) + bookend_right


class FixedWidthHeader(core.BaseHeader):
    """Fixed width table header reader.

    The key settable class attributes are:

    :param auto_format: format string for auto-generating column names
    :param start_line: None, int, or a function of ``lines`` that returns None or int
    :param comment: regular expression for comment lines
    :param splitter_class: Splitter class for splitting data lines into columns
    :param names: list of names corresponding to each data column
    :param include_names: list of names to include in output (default=None selects all names)
    :param exclude_names: list of names to exlude from output (applied after ``include_names``)
    :param position_line: row index of line that specifies position (default = 1)
    :param position_char: character used to write the position line (default = "-")
    :param col_starts: list of start positions for each column (0-based counting)
    :param col_ends: list of end positions (inclusive) for each column
    :param delimiter_pad: padding around delimiter when writing (default = None)
    :param bookend: put the delimiter at start and end of line when writing (default = False)
    """

    position_line = None   # secondary header line position

    def get_line(self, lines, index):
        for i, line in enumerate(self.process_lines(lines)):
            if i == index:
                break
        else: # No header line matching
            raise InconsistentTableError('No header line found in table')
        return line

    def get_cols(self, lines):
        """Initialize the header Column objects from the table ``lines``.

        Based on the previously set Header attributes find or create the column names.
        Sets ``self.cols`` with the list of Columns.  This list only includes the actual
        requested columns after filtering by the include_names and exclude_names
        attributes.  See ``self.names`` for the full list.

        :param lines: list of table lines
        :returns: None
        """

        # See "else" clause below for explanation of start_line and position_line
        start_line = core._get_line_index(self.start_line, self.process_lines(lines))
        position_line = core._get_line_index(self.position_line, self.process_lines(lines))

        # If start_line is none then there is no header line.  Column positions are
        # determined from first data line and column names are either supplied by user
        # or auto-generated.
        if start_line is None:
            if position_line is not None:
                raise ValueError("Cannot set position_line without also setting header_start")
            data_lines = self.data.process_lines(lines)
            if not data_lines:
                raise InconsistentTableError('No data lines found so cannot autogenerate column names')
            vals, starts, ends = self.get_fixedwidth_params(data_lines[0])

            if self.names is None:
                self.names = [self.auto_format % i for i in range(1, len(vals) + 1)]

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
                vals, self.col_starts, col_ends = self.get_fixedwidth_params(line)
                self.col_ends = [x - 1 for x in col_ends]

            # Get the header column names and column positions
            line = self.get_line(lines, start_line)
            vals, starts, ends = self.get_fixedwidth_params(line)

            # Possibly override the column names with user-supplied values
            if self.names is None:
                self.names = vals
        
        # Filter self.names using include_names and exclude_names, then create
        # the actual Column objects.
        self._set_cols_from_names()
        self.n_data_cols = len(self.cols)
        
        # Set column start and end positions.  Also re-index the cols because
        # the FixedWidthSplitter does NOT return the ignored cols (as is the
        # case for typical delimiter-based splitters)
        for i, col in enumerate(self.cols):
            col.start = starts[col.index]
            col.end = ends[col.index]
            col.index = i

    def get_fixedwidth_params(self, line):
        """Split ``line`` on the delimiter and determine column values and column
        start and end positions.  This might include null columns with zero length
        (e.g. for header row = "| col1 || col2 | col3 |" or
        header2_row = "-----    -------   -----").  The null columns are
        stripped out.  Returns the values between delimiters and the corresponding
        start and end positions.

        :param line: input line
        :returns: (vals, starts, ends)
        """

        # If column positions are already specified then just use those, otherwise
        # figure out positions between delimiters.
        if self.col_starts is not None and self.col_ends is not None:
            starts = list(self.col_starts)  # could be any iterable, e.g. np.array
            ends = [x + 1 for x in self.col_ends] # user supplies inclusive endpoint
            if len(starts) != len(ends):
                raise ValueError('Fixed width col_starts and col_ends must have the same length')
            vals = [line[start:end].strip() for start, end in zip(starts, ends)]
        else:
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
                raise InconsistentTableError('Error parsing fixed width header')

        return vals, starts, ends

    def write(self, lines):
        # Header line not written until data are formatted.  Until then it is
        # not known how wide each column will be for fixed width.
        pass


class FixedWidthData(core.BaseData):
    """Base table data reader.

    :param start_line: None, int, or a function of ``lines`` that returns None or int
    :param end_line: None, int, or a function of ``lines`` that returns None or int
    :param comment: Regular expression for comment lines
    :param splitter_class: Splitter class for splitting data lines into columns
    """

    splitter_class = FixedWidthSplitter

    def write(self, lines):
        formatters = []
        for col in self.cols:
            formatter = self.formats.get(col.name, self.default_formatter)
            if not hasattr(formatter, '__call__'):
                formatter = core._format_func(formatter)
            col.formatter = formatter
            
        vals_list = []
        # Col iterator does the formatting defined above so each val is a string
        # and vals is a tuple of strings for all columns of each row
        for vals in izip(*self.cols):
            vals_list.append(vals)
            
        for i, col in enumerate(self.cols):
            col.width = max([len(vals[i]) for vals in vals_list])
            if self.header.start_line is not None:
                col.width = max(col.width, len(col.name))

        widths = [col.width for col in self.cols]

        if self.header.start_line is not None:
            lines.append(self.splitter.join([col.name for col in self.cols], widths))

        if self.header.position_line is not None:
            char = self.header.position_char
            if len(char) != 1:
                raise ValueError('Position_char="%s" must be a single character' % char)
            vals = [char * col.width for col in self.cols]
            lines.append(self.splitter.join(vals, widths))

        for vals in vals_list:
            lines.append(self.splitter.join(vals, widths))

        return lines


class FixedWidth(core.BaseReader):
    """Read or write a fixed width table with a single header line that defines column
    names and positions.  Examples::

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

    See the :ref:`fixed_width_gallery` for specific usage examples.

    :param col_starts: list of start positions for each column (0-based counting)
    :param col_ends: list of end positions (inclusive) for each column
    :param delimiter_pad: padding around delimiter when writing (default = None)
    :param bookend: put the delimiter at start and end of line when writing (default = False)
    """
    def __init__(self, col_starts=None, col_ends=None, delimiter_pad=' ', bookend=True):
        core.BaseReader.__init__(self)

        self.header = FixedWidthHeader()
        self.data = FixedWidthData()
        self.data.header = self.header
        self.header.data = self.data

        self.header.splitter.delimiter = '|'
        self.data.splitter.delimiter = '|'
        self.data.splitter.delimiter_pad = delimiter_pad
        self.data.splitter.bookend = bookend
        self.header.start_line = 0
        self.data.start_line = 1
        self.header.comment = r'\s*#'
        self.header.write_comment = '# '
        self.data.comment = r'\s*#'
        self.data.write_comment = '# '
        self.header.col_starts = col_starts
        self.header.col_ends = col_ends


class FixedWidthNoHeader(FixedWidth):
    """Read or write a fixed width table which has no header line.  Column
    names are either input (``names`` keyword) or auto-generated.  Column
    positions are determined either by input (``col_starts`` and ``col_stops``
    keywords) or by splitting the first data line.  In the latter case a
    ``delimiter`` is required to split the data line.

    Examples::

      # Bar delimiter in header and data
      
      |  1.2  | hello there |     3 |
      |  2.4  | many words  |     7 |
      
      # Compact table having no delimiter and column positions specified as input

      1.2hello there3 
      2.4many words 7 

    This class is just a convenience wrapper around :class:`~asciitable.FixedWidth`
    but with ``header.start_line = None`` and ``data.start_line = 0``.

    See the :ref:`fixed_width_gallery` for specific usage examples.

    :param col_starts: list of start positions for each column (0-based counting)
    :param col_ends: list of end positions (inclusive) for each column
    :param delimiter_pad: padding around delimiter when writing (default = None)
    :param bookend: put the delimiter at start and end of line when writing (default = False)
    """
    def __init__(self, col_starts=None, col_ends=None, delimiter_pad=' ', bookend=True):
        FixedWidth.__init__(self, col_starts, col_ends,
                            delimiter_pad=delimiter_pad, bookend=bookend)
        self.header.start_line = None
        self.data.start_line = 0

        
class FixedWidthTwoLine(FixedWidth):
    """Read or write a fixed width table which has two header lines.  The first
    header line defines the column names and the second implicitly defines the
    column positions.  Examples::

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

    See the :ref:`fixed_width_gallery` for specific usage examples.

    :param position_line: row index of line that specifies position (default = 1)
    :param position_char: character used to write the position line (default = "-")
    :param delimiter_pad: padding around delimiter when writing (default = None)
    :param bookend: put the delimiter at start and end of line when writing (default = False)
    """
    def __init__(self, position_line=1, position_char='-', delimiter_pad=None, bookend=False):
        FixedWidth.__init__(self, delimiter_pad=delimiter_pad, bookend=bookend)
        self.header.position_line = position_line
        self.header.position_char = position_char
        self.data.start_line = position_line + 1
        self.header.splitter.delimiter = ' '
        self.data.splitter.delimiter = ' '

        
