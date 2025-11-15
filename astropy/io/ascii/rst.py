# Licensed under a 3-clause BSD style license
"""
:Author: Simon Gibbons (simongibbons@gmail.com)
"""


from .core import DefaultSplitter
from .fixedwidth import (
    FixedWidth,
    FixedWidthData,
    FixedWidthHeader,
    FixedWidthTwoLineDataSplitter,
)


class SimpleRSTHeader(FixedWidthHeader):
    position_line = 0
    start_line = 1
    splitter_class = DefaultSplitter
    position_char = "="

    def get_fixedwidth_params(self, line):
        vals, starts, ends = super().get_fixedwidth_params(line)
        # The right hand column can be unbounded
        ends[-1] = None
        return vals, starts, ends


class SimpleRSTData(FixedWidthData):
    end_line = -1
    splitter_class = FixedWidthTwoLineDataSplitter

    @property
    def start_line(self):
        """Calculate start line based on number of header rows."""
        # RST format has: position_line (0), header rows starting at line 1,
        # position_line repeat, then data
        # For default (single header row): position_line=0, name at 1, position_line=2, data=3
        header_rows = getattr(self, "header_rows", ["name"])
        return len(header_rows) + 2


class RST(FixedWidth):
    """reStructuredText simple format table.

    See: https://docutils.sourceforge.io/docs/ref/rst/restructuredtext.html#simple-tables

    Example::

        ==== ===== ======
        Col1  Col2  Col3
        ==== ===== ======
          1    2.3  Hello
          2    4.5  Worlds
        ==== ===== ======

    Currently there is no support for reading tables which utilize continuation lines,
    or for ones which define column spans through the use of an additional
    line of dashes in the header.

    """

    _format_name = "rst"
    _description = "reStructuredText simple table"
    data_class = SimpleRSTData
    header_class = SimpleRSTHeader

    def __init__(self, header_rows=None):
        super().__init__(delimiter_pad=None, bookend=False, header_rows=header_rows)

    def write(self, lines):
        lines = super().write(lines)
        # Determine the position line index (it comes after all header rows)
        header_rows = self.header.header_rows if hasattr(self.header, 'header_rows') else ["name"]
        position_line_index = len(header_rows)
        # RST format: separator line, header rows, separator line, data, separator line
        separator = lines[position_line_index]
        lines = [separator] + lines + [separator]
        return lines
