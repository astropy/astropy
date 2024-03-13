# Licensed under a 3-clause BSD style license
"""
:Author: Simon Gibbons (simongibbons@gmail.com).
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


class RST(FixedWidth):
    """reStructuredText simple format table.

    See: https://docutils.sourceforge.io/docs/ref/rst/restructuredtext.html#simple-tables

    Example::

      >>> from astropy.table import QTable
      >>> import astropy.units as u
      >>> import sys
      >>> tbl = QTable({"wave": [350, 950] * u.nm, "response": [0.7, 1.2] * u.count})
      >>> tbl.write(sys.stdout,  format="ascii.rst")
      ===== ========
       wave response
      ===== ========
      350.0      0.7
      950.0      1.2
      ===== ========

    Like other fixed-width formats, when writing a table you can provide ``header_rows``
    to specify a list of table rows to output as the header.  For example::

      >>> tbl.write(sys.stdout,  format="ascii.rst", header_rows=['name', 'unit'])
      ===== ========
       wave response
         nm       ct
      ===== ========
      350.0      0.7
      950.0      1.2
      ===== ========

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
        idx = len(self.header.header_rows)
        lines = [lines[idx]] + lines + [lines[idx]]
        return lines

    def read(self, table):
        self.data.start_line = 2 + len(self.header.header_rows)
        return super().read(table)
