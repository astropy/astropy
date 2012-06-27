"""Asciitable: an extensible ASCII table reader and writer.

besancon.py:
  Classes to read Besancon galaxy model

:Author: Adam Ginsburg (adam.g.ginsburg@gmail.com)
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

from . import core
import re

class Besancon(core.BaseReader):
    """Read a Besancon galaxy model.  See:
    http://model.obs-besancon.fr/::

    """
    def __init__(self):
        core.BaseReader.__init__(self)
        self.header = BesanconHeader()
        self.data = BesanconData()

    def write(self, table=None):
        """Not available for the Besancon class (raises NotImplementedError)"""
        raise NotImplementedError


class BesanconHeader(core.BaseHeader):
    """Besancon table header"""
    comment = r'#'
    splitter_class = core.BaseSplitter

    def __init__(self):
        self.splitter = self.__class__.splitter_class()
        self.splitter.process_line = None
        self.splitter.process_val = None
        self.splitter.delimiter = ' '

    def process_lines(self, lines):
        """Generator to yield Besancon header lines, i.e. those starting and ending with
        delimiter character."""
        delim = self.splitter.delimiter
        for linenum,line in enumerate(lines):
            if linenum >= 81:
                break
            yield line.strip(delim)

    def get_cols(self, lines):
        """Initialize the header Column objects from the table ``lines``.

        Based on the previously set Header attributes find or create the column names.
        Sets ``self.cols`` with the list of Columns.  This list only includes the actual
        requested columns after filtering by the include_names and exclude_names
        attributes.  See ``self.names`` for the full list.

        :param lines: list of table lines
        :returns: list of table Columns
        """
        header_lines = list(self.process_lines(lines))  # generator returning valid header lines
        header_vals = [vals for vals in self.splitter(header_lines)]
        #if len(header_vals) == 0:
        #    raise ValueError('At least one header line beginning and ending with delimiter required')
        #elif len(header_vals) > 4:
        #    raise ValueError('More than four header lines were found')

        # Generate column definitions
        cols = (header_lines[-1].split())
        start = 1
        for ii, name in enumerate(cols):
            col = core.Column(name=name, index=ii)
            #col.start = start
            #col.end = start + len(name)
            #if len(header_vals) > 1:
            #    col.raw_type = header_vals[1][i].strip(' -')
            #    col.type = self.get_col_type(col)
            #if len(header_vals) > 2:
            #    col.units = header_vals[2][i].strip() # Can't strip dashes here
            #if len(header_vals) > 3:
            #    null = header_vals[3][i].strip()
            #    if null.lower() != 'null':
            #        col.null = null  # Can't strip dashes here
            fillval = 'nan'
            col.type = core.FloatType
            col.null = 'null'
            self.data.fill_values.append((col.null, fillval, col.name))
            cols.append(col)

        # Standard column name filtering (include or exclude names)
        self.names = [x.name for x in cols]
        names = set(self.names)
        if self.include_names is not None:
            names.intersection_update(self.include_names)
        if self.exclude_names is not None:
            names.difference_update(self.exclude_names)

        # Generate final list of cols and re-index the cols because the
        # FixedWidthSplitter does NOT return the ignored cols (as is the
        # case for typical delimiter-based splitters)
        self.cols = [x for x in cols if x.name in names]
        for i, col in enumerate(self.cols):
            col.index = i

class BesanconData(core.BaseData):
    """Besacon table data reader"""
    splitter_class = core.BaseSplitter
    comment = r'#'

    def process_lines(self, lines):
        """Strip out comment lines and blank lines from list of ``lines``

        :param lines: all lines in table
        :returns: list of lines
        """
        nonblank_lines = (x for x in lines if x.strip())
        if self.comment:
            re_comment = re.compile(self.comment)
            lines = [x for x in nonblank_lines if not re_comment.match(x)]
        else:
            lines = [x for x in nonblank_lines]
        return lines[:-6]
