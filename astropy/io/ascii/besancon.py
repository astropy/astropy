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
from . import fixedwidth
import re

__all__ = ['BesanconFixed']

class BesanconFixed(fixedwidth.FixedWidth):
    """
    Read data from a Besancon galactic model table file.  Assumes a
    fixed-length header; it is possible that different parameters in the model
    will result in different length headers
    """

    def __init__(self, col_starts=None, col_ends=None, delimiter_pad=' ', bookend=True,
            header_line=80):
        super(BesanconFixed,self).__init__(col_starts=col_starts,
                col_ends=col_ends, delimiter_pad=delimiter_pad,
                bookend=bookend)

        self.header = BesanconFixedWidthHeader()

        self.data.header = self.header
        self.header.data = self.data

        # the first header line is at line 80, but fixedwidth ignores blank lines?
        # also, the header's widths do not match the data widths
        self.header.start_line = header_line
        self.header.position_line = header_line
        # this is a silly hack, I have NO idea why it's 71 and not 81
        self.data.start_line = header_line - 9 #71
        self.data.end_line = -6 # there are 6 lines at the end that should be excluded

        self.header.splitter.delimiter = ' '
        self.header.splitter.delimiter = ' '
        self.header.start_line = 0
        self.header.comment = r'\s*#'
        self.header.write_comment = '# '
        self.header.col_starts = col_starts
        self.header.col_ends = col_ends

class BesanconFixedWidthHeader(fixedwidth.FixedWidthHeader):
    def get_cols(self, lines):
        super(BesanconFixedWidthHeader,self).get_cols(lines)
        self.names = lines[self.position_line].split()

        data_lines = self.data.process_lines(lines)
        vals, starts, ends = self.get_fixedwidth_params(data_lines[0])
        self.n_data_cols = len(self.cols)

        self._set_cols_from_names()
        
        # Set column start and end positions.  Also re-index the cols because
        # the FixedWidthSplitter does NOT return the ignored cols (as is the
        # case for typical delimiter-based splitters)
        for i, col in enumerate(self.cols):
            col.start = starts[col.index]
            col.end = ends[col.index]
            col.index = i

            
class BesanconFixedWidthData(fixedwidth.FixedWidthData):
    def process_lines(self, lines):
        """Strip out comment lines from list of ``lines``
        (unlike the normal process_lines, does NOT exclude blank lines)

        :param lines: all lines in table
        :returns: list of lines
        """
        if self.comment:
            re_comment = re.compile(self.comment)
            return [x for x in lines if not re_comment.match(x)]
        else:
            return [x for x in lines]
