"""Asciitable: an extensible ASCII table reader and writer.

basic.py:
  Basic table read / write functionality for simple character
  delimited files with various options for column header definition.

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
import asciitable.core as core
from asciitable.core import io, next, izip, any

class Basic(core.BaseReader):
    """Read a character-delimited table with a single header line at the top
    followed by data lines to the end of the table.  Lines beginning with # as
    the first non-whitespace character are comments.  This reader is highly
    configurable.
    ::

        rdr = asciitable.get_reader(Reader=asciitable.Basic)
        rdr.header.splitter.delimiter = ' '
        rdr.data.splitter.delimiter = ' '
        rdr.header.start_line = 0
        rdr.data.start_line = 1
        rdr.data.end_line = None
        rdr.header.comment = r'\s*#'
        rdr.data.comment = r'\s*#'

    Example table::
    
      # Column definition is the first uncommented line
      # Default delimiter is the space character.
      apples oranges pears

      # Data starts after the header column definition, blank lines ignored
      1 2 3
      4 5 6
    """
    def __init__(self):
        core.BaseReader.__init__(self)
        self.header.splitter.delimiter = ' '
        self.data.splitter.delimiter = ' '
        self.header.start_line = 0
        self.data.start_line = 1
        self.header.comment = r'\s*#'
        self.header.write_comment = '# '
        self.data.comment = r'\s*#'
        self.data.write_comment = '# '

BasicReader = Basic

class NoHeader(BasicReader):
    """Read a table with no header line.  Columns are autonamed using
    header.auto_format which defaults to "col%d".  Otherwise this reader
    the same as the :class:`Basic` class from which it is derived.  Example::

      # Table data
      1 2 "hello there"
      3 4 world
    """
    def __init__(self):
        BasicReader.__init__(self)
        self.header.start_line = None
        self.data.start_line = 0

NoHeaderReader = NoHeader

class CommentedHeaderHeader(core.BaseHeader):
    """Header class for which the column definition line starts with the
    comment character.  See the :class:`CommentedHeader` class  for an example.
    """
    def process_lines(self, lines):
        """Return only lines that start with the comment regexp.  For these
        lines strip out the matching characters."""
        re_comment = re.compile(self.comment)
        for line in lines:
            match = re_comment.match(line)
            if match:
                yield line[match.end():]

    def write(self, lines):
        lines.append(self.write_comment + self.splitter.join([x.name for x in self.cols]))

class CommentedHeader(core.BaseReader):
    """Read a file where the column names are given in a line that begins with the
    header comment character.  The default delimiter is the <space> character.::

      # col1 col2 col3
      # Comment line
      1 2 3
      4 5 6
    """
    def __init__(self):
        core.BaseReader.__init__(self)
        self.header = CommentedHeaderHeader()
        self.header.data = self.data
        self.data.header = self.header
        self.header.splitter.delimiter = ' '
        self.data.splitter.delimiter = ' '
        self.header.start_line = 0
        self.data.start_line = 0
        self.header.comment = r'\s*#'
        self.header.write_comment = '# '
        self.data.comment = r'\s*#'
        self.data.write_comment = '# '

CommentedHeaderReader = CommentedHeader

class Tab(BasicReader):
    """Read a tab-separated file.  Unlike the :class:`Basic` reader, whitespace is 
    not stripped from the beginning and end of lines.  By default whitespace is
    still stripped from the beginning and end of individual column values.

    Example::

      col1 <tab> col2 <tab> col3
      # Comment line
      1 <tab> 2 <tab> 5
    """
    def __init__(self):
        BasicReader.__init__(self)
        self.header.splitter.delimiter = '\t'
        self.data.splitter.delimiter = '\t'
        # Don't strip line whitespace since that includes tabs
        self.header.splitter.process_line = None  
        self.data.splitter.process_line = None
        # Don't strip data value whitespace since that is significant in TSV tables
        self.data.splitter.process_val = None
        self.data.splitter.skipinitialspace = False

TabReader = Tab

class Rdb(TabReader):
    """Read a tab-separated file with an extra line after the column definition
    line.  The RDB format meets this definition.  Example::

      col1 <tab> col2 <tab> col3
      N <tab> S <tab> N
      1 <tab> 2 <tab> 5

    In this reader the second line is just ignored.
    """
    def __init__(self):
        TabReader.__init__(self)
        self.header = RdbHeader()
        self.header.start_line = 0
        self.header.comment = r'\s*#'
        self.header.write_comment = '# '
        self.header.splitter.delimiter = '\t'
        self.header.splitter.process_line = None
        self.header.data = self.data
        self.data.header = self.header
        self.data.start_line = 2

RdbReader = Rdb

class RdbHeader(core.BaseHeader):
    col_type_map = {'n': core.NumType,
                    's': core.StrType}

    def get_type_map_key(self, col):
        return col.raw_type[-1]

    def get_cols(self, lines):
        """Initialize the header Column objects from the table ``lines``.
        
        This is a specialized get_cols for the RDB type:
        Line 0: RDB col names
        Line 1: RDB col definitions
        Line 2+: RDB data rows

        :param lines: list of table lines
        :returns: None
        """
        header_lines = self.process_lines(lines)   # this is a generator
        header_vals_list = [hl for _, hl in zip(range(2), self.splitter(header_lines))]
        if len(header_vals_list) != 2:
            raise ValueError('RDB header requires 2 lines')
        self.names, raw_types = header_vals_list

        if len(self.names) != len(raw_types):
            raise ValueError('RDB header mismatch between number of column names and column types')
        
        if any(not re.match(r'\d*(N|S)$', x, re.IGNORECASE) for x in raw_types):
            raise ValueError('RDB types definitions do not all match [num](N|S): %s' % raw_types)
        
        self._set_cols_from_names()
        for col, raw_type in zip(self.cols, raw_types):
            col.raw_type = raw_type
            col.type = self.get_col_type(col)

    def write(self, lines):
        lines.append(self.splitter.join([x.name for x in self.cols]))
        rdb_types = []
        for col in self.cols:
            if issubclass(col.type, core.NumType):
                rdb_types.append('N')
            else:
                rdb_types.append('S')
        lines.append(self.splitter.join(rdb_types))

