# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""An extensible ASCII table reader and writer.

ipac.py:
  Classes to read IPAC table format

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
from ...utils import OrderedDict

from . import core
from . import fixedwidth

class Ipac(core.BaseReader):
    """Read an IPAC format table.  See
    http://irsa.ipac.caltech.edu/applications/DDGEN/Doc/ipac_tbl.html::

      \\name=value
      \\ Comment
      |  column1 |  column2 | column3 | column4  |    column5       |
      |  double  |  double  |   int   |   double |     char         |
      |   unit   |   unit   |   unit  |    unit  |     unit         |
      |   null   |   null   |   null  |    null  |     null         |
       2.0978     29.09056   73765     2.06000    B8IVpMnHg

    Or::

      |-----ra---|----dec---|---sao---|------v---|----sptype--------|
        2.09708   29.09056     73765    2.06000   B8IVpMnHg

    Parameters
    ----------
    definition : str, optional
        Specify the convention for characters in the data table that occur
        directly below the pipe (`|`) symbol in the header column definition:

          * 'ignore' - Any character beneath a pipe symbol is ignored (default)
          * 'right' - Character is associated with the column to the right
          * 'left' - Character is associated with the column to the left

    The comments and keywords defined in the header are available via the output
    table ``meta`` attribute::

      data = ascii.read('t/ipac.dat')
      print data.meta['comments']
      for name, keyword in data.meta['keywords'].items():
          print name, keyword['value']
    """
    def __init__(self, definition='ignore'):
        core.BaseReader.__init__(self)
        self.header = IpacHeader(definition=definition)
        self.data = IpacData()

    def write(self, table=None):
        """Not available for the Ipac class (raises NotImplementedError)"""
        raise NotImplementedError


class IpacHeader(core.BaseHeader):
    """IPAC table header"""
    comment = '\\'
    splitter_class = core.BaseSplitter
    col_type_map = {'int': core.IntType,
                    'integer': core.IntType,
                    'long': core.IntType,
                    'double': core.FloatType,
                    'float': core.FloatType,
                    'real': core.FloatType,
                    'char': core.StrType,
                    'date': core.StrType,
                    'i': core.IntType,
                    'l': core.IntType,
                    'd': core.FloatType,
                    'f': core.FloatType,
                    'r': core.FloatType,
                    'c': core.StrType}

    def __init__(self, definition='ignore'):
        self.splitter = self.__class__.splitter_class()
        self.splitter.process_line = None
        self.splitter.process_val = None
        self.splitter.delimiter = '|'
        if definition in ['ignore', 'left', 'right']:
            self.ipac_definition = definition
        else:
            raise ValueError("definition should be one of ignore/left/right")

    def process_lines(self, lines):
        """Generator to yield IPAC header lines, i.e. those starting and ending with
        delimiter character."""
        delim = self.splitter.delimiter
        for line in lines:
            if line.startswith(delim) and line.endswith(delim):
                yield line.strip(delim)

    def update_meta(self, lines, meta):
        """
        Extract table-level comments and keywords for IPAC table.  See:
        http://irsa.ipac.caltech.edu/applications/DDGEN/Doc/ipac_tbl.html#kw
        """
        def process_keyword_value(val):
            """
            Take a string value and convert to float, int or str, and strip quotes
            as needed.
            """
            val = val.strip()
            try:
                val = int(val)
            except:
                try:
                    val = float(val)
                except:
                    # Strip leading/trailing quote.  The spec says that a matched pair
                    # of quotes is required, but this code will allow a non-quoted value.
                    for quote in ('"', "'"):
                        if val.startswith(quote) and val.endswith(quote):
                            val = val[1:-1]
                            break
            return val

        table_meta = meta['table']
        table_meta['comments'] = []
        table_meta['keywords'] = OrderedDict()
        keywords = table_meta['keywords']

        re_keyword = re.compile(r'\\'
                                r'(?P<name> \w+)'
                                r'\s* = (?P<value> .+) $',
                                re.VERBOSE)
        for line in lines:
            # Keywords and comments start with "\".  Once the first non-slash
            # line is seen then bail out.
            if not line.startswith('\\'):
                break

            m = re_keyword.match(line)
            if m:
                name = m.group('name')
                val = process_keyword_value(m.group('value'))

                # IPAC allows for continuation keywords, e.g.
                # \SQL     = 'WHERE '
                # \SQL     = 'SELECT (25 column names follow in next row.)'
                if name in keywords and isinstance(val, basestring):
                    prev_val = keywords[name]['value']
                    if isinstance(prev_val, basestring):
                        val = prev_val + val

                table_meta['keywords'][name] = {'value': val}
            else:
                # Comment is required to start with "\ "
                if line.startswith('\\ '):
                    val = line[2:].strip()
                    if val:
                        table_meta['comments'].append(val)

    def get_cols(self, lines):
        """Initialize the header Column objects from the table ``lines``.

        Based on the previously set Header attributes find or create the column names.
        Sets ``self.cols`` with the list of Columns.  This list only includes the actual
        requested columns after filtering by the include_names and exclude_names
        attributes.  See ``self.names`` for the full list.

        :param lines: list of table lines
        :returns: list of table Columns
        """
        header_lines = self.process_lines(lines)  # generator returning valid header lines
        header_vals = [vals for vals in self.splitter(header_lines)]
        if len(header_vals) == 0:
            raise ValueError('At least one header line beginning and ending with '
                             'delimiter required')
        elif len(header_vals) > 4:
            raise ValueError('More than four header lines were found')

        # Generate column definitions
        cols = []
        start = 1
        for i, name in enumerate(header_vals[0]):
            col = core.Column(name=name.strip(' -'), index=i)
            col.start = start
            col.end = start + len(name)
            if len(header_vals) > 1:
                col.raw_type = header_vals[1][i].strip(' -')
                col.type = self.get_col_type(col)
            if len(header_vals) > 2:
                col.units = header_vals[2][i].strip() # Can't strip dashes here
            if len(header_vals) > 3:
                # The IPAC null value corresponds to the io.ascii bad_value.
                # In this case there isn't a fill_value defined, so just put
                # in the minimal entry that is sure to convert properly to the
                # required type.
                #
                # Strip spaces but not dashes (not allowed in NULL row per
                # https://github.com/astropy/astropy/issues/361)
                null = header_vals[3][i].strip()
                fillval = '' if issubclass(col.type, core.StrType) else '0'
                self.data.fill_values.append((null, fillval, col.name))
            start = col.end + 1
            cols.append(col)

            # Correct column start/end based on definition
            if self.ipac_definition == 'right':
                col.start -= 1
            elif self.ipac_definition == 'left':
                col.end += 1

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

class IpacData(core.BaseData):
    """IPAC table data reader"""
    splitter_class = fixedwidth.FixedWidthSplitter
    comment = r'[|\\]'
