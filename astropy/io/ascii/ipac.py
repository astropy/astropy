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

    Notes
    -----
    Caveats:

    * Data type, Units, and Null value specifications are ignored.
    * Keywords are ignored.
    * The IPAC spec requires the first two header lines but this reader only
      requires the initial column name definition line

    Overcoming these limitations would not be difficult, code contributions
    welcome from motivated users.

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
    comment = r'\\'
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
            raise ValueError('At least one header line beginning and ending with delimiter required')
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

        # Since the splitter returns only the actual requested columns, at this
        # point set self.n_data_cols to be the number of requested columns.  This
        # gets used later to validate the data as it gets read and split.
        self.n_data_cols = len(self.cols)

    def str_vals(self):

        if self.DBMS:
            IpacFormatE = IpacFormatErrorDBMS
        else:
            IpacFormatE = IpacFormatError

        namelist = [col.name for col in self.cols]
        if self.DBMS:
            countnamelist = defaultdict(int)
            for col in self.cols:
                countnamelist[col.name.lower()] += 1
            doublenames = [x for x in countnamelist if countnamelist[x] > 1]
            if doublenames != []:
                raise IpacFormatE('IPAC DBMS tables are not case sensitive. '
                                  'This causes duplicate column names: {0}'.format(doublenames))

        for name in namelist:
            m = re.match('\w+', name)
            if m.end() != len(name):
                raise IpacFormatE('{0} - Only alphanumaric characters and _ '
                                  'are allowed in column names.'.format(name))
            if self.DBMS and not(name[0].isalpha() or (name[0] == '_')):
                raise IpacFormatE('Column name cannot start with numbers: {}'.format(name))
            if self.DBMS:
                if name in ['x', 'y', 'z', 'X', 'Y', 'Z']:
                    raise IpacFormatE('{0} - x, y, z, X, Y, Z are reserved names and '
                                      'cannot be used as column names.'.format(name))
                if len(name) > 16:
                    raise IpacFormatE(
                        '{0} - Maximum length for column name is 16 characters'.format(name))
            else:
                if len(name) > 40:
                    raise IpacFormatE(
                        '{0} - Maximum length for column name is 40 characters.'.format(name))

        dtypelist = []
        unitlist = []
        for col in self.cols:
            if col.dtype.kind in ['i', 'u']:
                dtypelist.append('long')
            elif col.dtype.kind == 'f':
                dtypelist.append('double')
            else:
                dtypelist.append('char')
            if col.unit is None:
                unitlist.append('')
            else:
                unitlist.append(str(col.unit))
        nullist = [getattr(col, 'fill_value', 'null') for col in self.cols]
        return [namelist, dtypelist, unitlist, nullist]

    def write(self, lines, widths):
        '''Write header.

        The width of each column is determined in IpacData.write. Writing the header
        must be delayed until that time.
        This function is called from data, once the width information is
        available.'''

        for vals in self.str_vals():
            lines.append(self.splitter.join(vals, widths))
        return lines


class IpacData(fixedwidth.FixedWidthData):
    """IPAC table data reader"""
    splitter_class = fixedwidth.FixedWidthSplitter
    comment = r'[|\\]'
