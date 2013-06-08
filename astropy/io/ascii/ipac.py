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
## * Redistributions of source code must retain the above copyright
## notice, this list of conditions and the following disclaimer.
## * Redistributions in binary form must reproduce the above copyright
## notice, this list of conditions and the following disclaimer in the
## documentation and/or other materials provided with the distribution.
## * Neither the name of the Smithsonian Astrophysical Observatory nor the
## names of its contributors may be used to endorse or promote products
## derived from this software without specific prior written permission.
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
from collections import defaultdict
from textwrap import wrap
from warnings import warn

from . import core
from . import fixedwidth
from ...utils import OrderedDict



class IpacFormatErrorStrict(Exception):
    def __str__(self):
        return super(Exception, self).__str__() + '\nSee http://irsa.ipac.caltech.edu/applications/DDGEN/Doc/DBMSrestriction.html'


class IpacFormatError(Exception):
    def __str__(self):
        return super(Exception, self).__str__() + '\nSee http://irsa.ipac.caltech.edu/applications/DDGEN/Doc/ipac_tbl.html'
    

class Ipac(fixedwidth.FixedWidth):
    """Read or write an IPAC format table. See
    http://irsa.ipac.caltech.edu/applications/DDGEN/Doc/ipac_tbl.html::

        \\name=value
        \\ Comment
        | column1 | column2 | column3 | column4 | column5 |
        | double  | double  | int     | double  | char    |
        | unit    | unit    | unit    | unit    | unit    |
        | null    | null    | null    | null    | null    |
         2.0978    29.09056  73765     2.06000   B8IVpMnHg

    Or::

        |-----ra---|----dec---|---sao---|------v---|----sptype--------|
         2.09708    29.09056   73765     2.06000     B8IVpMnHg

    The comments and keywords defined in the header are available via the output
    table ``meta`` attribute::

        >>> from astropy.io import ascii
        >>> filename = os.path.join(ascii.__path__[0], 'tests/t/ipac.dat')
        >>> data = ascii.read(filename)
        >>> print data.meta['comments']
        ['This is an example of a valid comment']
        >>> for name, keyword in data.meta['keywords'].items():
        ... print name, keyword['value']
        ...
        intval 1
        floatval 2300.0
        date Wed Sp 20 09:48:36 1995
        key_continue IPAC keywords can continue across lines

    Parameters
    ----------
    definition : str, optional
    Specify the convention for characters in the data table that occur
    directly below the pipe (`|`) symbol in the header column definition:

        * 'ignore' - Any character beneath a pipe symbol is ignored (default)
        * 'right' - Character is associated with the column to the right
        * 'left' - Character is associated with the column to the left

    DBMS : bool, optional
        If true, this varifies that written tables adhere (semantically)
        to the `IPAC/DBMS <http://irsa.ipac.caltech.edu/applications/DDGEN/Doc/DBMSrestriction.html>`_
        definiton of IPAC tables. If 'False' it only checks for the (less strict)
        `IPAC <http://irsa.ipac.caltech.edu/applications/DDGEN/Doc/ipac_tbl.html>`_
        definition.
    """
    def __init__(self, definition='ignore', DBMS = True):
        super(fixedwidth.FixedWidth, self).__init__()
        self.header = IpacHeader(definition=definition)
        self.data = IpacData()
        self.data.header = self.header
        self.header.data = self.data
        self.header.DBMS = DBMS
        self.data.splitter.delimiter = ' '
        self.data.splitter.delimiter_pad = ''
        self.data.splitter.bookend = True

    def write(self, table):
        """Write ``table`` as list of strings.

        :param table: input table data (astropy.table.Table object)
        :returns: list of strings corresponding to ASCII table
        """
        # link information about the columns to the writer object (i.e. self)
        self.header.cols = table.cols
        self.data.cols = table.cols

        # Write header and data to lines list
        lines = []
        # Write meta information
        if 'comments' in table.meta:
            for comment in table.meta['comments']:
                if len(str(comment)) > 78:
                    warn('Comment string > 78 characters was automatically wrapped.',
                          UserWarning)
                for line in wrap(str(comment), 80, initial_indent='\\ ', subsequent_indent='\\ '):
                    lines.append(line)
        if 'keywords' in table.meta:
            keydict = table.meta['keywords']
            for keyword in keydict:
                try:
                    val = keydict[keyword]['value']
                    if isinstance(val, basestring): val = "'"+val+"'"
                    lines.append('\\{0}={1}'.format(keyword.strip(), val))
                    # meta is not standardized: Catch some common Errors.
                except TypeError:
                    pass

        #get header and data as strings to find width or each column
        for i, col in enumerate(table.cols):
            col.headwidth = max([len(vals[i]) for vals in self.header.str_vals()])
        # keep those because they take some time to make
        data_str_vals = self.data.str_vals()
        for i, col in enumerate(table.cols):
            col.width = max([len(vals[i]) for vals in data_str_vals])

        widths = [max(col.width, col.headwidth) for col in table.cols]

        # then write table
        self.header.write(lines, widths)
        self.data.write(lines, widths, data_str_vals)

        return lines



class IpacHeaderSplitter(core.BaseSplitter):
    '''Splitter for Ipac Headers.

    This splitter is similar its parent when reading, but supports a
    fixed width format (as required for Ipac table headers) for writing.
    '''
    process_line = None
    process_val = None
    delimiter = '|'
    delimiter_pad = ''
    skipinitialspace = False
    comment = r'\s*\\'
    write_comment = r'\\'
    col_starts = None
    col_ends = None

    def join(self, vals, widths):
        pad = self.delimiter_pad or ''
        delimiter = self.delimiter or ''
        padded_delim = pad + delimiter + pad
        bookend_left = delimiter + pad
        bookend_right = pad + delimiter
 
        vals = [' ' * (width - len(val)) + val for val, width in zip(vals, widths)]
        return bookend_left + padded_delim.join(vals) + bookend_right

class IpacHeader(fixedwidth.FixedWidthHeader):
    """IPAC table header"""
    splitter_class = IpacHeaderSplitter
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
       
        fixedwidth.FixedWidthHeader.__init__(self)
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
        Extract table-level comments and keywords for IPAC table. See:
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
                    # Strip leading/trailing quote. The spec says that a matched pair
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
            # Keywords and comments start with "\". Once the first non-slash
            # line is seen then bail out.
            if not line.startswith('\\'):
                break

            m = re_keyword.match(line)
            if m:
                name = m.group('name')
                val = process_keyword_value(m.group('value'))

                # IPAC allows for continuation keywords, e.g.
                # \SQL = 'WHERE '
                # \SQL = 'SELECT (25 column names follow in next row.)'
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
        Sets ``self.cols`` with the list of Columns. This list only includes the actual
        requested columns after filtering by the include_names and exclude_names
        attributes. See ``self.names`` for the full list.

        :param lines: list of table lines
        :returns: list of table Columns
        """
        header_lines = self.process_lines(lines) # generator returning valid header lines
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

    def str_vals(self):
        
        if self.DBMS:
            IpacFormatE = IpacFormatErrorStrict
        else:
            IpacFormatE = IpacFormatError
        
        namelist = [col.name for col in self.cols]
        if self.DBMS:
            countnamelist = defaultdict(int)
            for col in self.cols:
                countnamelist[col.name.lower()] += 1 
            doublenames = [x for x in countnamelist if countnamelist[x] > 1]
            if doublenames != []:
                raise IpacFormatE('IPAC DBMS tables are not case sensitive. This causes duplicate column names: {0}'.format(doublenames))

        for name in namelist:
            m = re.match('\w+', name)
            if m.end() != len(name):
                raise IpacFormatE('{0} - Only alphanumaric characters and _ are allowed in column names.'.format(name))
            if self.DBMS and not(name[0].isalpha() or (name[0] =='_')):
                raise IpacFormatE('Column names cannot stars with numbers: {}'.format(name))
            if self.DBMS:
                if name in ['x','y','z', 'X', 'Y','Z']:
                    raise IpacFormatE('{0} - x, y, z, X, Y, Z are reserved names and cannot be used as column names.'.format(name))
                if len(name) > 16:
                    raise IpacFormatE('{0} - Maximum length for column name is 16 characters'.format(name))
            else:
                if len(name) > 40:
                    raise IpacFormatE('{0} - Maximum length for column name is 40 characters.'.format(name))

        dtypelist = []
        unitlist = []
        for col in self.cols:
            if col.dtype.kind in ['i', 'u']:
                dtypelist.append('int')
            elif col.dtype.kind == 'f':
                dtypelist.append('float')
            else:
                dtypelist.append('char')
            if col.units is None:
                unitlist.append('unit')
            else:
                unitlist.append(str(col.units))
        nullist = [getattr(col, 'fill_value', 'null') for col in self.cols]
        return [namelist, dtypelist, unitlist, nullist]

    def write(self, lines, widths):
        '''Write header.

        The width of each column is determined in IpacData.write. Writing the header
        must be delayed until that time.
        This function is called from data, once the widht information is
        available.'''

        for vals in self.str_vals():
            lines.append(self.splitter.join(vals, widths))
        return lines

class IpacData(fixedwidth.FixedWidthData):
    """IPAC table data reader"""
    comment = r'[|\\]'


    def str_vals(self):
        '''return str vals for each in the table'''
        vals_list = []
        # just to make sure
        self._set_col_formats()
        col_str_iters = [col.iter_str_vals() for col in self.cols]
        for vals in core.izip(*col_str_iters):
            vals_list.append(vals)

        return vals_list


    def write(self, lines, widths, vals_list):
        """ IPAC writer, modified from FixedWidth writer """
        for vals in vals_list:
            lines.append(self.splitter.join(vals, widths))
        return lines
