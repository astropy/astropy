# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""An extensible ASCII table reader and writer.

ipac.py:
  Classes to read IPAC table format

:Copyright: Smithsonian Astrophysical Observatory (2011)
:Author: Tom Aldcroft (aldcroft@head.cfa.harvard.edu)
"""

from __future__ import absolute_import, division, print_function

import re
from collections import defaultdict
from textwrap import wrap
from warnings import warn

from ...extern import six
from ...extern.six.moves import zip

from . import core
from . import fixedwidth
from . import basic
from ...utils import OrderedDict
from ...utils.exceptions import AstropyUserWarning
from ...table.pprint import _format_funcs, _auto_format_func


class IpacFormatErrorDBMS(Exception):
    def __str__(self):
        return '{0}\nSee {1}'.format(
            super(Exception, self).__str__(),
            'http://irsa.ipac.caltech.edu/applications/DDGEN/Doc/DBMSrestriction.html')


class IpacFormatError(Exception):
    def __str__(self):
        return '{0}\nSee {1}'.format(
            super(Exception, self).__str__(),
            'http://irsa.ipac.caltech.edu/applications/DDGEN/Doc/ipac_tbl.html')


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
    definition='ignore'
    start_line = None

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
                if name in keywords and isinstance(val, six.string_types):
                    prev_val = keywords[name]['value']
                    if isinstance(prev_val, six.string_types):
                        val = prev_val + val

                keywords[name] = {'value': val}
            else:
                # Comment is required to start with "\ "
                if line.startswith('\\ '):
                    val = line[2:].strip()
                    if val:
                        table_meta['comments'].append(val)

    def get_cols(self, lines):
        """Initialize the header Column objects from the table ``lines``.

        Based on the previously set Header attributes find or create the column names.
        Sets ``self.cols`` with the list of Columns.

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
            col = core.Column(name=name.strip(' -'))
            col.start = start
            col.end = start + len(name)
            if len(header_vals) > 1:
                col.raw_type = header_vals[1][i].strip(' -')
                col.type = self.get_col_type(col)
            if len(header_vals) > 2:
                col.unit = header_vals[2][i].strip()  # Can't strip dashes here
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

        self.names = [x.name for x in cols]
        self.cols = cols



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
        nullist = []
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
            null = col.fill_values[core.masked]
            try:
                format_func = _format_funcs.get(col.format, _auto_format_func)
                nullist.append((format_func(col.format, null)).strip())
            except:
                # It is possible that null and the column values have different
                # data types (e.g. number und null = 'null' (i.e. a string).
                # This could cause all kinds of exceptions, so a catch all
                # block is needed here
                nullist.append(str(null).strip())

        return [namelist, dtypelist, unitlist, nullist]

    def write(self, lines, widths):
        '''Write header.

        The width of each column is determined in Ipac.write. Writing the header
        must be delayed until that time.
        This function is called from there, once the width information is
        available.'''

        for vals in self.str_vals():
            lines.append(self.splitter.join(vals, widths))
        return lines


class IpacDataSplitter(fixedwidth.FixedWidthSplitter):
    delimiter = ' '
    delimiter_pad = ''
    bookend = True


class IpacData(fixedwidth.FixedWidthData):
    """IPAC table data reader"""
    comment = r'[|\\]'
    start_line = 0
    splitter_class = IpacDataSplitter
    fill_values = [(core.masked, 'null')]

    def write(self, lines, widths, vals_list):
        """ IPAC writer, modified from FixedWidth writer """
        for vals in vals_list:
            lines.append(self.splitter.join(vals, widths))
        return lines


class Ipac(basic.Basic):
    """Read or write an IPAC format table.  See
    http://irsa.ipac.caltech.edu/applications/DDGEN/Doc/ipac_tbl.html::

      \\name=value
      \\ Comment
      |  column1 |   column2 | column3 | column4  |    column5    |
      |  double  |   double  |   int   |   double |     char      |
      |  unit    |   unit    |   unit  |    unit  |     unit      |
      |  null    |   null    |   null  |    null  |     null      |
       2.0978     29.09056    73765     2.06000    B8IVpMnHg

    Or::

      |-----ra---|----dec---|---sao---|------v---|----sptype--------|
        2.09708   29.09056     73765   2.06000    B8IVpMnHg

    The comments and keywords defined in the header are available via the output
    table ``meta`` attribute::

      >>> import os
      >>> from astropy.io import ascii
      >>> filename = os.path.join(ascii.__path__[0], 'tests/t/ipac.dat')
      >>> data = ascii.read(filename)
      >>> print(data.meta['comments'])
      ['This is an example of a valid comment']
      >>> for name, keyword in data.meta['keywords'].items():
      ...     print(name, keyword['value'])
      ...
      intval 1
      floatval 2300.0
      date Wed Sp 20 09:48:36 1995
      key_continue IPAC keywords can continue across lines

    Note that there are different conventions for characters occuring below the
    position of the ``|`` symbol in IPAC tables. By default, any character
    below a ``|`` will be ignored (since this is the current standard),
    but if you need to read files that assume characters below the ``|``
    symbols belong to the column before or after the ``|``, you can specify
    ``definition='left'`` or ``definition='right'`` respectively when reading
    the table (the default is ``definition='ignore'``). The following examples
    demonstrate the different conventions:

    * ``definition='ignore'``::

        |   ra  |  dec  |
        | float | float |
          1.2345  6.7890

    * ``definition='left'``::

        |   ra  |  dec  |
        | float | float |
           1.2345  6.7890

    * ``definition='right'``::

        |   ra  |  dec  |
        | float | float |
        1.2345  6.7890

    IPAC tables can specify a null value in the header that is shown in place
    of missing or bad data. On writing, this value defaults to ``null``.
    To specify a different null value, use the ``fill_values`` option to
    replace masked values with a string or number of your choice as
    described in :ref:`io_ascii_write_parameters`::

        >>> from astropy.io.ascii import masked
        >>> fill = [(masked, 'N/A', 'ra'), (masked, -999, 'sptype')]
        >>> ascii.write(data, format='ipac', fill_values=fill)
        \ This is an example of a valid comment
        ...
        |          ra|         dec|      sai|          v2|            sptype|
        |      double|      double|     long|      double|              char|
        |        unit|        unit|     unit|        unit|              ergs|
        |         N/A|        null|     null|        null|              -999|
                  N/A     29.09056      null         2.06               -999 
         2345678901.0 3456789012.0 456789012 4567890123.0 567890123456789012 


    Parameters
    ----------
    definition : str, optional
        Specify the convention for characters in the data table that occur
        directly below the pipe (``|``) symbol in the header column definition:

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
    _format_name = 'ipac'
    _io_registry_format_aliases = ['ipac']
    _io_registry_can_write = True
    _description = 'IPAC format table'

    data_class = IpacData
    header_class = IpacHeader

    def __init__(self, definition='ignore', DBMS=False):
        super(Ipac, self).__init__()
        # Usually the header is not defined in __init__, but here it need a keyword
        if definition in ['ignore', 'left', 'right']:
            self.header.ipac_definition = definition
        else:
            raise ValueError("definition should be one of ignore/left/right")
        self.header.DBMS = DBMS

    def write(self, table):
        """Write ``table`` as list of strings.

        :param table: input table data (astropy.table.Table object)
        :returns: list of strings corresponding to ASCII table
        """
        # Set a default null value for all columns by adding at the end, which
        # is the position with the lowest priority.
        # We have to do it this late, because the fill_value
        # defined in the class can be overwritten by ui.write
        self.data.fill_values.append((core.masked, 'null'))

        # Check column names before altering
        self.header.cols = list(six.itervalues(table.columns))
        self.header.check_column_names(self.names, self.strict_names, self.guessing)

        core._apply_include_exclude_names(table, self.names, self.include_names, self.exclude_names)

        # Now use altered columns
        new_cols = list(six.itervalues(table.columns))
        # link information about the columns to the writer object (i.e. self)
        self.header.cols = new_cols
        self.data.cols = new_cols

        # Write header and data to lines list
        lines = []
        # Write meta information
        if 'comments' in table.meta:
            for comment in table.meta['comments']:
                if len(str(comment)) > 78:
                    warn('Comment string > 78 characters was automatically wrapped.',
                         AstropyUserWarning)
                for line in wrap(str(comment), 80, initial_indent='\\ ', subsequent_indent='\\ '):
                    lines.append(line)
        if 'keywords' in table.meta:
            keydict = table.meta['keywords']
            for keyword in keydict:
                try:
                    val = keydict[keyword]['value']
                    lines.append('\\{0}={1!r}'.format(keyword.strip(), val))
                    # meta is not standardized: Catch some common Errors.
                except TypeError:
                    pass

        # Usually, this is done in data.write, but since the header is written
        # first, we need that here.
        self.data._set_fill_values(self.data.cols)

        # get header and data as strings to find width of each column
        for i, col in enumerate(table.columns.values()):
            col.headwidth = max([len(vals[i]) for vals in self.header.str_vals()])
        # keep data_str_vals because they take some time to make
        data_str_vals = []
        col_str_iters = self.data.str_vals()
        for vals in zip(*col_str_iters):
            data_str_vals.append(vals)

        for i, col in enumerate(table.columns.values()):
            # FIXME: In Python 3.4, use max([], default=0).
            # See: https://docs.python.org/3/library/functions.html#max
            if data_str_vals:
                col.width = max([len(vals[i]) for vals in data_str_vals])
            else:
                col.width = 0

        widths = [max(col.width, col.headwidth) for col in table.columns.values()]
        # then write table
        self.header.write(lines, widths)
        self.data.write(lines, widths, data_str_vals)

        return lines
