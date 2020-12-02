# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""An extensible ASCII table reader and writer.

basic.py:
  Basic table read / write functionality for simple character
  delimited files with various options for column header definition.

:Copyright: Smithsonian Astrophysical Observatory (2011)
:Author: Tom Aldcroft (aldcroft@head.cfa.harvard.edu)
"""


import re

from . import core


class BasicHeader(core.BaseHeader):
    """
    Basic table Header Reader

    Set a few defaults for common ascii table formats
    (start at line 0, comments begin with ``#`` and possibly white space)
    """
    start_line = 0
    comment = r'\s*#'
    write_comment = '# '


class BasicData(core.BaseData):
    """
    Basic table Data Reader

    Set a few defaults for common ascii table formats
    (start at line 1, comments begin with ``#`` and possibly white space)
    """
    start_line = 1
    comment = r'\s*#'
    write_comment = '# '


class Basic(core.BaseReader):
    r"""Character-delimited table with a single header line at the top.

    Lines beginning with a comment character (default='#') as the first
    non-whitespace character are comments.

    Example table::

      # Column definition is the first uncommented line
      # Default delimiter is the space character.
      apples oranges pears

      # Data starts after the header column definition, blank lines ignored
      1 2 3
      4 5 6
    """
    _format_name = 'basic'
    _description = 'Basic table with custom delimiters'
    _io_registry_format_aliases = ['ascii']

    header_class = BasicHeader
    data_class = BasicData


class NoHeaderHeader(BasicHeader):
    """
    Reader for table header without a header

    Set the start of header line number to `None`, which tells the basic
    reader there is no header line.
    """
    start_line = None


class NoHeaderData(BasicData):
    """
    Reader for table data without a header

    Data starts at first uncommented line since there is no header line.
    """
    start_line = 0


class NoHeader(Basic):
    """Character-delimited table with no header line.

    When reading, columns are autonamed using header.auto_format which defaults
    to "col%d".  Otherwise this reader the same as the :class:`Basic` class
    from which it is derived.  Example::

      # Table data
      1 2 "hello there"
      3 4 world

    """
    _format_name = 'no_header'
    _description = 'Basic table with no headers'
    header_class = NoHeaderHeader
    data_class = NoHeaderData


class CommentedHeaderHeader(BasicHeader):
    """
    Header class for which the column definition line starts with the
    comment character.  See the :class:`CommentedHeader` class  for an example.
    """

    def process_lines(self, lines):
        """
        Return only lines that start with the comment regexp.  For these
        lines strip out the matching characters.
        """
        re_comment = re.compile(self.comment)
        for line in lines:
            match = re_comment.match(line)
            if match:
                yield line[match.end():]

    def write(self, lines):
        lines.append(self.write_comment + self.splitter.join(self.colnames))


class CommentedHeader(Basic):
    """Character-delimited table with column names in a comment line.

    When reading, ``header_start`` can be used to specify the
    line index of column names, and it can be a negative index (for example -1
    for the last commented line).  The default delimiter is the <space>
    character.

    This matches the format produced by ``np.savetxt()``, with ``delimiter=','``,
    and ``header='<comma-delimited-column-names-list>'``.

    Example::

      # col1 col2 col3
      # Comment line
      1 2 3
      4 5 6

    """
    _format_name = 'commented_header'
    _description = 'Column names in a commented line'

    header_class = CommentedHeaderHeader
    data_class = NoHeaderData

    def read(self, table):
        """
        Read input data (file-like object, filename, list of strings, or
        single string) into a Table and return the result.
        """
        out = super().read(table)

        # Strip off the comment line set as the header line for
        # commented_header format (first by default).
        if 'comments' in out.meta:
            idx = self.header.start_line
            if idx < 0:
                idx = len(out.meta['comments']) + idx
            out.meta['comments'] = out.meta['comments'][:idx] + out.meta['comments'][idx + 1:]
            if not out.meta['comments']:
                del out.meta['comments']

        return out

    def write_header(self, lines, meta):
        """
        Write comment lines after, rather than before, the header.
        """
        self.header.write(lines)
        self.header.write_comments(lines, meta)


class TabHeaderSplitter(core.DefaultSplitter):
    """Split lines on tab and do not remove whitespace"""
    delimiter = '\t'
    process_line = None


class TabDataSplitter(TabHeaderSplitter):
    """
    Don't strip data value whitespace since that is significant in TSV tables
    """
    process_val = None
    skipinitialspace = False


class TabHeader(BasicHeader):
    """
    Reader for header of tables with tab separated header
    """
    splitter_class = TabHeaderSplitter


class TabData(BasicData):
    """
    Reader for data of tables with tab separated data
    """
    splitter_class = TabDataSplitter


class Tab(Basic):
    """Tab-separated table.

    Unlike the :class:`Basic` reader, whitespace is not stripped from the
    beginning and end of either lines or individual column values.

    Example::

      col1 <tab> col2 <tab> col3
      # Comment line
      1 <tab> 2 <tab> 5

    """
    _format_name = 'tab'
    _description = 'Basic table with tab-separated values'
    header_class = TabHeader
    data_class = TabData


class CsvSplitter(core.DefaultSplitter):
    """
    Split on comma for CSV (comma-separated-value) tables
    """
    delimiter = ','


class CsvHeader(BasicHeader):
    """
    Header that uses the :class:`astropy.io.ascii.basic.CsvSplitter`
    """
    splitter_class = CsvSplitter
    comment = None
    write_comment = None


class CsvData(BasicData):
    """
    Data that uses the :class:`astropy.io.ascii.basic.CsvSplitter`
    """
    splitter_class = CsvSplitter
    fill_values = [(core.masked, '')]
    comment = None
    write_comment = None


class Csv(Basic):
    """CSV (comma-separated-values) table.

    This file format may contain rows with fewer entries than the number of
    columns, a situation that occurs in output from some spreadsheet editors.
    The missing entries are marked as masked in the output table.

    Masked values (indicated by an empty '' field value when reading) are
    written out in the same way with an empty ('') field.  This is different
    from the typical default for `astropy.io.ascii` in which missing values are
    indicated by ``--``.

    Example::

      num,ra,dec,radius,mag
      1,32.23222,10.1211
      2,38.12321,-88.1321,2.2,17.0

    """
    _format_name = 'csv'
    _io_registry_format_aliases = ['csv']
    _io_registry_can_write = True
    _io_registry_suffix = '.csv'
    _description = 'Comma-separated-values'

    header_class = CsvHeader
    data_class = CsvData

    def inconsistent_handler(self, str_vals, ncols):
        """
        Adjust row if it is too short.

        If a data row is shorter than the header, add empty values to make it the
        right length.
        Note that this will *not* be called if the row already matches the header.

        Parameters
        ----------
        str_vals : list
            A list of value strings from the current row of the table.
        ncols : int
            The expected number of entries from the table header.

        Returns
        -------
        str_vals : list
            List of strings to be parsed into data entries in the output table.
        """
        if len(str_vals) < ncols:
            str_vals.extend((ncols - len(str_vals)) * [''])

        return str_vals


class RdbHeader(TabHeader):
    """
    Header for RDB tables
    """
    col_type_map = {'n': core.NumType,
                    's': core.StrType}

    def get_type_map_key(self, col):
        return col.raw_type[-1]

    def get_cols(self, lines):
        """
        Initialize the header Column objects from the table ``lines``.

        This is a specialized get_cols for the RDB type:
        Line 0: RDB col names
        Line 1: RDB col definitions
        Line 2+: RDB data rows


        Parameters
        ----------
        lines : list
            List of table lines

        Returns
        -------
        None

        """
        header_lines = self.process_lines(lines)   # this is a generator
        header_vals_list = [hl for _, hl in zip(range(2), self.splitter(header_lines))]
        if len(header_vals_list) != 2:
            raise ValueError('RDB header requires 2 lines')
        self.names, raw_types = header_vals_list

        if len(self.names) != len(raw_types):
            raise core.InconsistentTableError(
                'RDB header mismatch between number of column names and column types.')

        if any(not re.match(r'\d*(N|S)$', x, re.IGNORECASE) for x in raw_types):
            raise core.InconsistentTableError(
                f'RDB types definitions do not all match [num](N|S): {raw_types}')

        self._set_cols_from_names()
        for col, raw_type in zip(self.cols, raw_types):
            col.raw_type = raw_type
            col.type = self.get_col_type(col)

    def write(self, lines):
        lines.append(self.splitter.join(self.colnames))
        rdb_types = []
        for col in self.cols:
            # Check if dtype.kind is string or unicode.  See help(np.core.numerictypes)
            rdb_type = 'S' if col.info.dtype.kind in ('S', 'U') else 'N'
            rdb_types.append(rdb_type)

        lines.append(self.splitter.join(rdb_types))


class RdbData(TabData):
    """
    Data reader for RDB data. Starts reading at line 2.
    """
    start_line = 2


class Rdb(Tab):
    """Tab-separated file with an extra line after the column definition line that
    specifies either numeric (N) or string (S) data.

    See: https://www.drdobbs.com/rdb-a-unix-command-line-database/199101326

    Example::

      col1 <tab> col2 <tab> col3
      N <tab> S <tab> N
      1 <tab> 2 <tab> 5

    """
    _format_name = 'rdb'
    _io_registry_format_aliases = ['rdb']
    _io_registry_suffix = '.rdb'
    _description = 'Tab-separated with a type definition header line'

    header_class = RdbHeader
    data_class = RdbData
