# Licensed under a 3-clause BSD style license - see LICENSE.rst
""" An extensible ASCII table reader and writer.

core.py:
  Core base classes and functions for reading and writing tables.

:Copyright: Smithsonian Astrophysical Observatory (2010)
:Author: Tom Aldcroft (aldcroft@head.cfa.harvard.edu)
"""


import copy
import csv
import functools
import itertools
import operator
import os
import re
import warnings
import inspect
import fnmatch

from collections import OrderedDict
from contextlib import suppress
from io import StringIO

import numpy

from astropy.utils.exceptions import AstropyWarning

from astropy.table import Table
from astropy.utils.data import get_readable_fileobj
from . import connect
from .docs import READ_DOCSTRING, WRITE_DOCSTRING

# Global dictionary mapping format arg to the corresponding Reader class
FORMAT_CLASSES = {}

# Similar dictionary for fast readers
FAST_CLASSES = {}


def _check_multidim_table(table, max_ndim):
    """Check that ``table`` has only columns with ndim <= ``max_ndim``

    Currently ECSV is the only built-in format that supports output of arbitrary
    N-d columns, but HTML supports 2-d.
    """
    # No limit?
    if max_ndim is None:
        return

    # Check for N-d columns
    nd_names = [col.info.name for col in table.itercols() if len(col.shape) > max_ndim]
    if nd_names:
        raise ValueError(f'column(s) with dimension > {max_ndim} '
                         "cannot be be written with this format, try using 'ecsv' "
                         "(Enhanced CSV) format")


class CsvWriter:
    """
    Internal class to replace the csv writer ``writerow`` and ``writerows``
    functions so that in the case of ``delimiter=' '`` and
    ``quoting=csv.QUOTE_MINIMAL``, the output field value is quoted for empty
    fields (when value == '').

    This changes the API slightly in that the writerow() and writerows()
    methods return the output written string instead of the length of
    that string.

    Examples
    --------

    >>> from astropy.io.ascii.core import CsvWriter
    >>> writer = CsvWriter(delimiter=' ')
    >>> print(writer.writerow(['hello', '', 'world']))
    hello "" world
    """
    # Random 16-character string that gets injected instead of any
    # empty fields and is then replaced post-write with doubled-quotechar.
    # Created with:
    # ''.join(random.choice(string.printable[:90]) for _ in range(16))
    replace_sentinel = '2b=48Av%0-V3p>bX'

    def __init__(self, csvfile=None, **kwargs):
        self.csvfile = csvfile

        # Temporary StringIO for catching the real csv.writer() object output
        self.temp_out = StringIO()
        self.writer = csv.writer(self.temp_out, **kwargs)

        dialect = self.writer.dialect
        self.quotechar2 = dialect.quotechar * 2
        self.quote_empty = (dialect.quoting == csv.QUOTE_MINIMAL) and (dialect.delimiter == ' ')

    def writerow(self, values):
        """
        Similar to csv.writer.writerow but with the custom quoting behavior.
        Returns the written string instead of the length of that string.
        """
        has_empty = False

        # If QUOTE_MINIMAL and space-delimited then replace empty fields with
        # the sentinel value.
        if self.quote_empty:
            for i, value in enumerate(values):
                if value == '':
                    has_empty = True
                    values[i] = self.replace_sentinel

        return self._writerow(self.writer.writerow, values, has_empty)

    def writerows(self, values_list):
        """
        Similar to csv.writer.writerows but with the custom quoting behavior.
        Returns the written string instead of the length of that string.
        """
        has_empty = False

        # If QUOTE_MINIMAL and space-delimited then replace empty fields with
        # the sentinel value.
        if self.quote_empty:
            for values in values_list:
                for i, value in enumerate(values):
                    if value == '':
                        has_empty = True
                        values[i] = self.replace_sentinel

        return self._writerow(self.writer.writerows, values_list, has_empty)

    def _writerow(self, writerow_func, values, has_empty):
        """
        Call ``writerow_func`` (either writerow or writerows) with ``values``.
        If it has empty fields that have been replaced then change those
        sentinel strings back to quoted empty strings, e.g. ``""``.
        """
        # Clear the temporary StringIO buffer that self.writer writes into and
        # then call the real csv.writer().writerow or writerows with values.
        self.temp_out.seek(0)
        self.temp_out.truncate()
        writerow_func(values)

        row_string = self.temp_out.getvalue()

        if self.quote_empty and has_empty:
            row_string = re.sub(self.replace_sentinel, self.quotechar2, row_string)

        # self.csvfile is defined then write the output.  In practice the pure
        # Python writer calls with csvfile=None, while the fast writer calls with
        # a file-like object.
        if self.csvfile:
            self.csvfile.write(row_string)

        return row_string


class MaskedConstant(numpy.ma.core.MaskedConstant):
    """A trivial extension of numpy.ma.masked

    We want to be able to put the generic term ``masked`` into a dictionary.
    The constant ``numpy.ma.masked`` is not hashable (see
    https://github.com/numpy/numpy/issues/4660), so we need to extend it
    here with a hash value.

    See https://github.com/numpy/numpy/issues/11021 for rationale for
    __copy__ and __deepcopy__ methods.
    """

    def __hash__(self):
        '''All instances of this class shall have the same hash.'''
        # Any large number will do.
        return 1234567890

    def __copy__(self):
        """This is a singleton so just return self."""
        return self

    def __deepcopy__(self, memo):
        return self


masked = MaskedConstant()


class InconsistentTableError(ValueError):
    """
    Indicates that an input table is inconsistent in some way.

    The default behavior of ``BaseReader`` is to throw an instance of
    this class if a data row doesn't match the header.
    """


class OptionalTableImportError(ImportError):
    """
    Indicates that a dependency for table reading is not present.

    An instance of this class is raised whenever an optional reader
    with certain required dependencies cannot operate because of
    an ImportError.
    """


class ParameterError(NotImplementedError):
    """
    Indicates that a reader cannot handle a passed parameter.

    The C-based fast readers in ``io.ascii`` raise an instance of
    this error class upon encountering a parameter that the
    C engine cannot handle.
    """


class FastOptionsError(NotImplementedError):
    """
    Indicates that one of the specified options for fast
    reading is invalid.
    """


class NoType:
    """
    Superclass for ``StrType`` and ``NumType`` classes.

    This class is the default type of ``Column`` and provides a base
    class for other data types.
    """


class StrType(NoType):
    """
    Indicates that a column consists of text data.
    """


class NumType(NoType):
    """
    Indicates that a column consists of numerical data.
    """


class FloatType(NumType):
    """
    Describes floating-point data.
    """


class BoolType(NoType):
    """
    Describes boolean data.
    """


class IntType(NumType):
    """
    Describes integer data.
    """


class AllType(StrType, FloatType, IntType):
    """
    Subclass of all other data types.

    This type is returned by ``convert_numpy`` if the given numpy
    type does not match ``StrType``, ``FloatType``, or ``IntType``.
    """


class Column:
    """Table column.

    The key attributes of a Column object are:

    * **name** : column name
    * **type** : column type (NoType, StrType, NumType, FloatType, IntType)
    * **dtype** : numpy dtype (optional, overrides **type** if set)
    * **str_vals** : list of column values as strings
    * **fill_values** : dict of fill values
    * **shape** : list of element shape (default [] => scalar)
    * **data** : list of converted column values
    * **subtype** : actual datatype for columns serialized with JSON
    """

    def __init__(self, name):
        self.name = name
        self.type = NoType  # Generic type (Int, Float, Str etc)
        self.dtype = None  # Numpy dtype if available
        self.str_vals = []
        self.fill_values = {}
        self.shape = []
        self.subtype = None


class BaseInputter:
    """
    Get the lines from the table input and return a list of lines.

    """

    encoding = None
    """Encoding used to read the file"""

    def get_lines(self, table, newline=None):
        """
        Get the lines from the ``table`` input. The input table can be one of:

        * File name
        * String (newline separated) with all header and data lines (must have at least 2 lines)
        * File-like object with read() method
        * List of strings

        Parameters
        ----------
        table : str, file-like, list
            Can be either a file name, string (newline separated) with all header and data
            lines (must have at least 2 lines), a file-like object with a
            ``read()`` method, or a list of strings.
        newline :
            Line separator. If `None` use OS default from ``splitlines()``.

        Returns
        -------
        lines : list
            List of lines
        """
        try:
            if (hasattr(table, 'read')
                    or ('\n' not in table + '' and '\r' not in table + '')):
                with get_readable_fileobj(table,
                                          encoding=self.encoding) as fileobj:
                    table = fileobj.read()
            if newline is None:
                lines = table.splitlines()
            else:
                lines = table.split(newline)
        except TypeError:
            try:
                # See if table supports indexing, slicing, and iteration
                table[0]
                table[0:1]
                iter(table)
                if len(table) > 1:
                    lines = table
                else:
                    # treat single entry as if string had been passed directly
                    if newline is None:
                        lines = table[0].splitlines()
                    else:
                        lines = table[0].split(newline)

            except TypeError:
                raise TypeError(
                    'Input "table" must be a string (filename or data) or an iterable')

        return self.process_lines(lines)

    def process_lines(self, lines):
        """Process lines for subsequent use.  In the default case do nothing.
        This routine is not generally intended for removing comment lines or
        stripping whitespace.  These are done (if needed) in the header and
        data line processing.

        Override this method if something more has to be done to convert raw
        input lines to the table rows.  For example the
        ContinuationLinesInputter derived class accounts for continuation
        characters if a row is split into lines."""
        return lines


class BaseSplitter:
    """
    Base splitter that uses python's split method to do the work.

    This does not handle quoted values.  A key feature is the formulation of
    __call__ as a generator that returns a list of the split line values at
    each iteration.

    There are two methods that are intended to be overridden, first
    ``process_line()`` to do pre-processing on each input line before splitting
    and ``process_val()`` to do post-processing on each split string value.  By
    default these apply the string ``strip()`` function.  These can be set to
    another function via the instance attribute or be disabled entirely, for
    example::

      reader.header.splitter.process_val = lambda x: x.lstrip()
      reader.data.splitter.process_val = None

    """

    delimiter = None
    """ one-character string used to separate fields """

    def process_line(self, line):
        """Remove whitespace at the beginning or end of line.  This is especially useful for
        whitespace-delimited files to prevent spurious columns at the beginning or end."""
        return line.strip()

    def process_val(self, val):
        """Remove whitespace at the beginning or end of value."""
        return val.strip()

    def __call__(self, lines):
        if self.process_line:
            lines = (self.process_line(x) for x in lines)
        for line in lines:
            vals = line.split(self.delimiter)
            if self.process_val:
                yield [self.process_val(x) for x in vals]
            else:
                yield vals

    def join(self, vals):
        if self.delimiter is None:
            delimiter = ' '
        else:
            delimiter = self.delimiter
        return delimiter.join(str(x) for x in vals)


class DefaultSplitter(BaseSplitter):
    """Default class to split strings into columns using python csv.  The class
    attributes are taken from the csv Dialect class.

    Typical usage::

      # lines = ..
      splitter = ascii.DefaultSplitter()
      for col_vals in splitter(lines):
          for col_val in col_vals:
               ...

    """
    delimiter = ' '
    """ one-character string used to separate fields. """
    quotechar = '"'
    """ control how instances of *quotechar* in a field are quoted """
    doublequote = True
    """ character to remove special meaning from following character """
    escapechar = None
    """ one-character stringto quote fields containing special characters """
    quoting = csv.QUOTE_MINIMAL
    """ control when quotes are recognized by the reader """
    skipinitialspace = True
    """ ignore whitespace immediately following the delimiter """
    csv_writer = None
    csv_writer_out = StringIO()

    def process_line(self, line):
        """Remove whitespace at the beginning or end of line.  This is especially useful for
        whitespace-delimited files to prevent spurious columns at the beginning or end.
        If splitting on whitespace then replace unquoted tabs with space first"""
        if self.delimiter == r'\s':
            line = _replace_tab_with_space(line, self.escapechar, self.quotechar)
        return line.strip() + '\n'

    def process_val(self, val):
        """Remove whitespace at the beginning or end of value."""
        return val.strip(' \t')

    def __call__(self, lines):
        """Return an iterator over the table ``lines``, where each iterator output
        is a list of the split line values.

        Parameters
        ----------
        lines : list
            List of table lines

        Yields
        ------
        line : list of str
            Each line's split values.

        """
        if self.process_line:
            lines = [self.process_line(x) for x in lines]

        delimiter = ' ' if self.delimiter == r'\s' else self.delimiter

        csv_reader = csv.reader(lines,
                                delimiter=delimiter,
                                doublequote=self.doublequote,
                                escapechar=self.escapechar,
                                quotechar=self.quotechar,
                                quoting=self.quoting,
                                skipinitialspace=self.skipinitialspace
                                )
        for vals in csv_reader:
            if self.process_val:
                yield [self.process_val(x) for x in vals]
            else:
                yield vals

    def join(self, vals):

        delimiter = ' ' if self.delimiter is None else str(self.delimiter)

        if self.csv_writer is None:
            self.csv_writer = CsvWriter(delimiter=delimiter,
                                        doublequote=self.doublequote,
                                        escapechar=self.escapechar,
                                        quotechar=self.quotechar,
                                        quoting=self.quoting)
        if self.process_val:
            vals = [self.process_val(x) for x in vals]
        out = self.csv_writer.writerow(vals).rstrip('\r\n')

        return out


def _replace_tab_with_space(line, escapechar, quotechar):
    """Replace tabs with spaces in given string, preserving quoted substrings

    Parameters
    ----------
    line : str
        String containing tabs to be replaced with spaces.
    escapechar : str
        Character in ``line`` used to escape special characters.
    quotechar : str
        Character in ``line`` indicating the start/end of a substring.

    Returns
    -------
    line : str
        A copy of ``line`` with tabs replaced by spaces, preserving quoted substrings.
    """
    newline = []
    in_quote = False
    lastchar = 'NONE'
    for char in line:
        if char == quotechar and lastchar != escapechar:
            in_quote = not in_quote
        if char == '\t' and not in_quote:
            char = ' '
        lastchar = char
        newline.append(char)
    return ''.join(newline)


def _get_line_index(line_or_func, lines):
    """Return the appropriate line index, depending on ``line_or_func`` which
    can be either a function, a positive or negative int, or None.
    """

    if hasattr(line_or_func, '__call__'):
        return line_or_func(lines)
    elif line_or_func:
        if line_or_func >= 0:
            return line_or_func
        else:
            n_lines = sum(1 for line in lines)
            return n_lines + line_or_func
    else:
        return line_or_func


class BaseHeader:
    """
    Base table header reader
    """
    auto_format = 'col{}'
    """ format string for auto-generating column names """
    start_line = None
    """ None, int, or a function of ``lines`` that returns None or int """
    comment = None
    """ regular expression for comment lines """
    splitter_class = DefaultSplitter
    """ Splitter class for splitting data lines into columns """
    names = None
    """ list of names corresponding to each data column """
    write_comment = False
    write_spacer_lines = ['ASCII_TABLE_WRITE_SPACER_LINE']

    def __init__(self):
        self.splitter = self.splitter_class()

    def _set_cols_from_names(self):
        self.cols = [Column(name=x) for x in self.names]

    def update_meta(self, lines, meta):
        """
        Extract any table-level metadata, e.g. keywords, comments, column metadata, from
        the table ``lines`` and update the OrderedDict ``meta`` in place.  This base
        method extracts comment lines and stores them in ``meta`` for output.
        """
        if self.comment:
            re_comment = re.compile(self.comment)
            comment_lines = [x for x in lines if re_comment.match(x)]
        else:
            comment_lines = []
        comment_lines = [re.sub('^' + self.comment, '', x).strip()
                         for x in comment_lines]
        if comment_lines:
            meta.setdefault('table', {})['comments'] = comment_lines

    def get_cols(self, lines):
        """Initialize the header Column objects from the table ``lines``.

        Based on the previously set Header attributes find or create the column names.
        Sets ``self.cols`` with the list of Columns.

        Parameters
        ----------
        lines : list
            List of table lines

        """

        start_line = _get_line_index(self.start_line, self.process_lines(lines))
        if start_line is None:
            # No header line so auto-generate names from n_data_cols
            # Get the data values from the first line of table data to determine n_data_cols
            try:
                first_data_vals = next(self.data.get_str_vals())
            except StopIteration:
                raise InconsistentTableError('No data lines found so cannot autogenerate '
                                             'column names')
            n_data_cols = len(first_data_vals)
            self.names = [self.auto_format.format(i)
                          for i in range(1, n_data_cols + 1)]

        else:
            for i, line in enumerate(self.process_lines(lines)):
                if i == start_line:
                    break
            else:  # No header line matching
                raise ValueError('No header line found in table')

            self.names = next(self.splitter([line]))

        self._set_cols_from_names()

    def process_lines(self, lines):
        """Generator to yield non-blank and non-comment lines"""
        re_comment = re.compile(self.comment) if self.comment else None
        # Yield non-comment lines
        for line in lines:
            if line.strip() and (not self.comment or not re_comment.match(line)):
                yield line

    def write_comments(self, lines, meta):
        if self.write_comment not in (False, None):
            for comment in meta.get('comments', []):
                lines.append(self.write_comment + comment)

    def write(self, lines):
        if self.start_line is not None:
            for i, spacer_line in zip(range(self.start_line),
                                      itertools.cycle(self.write_spacer_lines)):
                lines.append(spacer_line)
            lines.append(self.splitter.join([x.info.name for x in self.cols]))

    @property
    def colnames(self):
        """Return the column names of the table"""
        return tuple(col.name if isinstance(col, Column) else col.info.name
                     for col in self.cols)

    def remove_columns(self, names):
        """
        Remove several columns from the table.

        Parameters
        ----------
        names : list
            A list containing the names of the columns to remove
        """
        colnames = self.colnames
        for name in names:
            if name not in colnames:
                raise KeyError(f"Column {name} does not exist")

        self.cols = [col for col in self.cols if col.name not in names]

    def rename_column(self, name, new_name):
        """
        Rename a column.

        Parameters
        ----------
        name : str
            The current name of the column.
        new_name : str
            The new name for the column
        """
        try:
            idx = self.colnames.index(name)
        except ValueError:
            raise KeyError(f"Column {name} does not exist")

        col = self.cols[idx]

        # For writing self.cols can contain cols that are not Column.  Raise
        # exception in that case.
        if isinstance(col, Column):
            col.name = new_name
        else:
            raise TypeError(f'got column type {type(col)} instead of required '
                            f'{Column}')

    def get_type_map_key(self, col):
        return col.raw_type

    def get_col_type(self, col):
        try:
            type_map_key = self.get_type_map_key(col)
            return self.col_type_map[type_map_key.lower()]
        except KeyError:
            raise ValueError('Unknown data type ""{}"" for column "{}"'.format(
                col.raw_type, col.name))

    def check_column_names(self, names, strict_names, guessing):
        """
        Check column names.

        This must be done before applying the names transformation
        so that guessing will fail appropriately if ``names`` is supplied.
        For instance if the basic reader is given a table with no column header
        row.

        Parameters
        ----------
        names : list
            User-supplied list of column names
        strict_names : bool
            Whether to impose extra requirements on names
        guessing : bool
            True if this method is being called while guessing the table format
        """
        if strict_names:
            # Impose strict requirements on column names (normally used in guessing)
            bads = [" ", ",", "|", "\t", "'", '"']
            for name in self.colnames:
                if (_is_number(name) or len(name) == 0
                        or name[0] in bads or name[-1] in bads):
                    raise InconsistentTableError(
                        f'Column name {name!r} does not meet strict name requirements')
        # When guessing require at least two columns, except for ECSV which can
        # reliably be guessed from the header requirements.
        if guessing and len(self.colnames) <= 1 and self.__class__.__name__ != 'EcsvHeader':
            raise ValueError('Table format guessing requires at least two columns, got {}'
                             .format(list(self.colnames)))

        if names is not None and len(names) != len(self.colnames):
            raise InconsistentTableError(
                'Length of names argument ({}) does not match number'
                ' of table columns ({})'.format(len(names), len(self.colnames)))


class BaseData:
    """
    Base table data reader.
    """
    start_line = None
    """ None, int, or a function of ``lines`` that returns None or int """
    end_line = None
    """ None, int, or a function of ``lines`` that returns None or int """
    comment = None
    """ Regular expression for comment lines """
    splitter_class = DefaultSplitter
    """ Splitter class for splitting data lines into columns """
    write_spacer_lines = ['ASCII_TABLE_WRITE_SPACER_LINE']
    fill_include_names = None
    fill_exclude_names = None
    fill_values = [(masked, '')]
    formats = {}

    def __init__(self):
        # Need to make sure fill_values list is instance attribute, not class attribute.
        # On read, this will be overwritten by the default in the ui.read (thus, in
        # the current implementation there can be no different default for different
        # Readers). On write, ui.py does not specify a default, so this line here matters.
        self.fill_values = copy.copy(self.fill_values)
        self.formats = copy.copy(self.formats)
        self.splitter = self.splitter_class()

    def process_lines(self, lines):
        """
        READ: Strip out comment lines and blank lines from list of ``lines``

        Parameters
        ----------
        lines : list
            All lines in table

        Returns
        -------
        lines : list
            List of lines

        """
        nonblank_lines = (x for x in lines if x.strip())
        if self.comment:
            re_comment = re.compile(self.comment)
            return [x for x in nonblank_lines if not re_comment.match(x)]
        else:
            return [x for x in nonblank_lines]

    def get_data_lines(self, lines):
        """READ: Set ``data_lines`` attribute to lines slice comprising table data values.
        """
        data_lines = self.process_lines(lines)
        start_line = _get_line_index(self.start_line, data_lines)
        end_line = _get_line_index(self.end_line, data_lines)

        if start_line is not None or end_line is not None:
            self.data_lines = data_lines[slice(start_line, end_line)]
        else:  # Don't copy entire data lines unless necessary
            self.data_lines = data_lines

    def get_str_vals(self):
        """Return a generator that returns a list of column values (as strings)
        for each data line."""
        return self.splitter(self.data_lines)

    def masks(self, cols):
        """READ: Set fill value for each column and then apply that fill value

        In the first step it is evaluated with value from ``fill_values`` applies to
        which column using ``fill_include_names`` and ``fill_exclude_names``.
        In the second step all replacements are done for the appropriate columns.
        """
        if self.fill_values:
            self._set_fill_values(cols)
            self._set_masks(cols)

    def _set_fill_values(self, cols):
        """READ, WRITE: Set fill values of individual cols based on fill_values of BaseData

        fill values has the following form:
        <fill_spec> = (<bad_value>, <fill_value>, <optional col_name>...)
        fill_values = <fill_spec> or list of <fill_spec>'s

        """
        if self.fill_values:
            # when we write tables the columns may be astropy.table.Columns
            # which don't carry a fill_values by default
            for col in cols:
                if not hasattr(col, 'fill_values'):
                    col.fill_values = {}

            # if input is only one <fill_spec>, then make it a list
            with suppress(TypeError):
                self.fill_values[0] + ''
                self.fill_values = [self.fill_values]

            # Step 1: Set the default list of columns which are affected by
            # fill_values
            colnames = set(self.header.colnames)
            if self.fill_include_names is not None:
                colnames.intersection_update(self.fill_include_names)
            if self.fill_exclude_names is not None:
                colnames.difference_update(self.fill_exclude_names)

            # Step 2a: Find out which columns are affected by this tuple
            # iterate over reversed order, so last condition is set first and
            # overwritten by earlier conditions
            for replacement in reversed(self.fill_values):
                if len(replacement) < 2:
                    raise ValueError("Format of fill_values must be "
                                     "(<bad>, <fill>, <optional col1>, ...)")
                elif len(replacement) == 2:
                    affect_cols = colnames
                else:
                    affect_cols = replacement[2:]

                for i, key in ((i, x) for i, x in enumerate(self.header.colnames)
                               if x in affect_cols):
                    cols[i].fill_values[replacement[0]] = str(replacement[1])

    def _set_masks(self, cols):
        """READ: Replace string values in col.str_vals and set masks"""
        if self.fill_values:
            for col in (col for col in cols if col.fill_values):
                col.mask = numpy.zeros(len(col.str_vals), dtype=bool)
                for i, str_val in ((i, x) for i, x in enumerate(col.str_vals)
                                   if x in col.fill_values):
                    col.str_vals[i] = col.fill_values[str_val]
                    col.mask[i] = True

    def _replace_vals(self, cols):
        """WRITE: replace string values in col.str_vals"""
        if self.fill_values:
            for col in (col for col in cols if col.fill_values):
                for i, str_val in ((i, x) for i, x in enumerate(col.str_vals)
                                   if x in col.fill_values):
                    col.str_vals[i] = col.fill_values[str_val]
                if masked in col.fill_values and hasattr(col, 'mask'):
                    mask_val = col.fill_values[masked]
                    for i in col.mask.nonzero()[0]:
                        col.str_vals[i] = mask_val

    def str_vals(self):
        """WRITE: convert all values in table to a list of lists of strings

        This sets the fill values and possibly column formats from the input
        formats={} keyword, then ends up calling table.pprint._pformat_col_iter()
        by a circuitous path. That function does the real work of formatting.
        Finally replace anything matching the fill_values.

        Returns
        -------
        values : list of list of str
        """
        self._set_fill_values(self.cols)
        self._set_col_formats()
        for col in self.cols:
            col.str_vals = list(col.info.iter_str_vals())
        self._replace_vals(self.cols)
        return [col.str_vals for col in self.cols]

    def write(self, lines):
        """Write ``self.cols`` in place to ``lines``.

        Parameters
        ----------
        lines : list
            List for collecting output of writing self.cols.
        """
        if hasattr(self.start_line, '__call__'):
            raise TypeError('Start_line attribute cannot be callable for write()')
        else:
            data_start_line = self.start_line or 0

        while len(lines) < data_start_line:
            lines.append(itertools.cycle(self.write_spacer_lines))

        col_str_iters = self.str_vals()
        for vals in zip(*col_str_iters):
            lines.append(self.splitter.join(vals))

    def _set_col_formats(self):
        """WRITE: set column formats."""
        for col in self.cols:
            if col.info.name in self.formats:
                col.info.format = self.formats[col.info.name]


def convert_numpy(numpy_type):
    """Return a tuple containing a function which converts a list into a numpy
    array and the type produced by the converter function.

    Parameters
    ----------
    numpy_type : numpy data-type
        The numpy type required of an array returned by ``converter``. Must be a
        valid `numpy type <https://numpy.org/doc/stable/user/basics.types.html>`_
        (e.g., numpy.uint, numpy.int8, numpy.int64, numpy.float64) or a python
        type covered by a numpy type (e.g., int, float, str, bool).

    Returns
    -------
    converter : callable
        ``converter`` is a function which accepts a list and converts it to a
        numpy array of type ``numpy_type``.
    converter_type : type
        ``converter_type`` tracks the generic data type produced by the
        converter function.

    Raises
    ------
    ValueError
        Raised by ``converter`` if the list elements could not be converted to
        the required type.
    """

    # Infer converter type from an instance of numpy_type.
    type_name = numpy.array([], dtype=numpy_type).dtype.name
    if 'int' in type_name:
        converter_type = IntType
    elif 'float' in type_name:
        converter_type = FloatType
    elif 'bool' in type_name:
        converter_type = BoolType
    elif 'str' in type_name:
        converter_type = StrType
    else:
        converter_type = AllType

    def bool_converter(vals):
        """
        Convert values "False" and "True" to bools.  Raise an exception
        for any other string values.
        """
        if len(vals) == 0:
            return numpy.array([], dtype=bool)

        # Try a smaller subset first for a long array
        if len(vals) > 10000:
            svals = numpy.asarray(vals[:1000])
            if not numpy.all((svals == 'False')
                             | (svals == 'True')
                             | (svals == '0')
                             | (svals == '1')):
                raise ValueError('bool input strings must be False, True, 0, 1, or ""')
        vals = numpy.asarray(vals)

        trues = (vals == 'True') | (vals == '1')
        falses = (vals == 'False') | (vals == '0')
        if not numpy.all(trues | falses):
            raise ValueError('bool input strings must be only False, True, 0, 1, or ""')

        return trues

    def generic_converter(vals):
        return numpy.array(vals, numpy_type)

    converter = bool_converter if converter_type is BoolType else generic_converter

    return converter, converter_type


class BaseOutputter:
    """Output table as a dict of column objects keyed on column name.  The
    table data are stored as plain python lists within the column objects.
    """
    # User-defined converters which gets set in ascii.ui if a `converter` kwarg
    # is supplied.
    converters = {}

    # Derived classes must define default_converters and __call__

    @staticmethod
    def _validate_and_copy(col, converters):
        """Validate the format for the type converters and then copy those
        which are valid converters for this column (i.e. converter type is
        a subclass of col.type)"""
        # Allow specifying a single converter instead of a list of converters.
        # The input `converters` must be a ``type`` value that can init np.dtype.
        try:
            # Don't allow list-like things that dtype accepts
            assert type(converters) is type
            converters = [numpy.dtype(converters)]
        except (AssertionError, TypeError):
            pass

        converters_out = []
        try:
            for converter in converters:
                try:
                    converter_func, converter_type = converter
                except TypeError as err:
                    if str(err).startswith('cannot unpack'):
                        converter_func, converter_type = convert_numpy(converter)
                    else:
                        raise
                if not issubclass(converter_type, NoType):
                    raise ValueError('converter_type must be a subclass of NoType')
                if issubclass(converter_type, col.type):
                    converters_out.append((converter_func, converter_type))

        except (ValueError, TypeError) as err:
            raise ValueError('Error: invalid format for converters, see '
                             f'documentation\n{converters}: {err}')
        return converters_out

    def _convert_vals(self, cols):
        for col in cols:
            for key, converters in self.converters.items():
                if fnmatch.fnmatch(col.name, key):
                    break
            else:
                if col.dtype is not None:
                    converters = [convert_numpy(col.dtype)]
                else:
                    converters = self.default_converters

            col.converters = self._validate_and_copy(col, converters)

            # Catch the last error in order to provide additional information
            # in case all attempts at column conversion fail.  The initial
            # value of of last_error will apply if no converters are defined
            # and the first col.converters[0] access raises IndexError.
            last_err = 'no converters defined'

            while not hasattr(col, 'data'):
                # Try converters, popping the unsuccessful ones from the list.
                # If there are no converters left here then fail.
                if not col.converters:
                    raise ValueError(f'Column {col.name} failed to convert: {last_err}')

                converter_func, converter_type = col.converters[0]
                if not issubclass(converter_type, col.type):
                    raise TypeError('converter type does not match column type')

                try:
                    col.data = converter_func(col.str_vals)
                    col.type = converter_type
                except (TypeError, ValueError) as err:
                    col.converters.pop(0)
                    last_err = err
                except OverflowError as err:
                    # Overflow during conversion (most likely an int that
                    # doesn't fit in native C long). Put string at the top of
                    # the converters list for the next while iteration.
                    warnings.warn(
                        "OverflowError converting to {} in column {}, reverting to String."
                        .format(converter_type.__name__, col.name), AstropyWarning)
                    col.converters.insert(0, convert_numpy(numpy.str))
                    last_err = err


def _deduplicate_names(names):
    """Ensure there are no duplicates in ``names``

    This is done by iteratively adding ``_<N>`` to the name for increasing N
    until the name is unique.
    """
    new_names = []
    existing_names = set()

    for name in names:
        base_name = name + '_'
        i = 1
        while name in existing_names:
            # Iterate until a unique name is found
            name = base_name + str(i)
            i += 1
        new_names.append(name)
        existing_names.add(name)

    return new_names


class TableOutputter(BaseOutputter):
    """
    Output the table as an astropy.table.Table object.
    """

    default_converters = [convert_numpy(int),
                          convert_numpy(float),
                          convert_numpy(str)]

    def __call__(self, cols, meta):
        # Sets col.data to numpy array and col.type to io.ascii Type class (e.g.
        # FloatType) for each col.
        self._convert_vals(cols)

        t_cols = [numpy.ma.MaskedArray(x.data, mask=x.mask)
                  if hasattr(x, 'mask') and numpy.any(x.mask)
                  else x.data for x in cols]
        out = Table(t_cols, names=[x.name for x in cols], meta=meta['table'])

        for col, out_col in zip(cols, out.columns.values()):
            for attr in ('format', 'unit', 'description'):
                if hasattr(col, attr):
                    setattr(out_col, attr, getattr(col, attr))
            if hasattr(col, 'meta'):
                out_col.meta.update(col.meta)

        return out


class MetaBaseReader(type):
    def __init__(cls, name, bases, dct):
        super().__init__(name, bases, dct)

        format = dct.get('_format_name')
        if format is None:
            return

        fast = dct.get('_fast')
        if fast is not None:
            FAST_CLASSES[format] = cls

        FORMAT_CLASSES[format] = cls

        io_formats = ['ascii.' + format] + dct.get('_io_registry_format_aliases', [])

        if dct.get('_io_registry_suffix'):
            func = functools.partial(connect.io_identify, dct['_io_registry_suffix'])
            connect.io_registry.register_identifier(io_formats[0], Table, func)

        for io_format in io_formats:
            func = functools.partial(connect.io_read, io_format)
            header = f"ASCII reader '{io_format}' details\n"
            func.__doc__ = (inspect.cleandoc(READ_DOCSTRING).strip() + '\n\n'
                            + header + re.sub('.', '=', header) + '\n')
            func.__doc__ += inspect.cleandoc(cls.__doc__).strip()
            connect.io_registry.register_reader(io_format, Table, func)

            if dct.get('_io_registry_can_write', True):
                func = functools.partial(connect.io_write, io_format)
                header = f"ASCII writer '{io_format}' details\n"
                func.__doc__ = (inspect.cleandoc(WRITE_DOCSTRING).strip() + '\n\n'
                                + header + re.sub('.', '=', header) + '\n')
                func.__doc__ += inspect.cleandoc(cls.__doc__).strip()
                connect.io_registry.register_writer(io_format, Table, func)


def _is_number(x):
    with suppress(ValueError):
        x = float(x)
        return True
    return False


def _apply_include_exclude_names(table, names, include_names, exclude_names):
    """
    Apply names, include_names and exclude_names to a table or BaseHeader.

    For the latter this relies on BaseHeader implementing ``colnames``,
    ``rename_column``, and ``remove_columns``.

    Parameters
    ----------
    table : `~astropy.table.Table`, `~astropy.io.ascii.BaseHeader`
        Input table or BaseHeader subclass instance
    names : list
        List of names to override those in table (set to None to use existing names)
    include_names : list
        List of names to include in output
    exclude_names : list
        List of names to exclude from output (applied after ``include_names``)

    """
    def rename_columns(table, names):
        # Rename table column names to those passed by user
        # Temporarily rename with names that are not in `names` or `table.colnames`.
        # This ensures that rename succeeds regardless of existing names.
        xxxs = 'x' * max(len(name) for name in list(names) + list(table.colnames))
        for ii, colname in enumerate(table.colnames):
            table.rename_column(colname, xxxs + str(ii))

        for ii, name in enumerate(names):
            table.rename_column(xxxs + str(ii), name)

    if names is not None:
        rename_columns(table, names)
    else:
        colnames_uniq = _deduplicate_names(table.colnames)
        if colnames_uniq != list(table.colnames):
            rename_columns(table, colnames_uniq)

    names_set = set(table.colnames)

    if include_names is not None:
        names_set.intersection_update(include_names)
    if exclude_names is not None:
        names_set.difference_update(exclude_names)
    if names_set != set(table.colnames):
        remove_names = set(table.colnames) - names_set
        table.remove_columns(remove_names)


class BaseReader(metaclass=MetaBaseReader):
    """Class providing methods to read and write an ASCII table using the specified
    header, data, inputter, and outputter instances.

    Typical usage is to instantiate a Reader() object and customize the
    ``header``, ``data``, ``inputter``, and ``outputter`` attributes.  Each
    of these is an object of the corresponding class.

    There is one method ``inconsistent_handler`` that can be used to customize the
    behavior of ``read()`` in the event that a data row doesn't match the header.
    The default behavior is to raise an InconsistentTableError.

    """

    names = None
    include_names = None
    exclude_names = None
    strict_names = False
    guessing = False
    encoding = None

    header_class = BaseHeader
    data_class = BaseData
    inputter_class = BaseInputter
    outputter_class = TableOutputter

    # Max column dimension that writer supports for this format. Exceptions
    # include ECSV (no limit) and HTML (max_ndim=2).
    max_ndim = 1

    def __init__(self):
        self.header = self.header_class()
        self.data = self.data_class()
        self.inputter = self.inputter_class()
        self.outputter = self.outputter_class()
        # Data and Header instances benefit from a little cross-coupling.  Header may need to
        # know about number of data columns for auto-column name generation and Data may
        # need to know about header (e.g. for fixed-width tables where widths are spec'd in header.
        self.data.header = self.header
        self.header.data = self.data

        # Metadata, consisting of table-level meta and column-level meta.  The latter
        # could include information about column type, description, formatting, etc,
        # depending on the table meta format.
        self.meta = OrderedDict(table=OrderedDict(),
                                cols=OrderedDict())

    def _check_multidim_table(self, table):
        """Check that the dimensions of columns in ``table`` are acceptable.

        The reader class attribute ``max_ndim`` defines the maximum dimension of
        columns that can be written using this format. The base value is ``1``,
        corresponding to normal scalar columns with just a length.

        Parameters
        ----------
        table : `~astropy.table.Table`
            Input table.

        Raises
        ------
        ValueError
            If any column exceeds the number of allowed dimensions
        """
        _check_multidim_table(table, self.max_ndim)

    def read(self, table):
        """Read the ``table`` and return the results in a format determined by
        the ``outputter`` attribute.

        The ``table`` parameter is any string or object that can be processed
        by the instance ``inputter``.  For the base Inputter class ``table`` can be
        one of:

        * File name
        * File-like object
        * String (newline separated) with all header and data lines (must have at least 2 lines)
        * List of strings

        Parameters
        ----------
        table : str, file-like, list
            Input table.

        Returns
        -------
        table : `~astropy.table.Table`
            Output table

        """
        # If ``table`` is a file then store the name in the ``data``
        # attribute. The ``table`` is a "file" if it is a string
        # without the new line specific to the OS.
        with suppress(TypeError):
            # Strings only
            if os.linesep not in table + '':
                self.data.table_name = os.path.basename(table)

        # If one of the newline chars is set as field delimiter, only
        # accept the other one as line splitter
        if self.header.splitter.delimiter == '\n':
            newline = '\r'
        elif self.header.splitter.delimiter == '\r':
            newline = '\n'
        else:
            newline = None

        # Get a list of the lines (rows) in the table
        self.lines = self.inputter.get_lines(table, newline=newline)

        # Set self.data.data_lines to a slice of lines contain the data rows
        self.data.get_data_lines(self.lines)

        # Extract table meta values (e.g. keywords, comments, etc).  Updates self.meta.
        self.header.update_meta(self.lines, self.meta)

        # Get the table column definitions
        self.header.get_cols(self.lines)

        # Make sure columns are valid
        self.header.check_column_names(self.names, self.strict_names, self.guessing)

        self.cols = cols = self.header.cols
        self.data.splitter.cols = cols
        n_cols = len(cols)

        for i, str_vals in enumerate(self.data.get_str_vals()):
            if len(str_vals) != n_cols:
                str_vals = self.inconsistent_handler(str_vals, n_cols)

                # if str_vals is None, we skip this row
                if str_vals is None:
                    continue

                # otherwise, we raise an error only if it is still inconsistent
                if len(str_vals) != n_cols:
                    errmsg = ('Number of header columns ({}) inconsistent with'
                              ' data columns ({}) at data line {}\n'
                              'Header values: {}\n'
                              'Data values: {}'.format(
                                  n_cols, len(str_vals), i,
                                  [x.name for x in cols], str_vals))

                    raise InconsistentTableError(errmsg)

            for j, col in enumerate(cols):
                col.str_vals.append(str_vals[j])

        self.data.masks(cols)
        if hasattr(self.header, 'table_meta'):
            self.meta['table'].update(self.header.table_meta)

        _apply_include_exclude_names(self.header, self.names,
                                     self.include_names, self.exclude_names)

        table = self.outputter(self.header.cols, self.meta)
        self.cols = self.header.cols

        return table

    def inconsistent_handler(self, str_vals, ncols):
        """
        Adjust or skip data entries if a row is inconsistent with the header.

        The default implementation does no adjustment, and hence will always trigger
        an exception in read() any time the number of data entries does not match
        the header.

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
            List of strings to be parsed into data entries in the output table. If
            the length of this list does not match ``ncols``, an exception will be
            raised in read().  Can also be None, in which case the row will be
            skipped.
        """
        # an empty list will always trigger an InconsistentTableError in read()
        return str_vals

    @property
    def comment_lines(self):
        """Return lines in the table that match header.comment regexp"""
        if not hasattr(self, 'lines'):
            raise ValueError('Table must be read prior to accessing the header comment lines')
        if self.header.comment:
            re_comment = re.compile(self.header.comment)
            comment_lines = [x for x in self.lines if re_comment.match(x)]
        else:
            comment_lines = []
        return comment_lines

    def update_table_data(self, table):
        """
        Update table columns in place if needed.

        This is a hook to allow updating the table columns after name
        filtering but before setting up to write the data.  This is currently
        only used by ECSV and is otherwise just a pass-through.

        Parameters
        ----------
        table : `astropy.table.Table`
            Input table for writing

        Returns
        -------
        table : `astropy.table.Table`
            Output table for writing
        """
        return table

    def write_header(self, lines, meta):
        self.header.write_comments(lines, meta)
        self.header.write(lines)

    def write(self, table):
        """
        Write ``table`` as list of strings.

        Parameters
        ----------
        table : `~astropy.table.Table`
            Input table data.

        Returns
        -------
        lines : list
            List of strings corresponding to ASCII table

        """

        # Check column names before altering
        self.header.cols = list(table.columns.values())
        self.header.check_column_names(self.names, self.strict_names, False)

        # In-place update of columns in input ``table`` to reflect column
        # filtering.  Note that ``table`` is guaranteed to be a copy of the
        # original user-supplied table.
        _apply_include_exclude_names(table, self.names, self.include_names, self.exclude_names)

        # This is a hook to allow updating the table columns after name
        # filtering but before setting up to write the data.  This is currently
        # only used by ECSV and is otherwise just a pass-through.
        table = self.update_table_data(table)

        # Check that table column dimensions are supported by this format class.
        # Most formats support only 1-d columns, but some like ECSV support N-d.
        self._check_multidim_table(table)

        # Now use altered columns
        new_cols = list(table.columns.values())
        # link information about the columns to the writer object (i.e. self)
        self.header.cols = new_cols
        self.data.cols = new_cols
        self.header.table_meta = table.meta

        # Write header and data to lines list
        lines = []
        self.write_header(lines, table.meta)
        self.data.write(lines)

        return lines


class ContinuationLinesInputter(BaseInputter):
    """Inputter where lines ending in ``continuation_char`` are joined
    with the subsequent line.  Example::

      col1 col2 col3
      1 \
      2 3
      4 5 \
      6
    """

    continuation_char = '\\'
    replace_char = ' '
    # If no_continue is not None then lines matching this regex are not subject
    # to line continuation.  The initial use case here is Daophot.  In this
    # case the continuation character is just replaced with replace_char.
    no_continue = None

    def process_lines(self, lines):
        re_no_continue = re.compile(self.no_continue) if self.no_continue else None

        parts = []
        outlines = []
        for line in lines:
            if re_no_continue and re_no_continue.match(line):
                line = line.replace(self.continuation_char, self.replace_char)
            if line.endswith(self.continuation_char):
                parts.append(line.replace(self.continuation_char, self.replace_char))
            else:
                parts.append(line)
                outlines.append(''.join(parts))
                parts = []

        return outlines


class WhitespaceSplitter(DefaultSplitter):
    def process_line(self, line):
        """Replace tab with space within ``line`` while respecting quoted substrings"""
        newline = []
        in_quote = False
        lastchar = None
        for char in line:
            if char == self.quotechar and (self.escapechar is None
                                           or lastchar != self.escapechar):
                in_quote = not in_quote
            if char == '\t' and not in_quote:
                char = ' '
            lastchar = char
            newline.append(char)

        return ''.join(newline)


extra_reader_pars = ('Reader', 'Inputter', 'Outputter',
                     'delimiter', 'comment', 'quotechar', 'header_start',
                     'data_start', 'data_end', 'converters', 'encoding',
                     'data_Splitter', 'header_Splitter',
                     'names', 'include_names', 'exclude_names', 'strict_names',
                     'fill_values', 'fill_include_names', 'fill_exclude_names')


def _get_reader(Reader, Inputter=None, Outputter=None, **kwargs):
    """Initialize a table reader allowing for common customizations.  See ui.get_reader()
    for param docs.  This routine is for internal (package) use only and is useful
    because it depends only on the "core" module.
    """

    from .fastbasic import FastBasic
    if issubclass(Reader, FastBasic):  # Fast readers handle args separately
        if Inputter is not None:
            kwargs['Inputter'] = Inputter
        return Reader(**kwargs)

    # If user explicitly passed a fast reader with enable='force'
    # (e.g. by passing non-default options), raise an error for slow readers
    if 'fast_reader' in kwargs:
        if kwargs['fast_reader']['enable'] == 'force':
            raise ParameterError('fast_reader required with '
                                 '{}, but this is not a fast C reader: {}'
                                 .format(kwargs['fast_reader'], Reader))
        else:
            del kwargs['fast_reader']  # Otherwise ignore fast_reader parameter

    reader_kwargs = dict([k, v] for k, v in kwargs.items() if k not in extra_reader_pars)
    reader = Reader(**reader_kwargs)

    if Inputter is not None:
        reader.inputter = Inputter()

    if Outputter is not None:
        reader.outputter = Outputter()

    # Issue #855 suggested to set data_start to header_start + default_header_length
    # Thus, we need to retrieve this from the class definition before resetting these numbers.
    try:
        default_header_length = reader.data.start_line - reader.header.start_line
    except TypeError:  # Start line could be None or an instancemethod
        default_header_length = None

    # csv.reader is hard-coded to recognise either '\r' or '\n' as end-of-line,
    # therefore DefaultSplitter cannot handle these as delimiters.
    if 'delimiter' in kwargs:
        if kwargs['delimiter'] in ('\n', '\r', '\r\n'):
            reader.header.splitter = BaseSplitter()
            reader.data.splitter = BaseSplitter()
        reader.header.splitter.delimiter = kwargs['delimiter']
        reader.data.splitter.delimiter = kwargs['delimiter']
    if 'comment' in kwargs:
        reader.header.comment = kwargs['comment']
        reader.data.comment = kwargs['comment']
    if 'quotechar' in kwargs:
        reader.header.splitter.quotechar = kwargs['quotechar']
        reader.data.splitter.quotechar = kwargs['quotechar']
    if 'data_start' in kwargs:
        reader.data.start_line = kwargs['data_start']
    if 'data_end' in kwargs:
        reader.data.end_line = kwargs['data_end']
    if 'header_start' in kwargs:
        if (reader.header.start_line is not None):
            reader.header.start_line = kwargs['header_start']
            # For FixedWidthTwoLine the data_start is calculated relative to the position line.
            # However, position_line is given as absolute number and not relative to header_start.
            # So, ignore this Reader here.
            if (('data_start' not in kwargs) and (default_header_length is not None)
                    and reader._format_name not in ['fixed_width_two_line', 'commented_header']):
                reader.data.start_line = reader.header.start_line + default_header_length
        elif kwargs['header_start'] is not None:
            # User trying to set a None header start to some value other than None
            raise ValueError('header_start cannot be modified for this Reader')
    if 'converters' in kwargs:
        reader.outputter.converters = kwargs['converters']
    if 'data_Splitter' in kwargs:
        reader.data.splitter = kwargs['data_Splitter']()
    if 'header_Splitter' in kwargs:
        reader.header.splitter = kwargs['header_Splitter']()
    if 'names' in kwargs:
        reader.names = kwargs['names']
        if None in reader.names:
            raise TypeError('Cannot have None for column name')
        if len(set(reader.names)) != len(reader.names):
            raise ValueError('Duplicate column names')
    if 'include_names' in kwargs:
        reader.include_names = kwargs['include_names']
    if 'exclude_names' in kwargs:
        reader.exclude_names = kwargs['exclude_names']
    # Strict names is normally set only within the guessing process to
    # indicate that column names cannot be numeric or have certain
    # characters at the beginning or end.  It gets used in
    # BaseHeader.check_column_names().
    if 'strict_names' in kwargs:
        reader.strict_names = kwargs['strict_names']
    if 'fill_values' in kwargs:
        reader.data.fill_values = kwargs['fill_values']
    if 'fill_include_names' in kwargs:
        reader.data.fill_include_names = kwargs['fill_include_names']
    if 'fill_exclude_names' in kwargs:
        reader.data.fill_exclude_names = kwargs['fill_exclude_names']
    if 'encoding' in kwargs:
        reader.encoding = kwargs['encoding']
        reader.inputter.encoding = kwargs['encoding']

    return reader


extra_writer_pars = ('delimiter', 'comment', 'quotechar', 'formats',
                     'strip_whitespace',
                     'names', 'include_names', 'exclude_names',
                     'fill_values', 'fill_include_names',
                     'fill_exclude_names')


def _get_writer(Writer, fast_writer, **kwargs):
    """Initialize a table writer allowing for common customizations. This
    routine is for internal (package) use only and is useful because it depends
    only on the "core" module."""

    from .fastbasic import FastBasic

    # A value of None for fill_values imply getting the default string
    # representation of masked values (depending on the writer class), but the
    # machinery expects a list.  The easiest here is to just pop the value off,
    # i.e. fill_values=None is the same as not providing it at all.
    if 'fill_values' in kwargs and kwargs['fill_values'] is None:
        del kwargs['fill_values']

    if issubclass(Writer, FastBasic):  # Fast writers handle args separately
        return Writer(**kwargs)
    elif fast_writer and f'fast_{Writer._format_name}' in FAST_CLASSES:
        # Switch to fast writer
        kwargs['fast_writer'] = fast_writer
        return FAST_CLASSES[f'fast_{Writer._format_name}'](**kwargs)

    writer_kwargs = dict([k, v] for k, v in kwargs.items() if k not in extra_writer_pars)
    writer = Writer(**writer_kwargs)

    if 'delimiter' in kwargs:
        writer.header.splitter.delimiter = kwargs['delimiter']
        writer.data.splitter.delimiter = kwargs['delimiter']
    if 'comment' in kwargs:
        writer.header.write_comment = kwargs['comment']
        writer.data.write_comment = kwargs['comment']
    if 'quotechar' in kwargs:
        writer.header.splitter.quotechar = kwargs['quotechar']
        writer.data.splitter.quotechar = kwargs['quotechar']
    if 'formats' in kwargs:
        writer.data.formats = kwargs['formats']
    if 'strip_whitespace' in kwargs:
        if kwargs['strip_whitespace']:
            # Restore the default SplitterClass process_val method which strips
            # whitespace.  This may have been changed in the Writer
            # initialization (e.g. Rdb and Tab)
            writer.data.splitter.process_val = operator.methodcaller('strip', ' \t')
        else:
            writer.data.splitter.process_val = None
    if 'names' in kwargs:
        writer.header.names = kwargs['names']
    if 'include_names' in kwargs:
        writer.include_names = kwargs['include_names']
    if 'exclude_names' in kwargs:
        writer.exclude_names = kwargs['exclude_names']
    if 'fill_values' in kwargs:
        # Prepend user-specified values to the class default.
        with suppress(TypeError, IndexError):
            # Test if it looks like (match, replace_string, optional_colname),
            # in which case make it a list
            kwargs['fill_values'][1] + ''
            kwargs['fill_values'] = [kwargs['fill_values']]
        writer.data.fill_values = kwargs['fill_values'] + writer.data.fill_values
    if 'fill_include_names' in kwargs:
        writer.data.fill_include_names = kwargs['fill_include_names']
    if 'fill_exclude_names' in kwargs:
        writer.data.fill_exclude_names = kwargs['fill_exclude_names']
    return writer
