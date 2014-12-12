# Licensed under a 3-clause BSD style license - see LICENSE.rst
""" An extensible ASCII table reader and writer.

core.py:
  Core base classes and functions for reading and writing tables.

:Copyright: Smithsonian Astrophysical Observatory (2010)
:Author: Tom Aldcroft (aldcroft@head.cfa.harvard.edu)
"""

from __future__ import absolute_import, division, print_function

import os
import re
import csv
import itertools
import functools
import numpy
import warnings
import copy

from ...extern import six
from ...extern.six.moves import zip
from ...extern.six.moves import cStringIO as StringIO
from ...utils.exceptions import AstropyWarning

from ...table import Table, col_getattr, col_setattr, col_iter_str_vals
from ...utils.compat import ignored
from ...utils.data import get_readable_fileobj
from ...utils import OrderedDict
from . import connect

# Global dictionary mapping format arg to the corresponding Reader class
FORMAT_CLASSES = {}

# Similar dictionary for fast readers
FAST_CLASSES = {}

class MaskedConstant(numpy.ma.core.MaskedConstant):
    """A trivial extension of numpy.ma.masked

    We want to be able to put the generic term ``masked`` into a dictionary.
    In python 2.7 we can just use ``numpy.ma.masked``, but in python 3.1 and 3.2 that
    is not hashable, see https://github.com/numpy/numpy/issues/4660
    So, we need to extend it here with a hash value.
    """
    def __hash__(self):
        '''All instances of this class shall have the same hash.'''
        # Any large number will do.
        return 1234567890

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

class NoType(object):
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


class Column(object):
    """Table column.

    The key attributes of a Column object are:

    * **name** : column name
    * **type** : column type (NoType, StrType, NumType, FloatType, IntType)
    * **str_vals** : list of column values as strings
    * **data** : list of converted column values
    """
    def __init__(self, name):
        self.name = name
        self.type = NoType
        self.str_vals = []
        self.fill_values = {}


class BaseInputter(object):
    """
    Get the lines from the table input and return a list of lines.  The input
    table can be one of:

    * File name
    * String (newline separated) with all header and data lines (must have at least 2 lines)
    * File-like object with read() method
    * List of strings
    """
    def get_lines(self, table):
        """Get the lines from the ``table`` input.

        :param table: table input
        :returns: list of lines
        """
        try:
            if (hasattr(table, 'read') or
                    ('\n' not in table + '' and '\r' not in table + '')):
                with get_readable_fileobj(table) as file_obj:
                    table = file_obj.read()
            lines = table.splitlines()
        except TypeError:
            try:
                # See if table supports indexing, slicing, and iteration
                table[0]
                table[0:1]
                iter(table)
                lines = table
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


class BaseSplitter(object):
    """Base splitter that uses python's split method to do the work.

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

    :param delimiter: one-character string used to separate fields
    """
    delimiter = None

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

    :param delimiter: one-character string used to separate fields.
    :param doublequote:  control how instances of *quotechar* in a field are quoted
    :param escapechar: character to remove special meaning from following character
    :param quotechar: one-character stringto quote fields containing special characters
    :param quoting: control when quotes are recognised by the reader
    :param skipinitialspace: ignore whitespace immediately following the delimiter
    """
    delimiter = ' '
    quotechar = '"'
    doublequote = True
    escapechar = None
    quoting = csv.QUOTE_MINIMAL
    skipinitialspace = True
    csv_writer = None
    csv_writer_out = StringIO()


    def process_line(self, line):
        """Remove whitespace at the beginning or end of line.  This is especially useful for
        whitespace-delimited files to prevent spurious columns at the beginning or end.
        If splitting on whitespace then replace unquoted tabs with space first"""
        if self.delimiter == '\s':
            line = _replace_tab_with_space(line, self.escapechar, self.quotechar)
        return line.strip()


    def __call__(self, lines):
        """Return an iterator over the table ``lines``, where each iterator output
        is a list of the split line values.

        :param lines: list of table lines
        :returns: iterator
        """
        if self.process_line:
            lines = [self.process_line(x) for x in lines]

        # In Python 2.x the inputs to csv cannot be unicode.  In Python 3 these
        # lines do nothing.
        escapechar = None if self.escapechar is None else str(self.escapechar)
        quotechar = None if self.quotechar is None else str(self.quotechar)
        delimiter = None if self.delimiter is None else str(self.delimiter)

        if delimiter == '\s':
            delimiter = ' '

        csv_reader = csv.reader(lines,
                                delimiter=delimiter,
                                doublequote=self.doublequote,
                                escapechar=escapechar,
                                quotechar=quotechar,
                                quoting=self.quoting,
                                skipinitialspace=self.skipinitialspace
                                )
        for vals in csv_reader:
            if self.process_val:
                yield [self.process_val(x) for x in vals]
            else:
                yield vals

    def join(self, vals):

        # In Python 2.x the inputs to csv cannot be unicode
        escapechar = None if self.escapechar is None else str(self.escapechar)
        quotechar = None if self.quotechar is None else str(self.quotechar)
        delimiter = ' ' if self.delimiter is None else str(self.delimiter)

        if self.csv_writer is None:
            self.csv_writer = csv.writer(self.csv_writer_out,
                                         delimiter=delimiter,
                                         doublequote=self.doublequote,
                                         escapechar=escapechar,
                                         quotechar=quotechar,
                                         quoting=self.quoting,
                                         lineterminator='',
                                         )
        self.csv_writer_out.seek(0)
        self.csv_writer_out.truncate()
        if self.process_val:
            vals = [self.process_val(x) for x in vals]
        self.csv_writer.writerow(vals)

        return self.csv_writer_out.getvalue()


def _replace_tab_with_space(line, escapechar, quotechar):
    """Replace tab with space within ``line`` while respecting quoted substrings"""
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


class BaseHeader(object):
    """Base table header reader

    :param auto_format: format string for auto-generating column names
    :param start_line: None, int, or a function of ``lines`` that returns None or int
    :param comment: regular expression for comment lines
    :param splitter_class: Splitter class for splitting data lines into columns
    :param names: list of names corresponding to each data column
    """
    auto_format = 'col%d'
    start_line = None
    comment = None
    splitter_class = DefaultSplitter
    names = None
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

        :param lines: list of table lines
        :returns: None
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
            self.names = [self.auto_format % i for i in range(1, n_data_cols + 1)]

        else:
            for i, line in enumerate(self.process_lines(lines)):
                if i == start_line:
                    break
            else:  # No header line matching
                raise ValueError('No header line found in table')

            self.names = next(self.splitter([line]))

        self._set_cols_from_names()

    def process_lines(self, lines):
        """Generator to yield non-comment lines"""
        if self.comment:
            re_comment = re.compile(self.comment)
        # Yield non-comment lines
        for line in lines:
            if line and (not self.comment or not re_comment.match(line)):
                yield line

    def write_comments(self, lines, meta):
        if self.write_comment is not False:
            for comment in meta.get('comments', []):
                lines.append(self.write_comment + comment)

    def write(self, lines):
        if self.start_line is not None:
            for i, spacer_line in zip(range(self.start_line),
                                      itertools.cycle(self.write_spacer_lines)):
                lines.append(spacer_line)
            lines.append(self.splitter.join([col_getattr(x, 'name') for x in self.cols]))

    @property
    def colnames(self):
        """Return the column names of the table"""
        return tuple(col.name if isinstance(col, Column) else col_getattr(col, 'name')
                     for col in self.cols)

    def get_type_map_key(self, col):
        return col.raw_type

    def get_col_type(self, col):
        try:
            type_map_key = self.get_type_map_key(col)
            return self.col_type_map[type_map_key.lower()]
        except KeyError:
            raise ValueError('Unknown data type ""%s"" for column "%s"' % (
                col.raw_type, col.name))

    def check_column_names(self, names, strict_names, guessing):
        """Check column names.

        This must be done before applying the names transformation
        so that guessing will fail appropriately if `names` is supplied.
        For instance if the basic reader is given a table with no column header
        row.

        :param names: user-supplied list of column names
        :param strict_names: whether to impose extra requirements on names
        """
        if strict_names:
            # Impose strict requirements on column names (normally used in guessing)
            bads = [" ", ",", "|", "\t", "'", '"']
            for name in self.colnames:
                if (_is_number(name) or
                    len(name) == 0 or
                    name[0] in bads or
                    name[-1] in bads):
                    raise ValueError('Column name {0!r} does not meet strict name requirements'
                                     .format(name))
        # When guessing require at least two columns
        if guessing and len(self.colnames) <= 1:
            raise ValueError

        if names is not None and len(names) != len(self.colnames):
            raise ValueError('Length of names argument ({0}) does not match number'
                             ' of table columns ({1})'.format(len(names), len(self.colnames)))


class BaseData(object):
    """Base table data reader.

    :param start_line: None, int, or a function of ``lines`` that returns None or int
    :param end_line: None, int, or a function of ``lines`` that returns None or int
    :param comment: Regular expression for comment lines
    :param splitter_class: Splitter class for splitting data lines into columns
    """
    start_line = None
    end_line = None
    comment = None
    splitter_class = DefaultSplitter
    write_spacer_lines = ['ASCII_TABLE_WRITE_SPACER_LINE']
    fill_include_names = None
    fill_exclude_names = None
    # Currently, the default matches the numpy default for masked values.
    fill_values = [(masked, '--')]
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
        """Strip out comment lines and blank lines from list of ``lines``

        :param lines: all lines in table
        :returns: list of lines
        """
        nonblank_lines = (x for x in lines if x.strip())
        if self.comment:
            re_comment = re.compile(self.comment)
            return [x for x in nonblank_lines if not re_comment.match(x)]
        else:
            return [x for x in nonblank_lines]

    def get_data_lines(self, lines):
        """Set the ``data_lines`` attribute to the lines slice comprising the
        table data values."""
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
        """Set fill value for each column and then apply that fill value

        In the first step it is evaluated with value from ``fill_values`` applies to
        which column using ``fill_include_names`` and ``fill_exclude_names``.
        In the second step all replacements are done for the appropriate columns.
        """
        if self.fill_values:
            self._set_fill_values(cols)
            self._set_masks(cols)

    def _set_fill_values(self, cols):
        """Set the fill values of the individual cols based on fill_values of BaseData

        fill values has the following form:
        <fill_spec> = (<bad_value>, <fill_value>, <optional col_name>...)
        fill_values = <fill_spec> or list of <fill_spec>'s

        """
        if self.fill_values:
            # when we write tables the columns may be astropy.table.Columns
            # which don't carry a fill_values by default
            for col in cols:
                if ~hasattr(col, 'fill_values'):
                    col.fill_values = {}

            # if input is only one <fill_spec>, then make it a list
            with ignored(TypeError):
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
        """Replace string values in col.str_vals and set masks"""
        if self.fill_values:
            for col in (col for col in cols if col.fill_values):
                col.mask = numpy.zeros(len(col.str_vals), dtype=numpy.bool)
                for i, str_val in ((i, x) for i, x in enumerate(col.str_vals)
                                   if x in col.fill_values):
                    col.str_vals[i] = col.fill_values[str_val]
                    col.mask[i] = True

    def _replace_vals(self, cols):
        """Replace string values in col.str_vals"""
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
        '''convert all values in table to a list of lists of strings'''
        self._set_fill_values(self.cols)
        self._set_col_formats()
        for col in self.cols:
            col.str_vals = list(col_iter_str_vals())
        self._replace_vals(self.cols)
        return [col.str_vals for col in self.cols]

    def write(self, lines):
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
        """
        """
        for col in self.cols:
            if col_getattr(col, 'name') in self.formats:
                col_setattr(col, 'format', self.formats[col.name])


def convert_numpy(numpy_type):
    """Return a tuple ``(converter_func, converter_type)``.  The converter
    function converts a list into a numpy array of the given ``numpy_type``.
    This type must be a valid `numpy type
    <http://docs.scipy.org/doc/numpy/user/basics.types.html>`_, e.g.
    numpy.int, numpy.uint, numpy.int8, numpy.int64, numpy.float, numpy.float64,
    numpy.str.  The converter type is used to track the generic data type (int,
    float, str) that is produced by the converter function.
    """

    # Infer converter type from an instance of numpy_type.
    type_name = numpy.array([], dtype=numpy_type).dtype.name
    if 'int' in type_name:
        converter_type = IntType
    elif 'float' in type_name:
        converter_type = FloatType
    elif 'str' in type_name:
        converter_type = StrType
    else:
        converter_type = AllType

    def converter(vals):
        return numpy.array(vals, numpy_type)
    return converter, converter_type


class BaseOutputter(object):
    """Output table as a dict of column objects keyed on column name.  The
    table data are stored as plain python lists within the column objects.
    """
    converters = {}
    # Derived classes must define default_converters and __call__

    @staticmethod
    def _validate_and_copy(col, converters):
        """Validate the format for the type converters and then copy those
        which are valid converters for this column (i.e. converter type is
        a subclass of col.type)"""
        converters_out = []
        try:
            for converter in converters:
                converter_func, converter_type = converter
                if not issubclass(converter_type, NoType):
                    raise ValueError()
                if issubclass(converter_type, col.type):
                    converters_out.append((converter_func, converter_type))

        except (ValueError, TypeError):
            raise ValueError('Error: invalid format for converters, see documentation\n%s' %
                             converters)
        return converters_out

    def _convert_vals(self, cols):
        for col in cols:
            converters = self.converters.get(col.name,
                                             self.default_converters)
            col.converters = self._validate_and_copy(col, converters)

            while not hasattr(col, 'data'):
                try:
                    converter_func, converter_type = col.converters[0]
                    if not issubclass(converter_type, col.type):
                        raise TypeError()
                    col.data = converter_func(col.str_vals)
                    col.type = converter_type
                except (TypeError, ValueError):
                    col.converters.pop(0)
                except OverflowError:
                    # Overflow during conversion (most likely an int that doesn't fit in native C long).
                    # Put string at the top of the converters list for the next while iteration.
                    warnings.warn("OverflowError converting to {0} for column {1}, using string instead."
                                  .format(converter_type.__name__, col.name), AstropyWarning)
                    col.converters.insert(0, convert_numpy(numpy.str))
                except IndexError:
                    raise ValueError('Column %s failed to convert' % col.name)


class TableOutputter(BaseOutputter):
    """
    Output the table as an astropy.table.Table object.
    """

    default_converters = [convert_numpy(numpy.int),
                          convert_numpy(numpy.float),
                          convert_numpy(numpy.str)]

    def __call__(self, cols, meta):
        self._convert_vals(cols)

        # If there are any values that were filled and tagged with a mask bit then this
        # will be a masked table.  Otherwise use a plain table.
        masked = any(hasattr(col, 'mask') and numpy.any(col.mask) for col in cols)

        out = Table([x.data for x in cols], names=[x.name for x in cols], masked=masked,
                    meta=meta['table'])
        for col, out_col in zip(cols, out.columns.values()):
            if masked and hasattr(col, 'mask'):
                out_col.data.mask = col.mask
            for attr in ('format', 'unit', 'description'):
                if hasattr(col, attr):
                    setattr(out_col, attr, getattr(col, attr))

        return out


class MetaBaseReader(type):
    def __init__(cls, name, bases, dct):
        super(MetaBaseReader, cls).__init__(name, bases, dct)

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
            connect.io_registry.register_reader(io_format, Table, func)

            if dct.get('_io_registry_can_write', True):
                func = functools.partial(connect.io_write, io_format)
                connect.io_registry.register_writer(io_format, Table, func)


def _is_number(x):
    with ignored(ValueError):
        x = float(x)
        return True
    return False

def _apply_include_exclude_names(table, names, include_names, exclude_names):
    """Apply names, include_names and exclude_names to a table.

    :param table: input table (Reader object, NumPy struct array, list of lists, etc)
    :param names: list of names to override those in table (default=None uses existing names)
    :param include_names: list of names to include in output (default=None selects all names)
    :param exclude_names: list of names to exlude from output (applied after ``include_names``)
    """

    if names is not None:
        # Rename table column names to those passed by user
        # Temporarily rename with names that are not in `names` or `table.colnames`.
        # This ensures that rename succeeds regardless of existing names.
        xxxs = 'x' * max(len(name) for name in list(names) + list(table.colnames))
        for ii, colname in enumerate(table.colnames):
            table.rename_column(colname, xxxs + str(ii))

        for ii, name in enumerate(names):
            table.rename_column(xxxs + str(ii), name)

    names = set(table.colnames)
    if include_names is not None:
        names.intersection_update(include_names)
    if exclude_names is not None:
        names.difference_update(exclude_names)
    if names != set(table.colnames):
        remove_names = set(table.colnames) - set(names)
        table.remove_columns(remove_names)


@six.add_metaclass(MetaBaseReader)
class BaseReader(object):
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

    header_class = BaseHeader
    data_class = BaseData
    inputter_class = BaseInputter
    outputter_class = TableOutputter

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

        :param table: table input
        :returns: output table
        """
        # If ``table`` is a file then store the name in the ``data``
        # attribute. The ``table`` is a "file" if it is a string
        # without the new line specific to the OS.
        with ignored(TypeError):
            # Strings only
            if os.linesep not in table + '':
                self.data.table_name = os.path.basename(table)

        # Get a list of the lines (rows) in the table
        self.lines = self.inputter.get_lines(table)

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
                    errmsg = ('Number of header columns (%d) inconsistent with '
                              'data columns (%d) at data line %d\n'
                              'Header values: %s\n'
                              'Data values: %s' % (n_cols, len(str_vals), i,
                                                   [x.name for x in cols], str_vals))
                    raise InconsistentTableError(errmsg)

            for j, col in enumerate(cols):
                col.str_vals.append(str_vals[j])

        self.data.masks(cols)
        table = self.outputter(cols, self.meta)

        _apply_include_exclude_names(table, self.names, self.include_names, self.exclude_names)

        return table

    def inconsistent_handler(self, str_vals, ncols):
        """Adjust or skip data entries if a row is inconsistent with the header.

        The default implementation does no adjustment, and hence will always trigger
        an exception in read() any time the number of data entries does not match
        the header.

        Note that this will *not* be called if the row already matches the header.

        :param str_vals: A list of value strings from the current row of the table.
        :param ncols: The expected number of entries from the table header.
        :returns:
            list of strings to be parsed into data entries in the output table. If
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

    def write_header(self, lines, meta):
        self.header.write_comments(lines, meta)
        self.header.write(lines)

    def write(self, table):
        """Write ``table`` as list of strings.

        :param table: input table data (astropy.table.Table object)
        :returns: list of strings corresponding to ASCII table
        """

        # Check column names before altering
        self.header.cols = list(six.itervalues(table.columns))
        self.header.check_column_names(self.names, self.strict_names, False)

        _apply_include_exclude_names(table, self.names, self.include_names, self.exclude_names)

        # Now use altered columns
        new_cols = list(six.itervalues(table.columns))
        # link information about the columns to the writer object (i.e. self)
        self.header.cols = new_cols
        self.data.cols = new_cols

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
            if char == self.quotechar and (self.escapechar is None or
                                           lastchar != self.escapechar):
                in_quote = not in_quote
            if char == '\t' and not in_quote:
                char = ' '
            lastchar = char
            newline.append(char)

        return ''.join(newline)

extra_reader_pars = ('Reader', 'Inputter', 'Outputter',
                     'delimiter', 'comment', 'quotechar', 'header_start',
                     'data_start', 'data_end', 'converters',
                     'data_Splitter', 'header_Splitter',
                     'names', 'include_names', 'exclude_names', 'strict_names',
                     'fill_values', 'fill_include_names', 'fill_exclude_names')


def _get_reader(Reader, Inputter=None, Outputter=None, **kwargs):
    """Initialize a table reader allowing for common customizations.  See ui.get_reader()
    for param docs.  This routine is for internal (package) use only and is useful
    because it depends only on the "core" module.
    """

    from .fastbasic import FastBasic
    if issubclass(Reader, FastBasic): # Fast readers handle args separately
        if Inputter is not None:
            kwargs['Inputter'] = Inputter
        return Reader(**kwargs)

    if 'fast_reader' in kwargs:
        del kwargs['fast_reader'] # ignore fast_reader parameter for slow readers
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

    if 'delimiter' in kwargs:
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

    return reader

extra_writer_pars = ('delimiter', 'comment', 'quotechar', 'formats',
                     'strip_whitespace',
                     'names', 'include_names', 'exclude_names',
                     'fill_values', 'fill_include_names',
                     'fill_exclude_names')


def _get_writer(Writer, fast_writer, **kwargs):
    """Initialize a table writer allowing for common customizations. This
    routine is for internal (package) use only and is useful because it depends
    only on the "core" module. """

    from .fastbasic import FastBasic

    if issubclass(Writer, FastBasic): # Fast writers handle args separately
        return Writer(**kwargs)
    elif fast_writer and 'fast_{0}'.format(Writer._format_name) in FAST_CLASSES:
        # Switch to fast writer
        kwargs['fast_writer'] = fast_writer
        return FAST_CLASSES['fast_{0}'.format(Writer._format_name)](**kwargs)

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
            writer.data.splitter.process_val = lambda x: x.strip()
        else:
            writer.data.splitter.process_val = None
    if 'names' in kwargs:
        writer.header.names = kwargs['names']
    if 'include_names' in kwargs:
        writer.include_names = kwargs['include_names']
    if 'exclude_names' in kwargs:
        writer.exclude_names = kwargs['exclude_names']
    if 'fill_values' in kwargs:
        writer.data.fill_values = kwargs['fill_values']
    if 'fill_include_names' in kwargs:
        writer.data.fill_include_names = kwargs['fill_include_names']
    if 'fill_exclude_names' in kwargs:
        writer.data.fill_exclude_names = kwargs['fill_exclude_names']
    return writer
