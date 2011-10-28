""" Asciitable: an extensible ASCII table reader and writer.

core.py:
  Core base classes and functions for reading and writing tables.

:Copyright: Smithsonian Astrophysical Observatory (2010)
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

import os
import sys
import re
import csv
import itertools

try:
    import numpy
    has_numpy = True
except ImportError:
    has_numpy = False

class InconsistentTableError(ValueError):
    pass

# Python 3 compatibility tweaks.  Should work back through 2.4.
try:
    import cStringIO as io
except ImportError:
    import io

try:
    next = next
except NameError:
    next = lambda x: x.next()

try:
    izip = itertools.izip
except AttributeError:
    izip = zip

try:
    long = long
except NameError:
    long = int

try:
    unicode = unicode
except NameError:
    unicode = str

# Python 2.4 comptability: any() function is built-in only for 2.5 onward
try:
    any = any
except NameError:
    def any(vals):
        for val in vals:
            if val:
                return True
        return False

class Keyword(object):
    """Table keyword"""
    def __init__(self, name, value, units=None, comment=None, format=None):
        self.name = name
        self.value = value
        self.units = units
        self.comment = comment
        self.format = format

class NoType(object):
    pass

class StrType(NoType):
    pass

class NumType(NoType):
    pass

class FloatType(NumType):
    pass

class IntType(NumType):
    pass

class AllType(StrType, FloatType, IntType):
    pass

class Column(object):
    """Table column.

    The key attributes of a Column object are:

    * **name** : column name
    * **index** : column index (first column has index=0, second has index=1, etc)
    * **type** : column type (NoType, StrType, NumType, FloatType, IntType)
    * **str_vals** : list of column values as strings
    * **data** : list of converted column values
    """
    def __init__(self, name, index):
        self.name = name
        self.index = index
        self.type = NoType
        self.str_vals = []
        self.fill_values = {}
        self.formatter = None
    
    def __iter__(self):
        '''iterate over formated column values
        
        Each value is paased through self.formatter.
        If str(self.formatter(value)) is found in the fill_values specification,
        the corresponding fill_value is returned, otherwise the formated value.
        '''
        for val in self.data:
            yield self.fill_values.get(str(self.formatter(val)).strip(), self.formatter(val))

class BaseInputter(object):
    """Get the lines from the table input and return a list of lines.  The input table can be one of:

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
            if hasattr(table, 'read'):
                table = table.read()
            elif '\n' not in table and '\r' not in table + '':
                table = open(table, 'r').read()
            lines = table.splitlines()
        except TypeError:
            try:
                # See if table supports indexing, slicing, and iteration
                table[0]
                table[0:1]
                iter(table)
                lines = table
            except TypeError:
                raise TypeError('Input "table" must be a string (filename or data) or an iterable')

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
      splitter = asciitable.DefaultSplitter()
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
    
    def process_line(self, line):
        """Remove whitespace at the beginning or end of line.  This is especially useful for
        whitespace-delimited files to prevent spurious columns at the beginning or end. 
        If splitting on whitespace then replace unquoted tabs with space first"""
        if self.delimiter == '\s':
            line = _replace_tab_with_space(line, self.escapechar, self.quotechar)
        return line.strip()

    def __init__(self):
        self.csv_writer = None
        self.csv_writer_out = io.StringIO()

    def __call__(self, lines):
        """Return an iterator over the table ``lines``, where each iterator output
        is a list of the split line values.

        :param lines: list of table lines
        :returns: iterator
        """
        if self.process_line:
            lines = [self.process_line(x) for x in lines]

        if self.delimiter == '\s':
            delimiter = ' '
        else:
            delimiter = self.delimiter

        csv_reader = csv.reader(lines,
                                delimiter = delimiter,
                                doublequote = self.doublequote,
                                escapechar =self.escapechar,
                                quotechar = self.quotechar,
                                quoting = self.quoting,
                                skipinitialspace = self.skipinitialspace
                                )
        for vals in csv_reader:
            if self.process_val:
                yield [self.process_val(x) for x in vals]
            else:
                yield vals
            
    def join(self, vals):
        if self.delimiter is None:
            delimiter = ' '
        else:
            delimiter = self.delimiter

        if self.csv_writer is None:
            self.csv_writer = csv.writer(self.csv_writer_out,
                                         delimiter = self.delimiter,
                                         doublequote = self.doublequote,
                                         escapechar = self.escapechar,
                                         quotechar = self.quotechar,
                                         quoting = self.quoting,
                                         lineterminator = '',
                                         )
        self.csv_writer_out.seek(0)
        self.csv_writer_out.truncate()
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
    if hasattr(line_or_func, '__call__'):
        return line_or_func(lines)
    else:
        return line_or_func

class BaseHeader(object):
    """Base table header reader

    :param auto_format: format string for auto-generating column names
    :param start_line: None, int, or a function of ``lines`` that returns None or int
    :param comment: regular expression for comment lines
    :param splitter_class: Splitter class for splitting data lines into columns
    :param names: list of names corresponding to each data column
    :param include_names: list of names to include in output (default=None selects all names)
    :param exclude_names: list of names to exlude from output (applied after ``include_names``)
    """
    auto_format = 'col%d'
    start_line = None
    comment = None
    splitter_class = DefaultSplitter
    names = None
    include_names = None
    exclude_names = None
    write_spacer_lines = ['ASCIITABLE_WRITE_SPACER_LINE']

    def __init__(self):
        self.splitter = self.__class__.splitter_class()
       
    def _set_cols_from_names(self):
        # Filter full list of non-null column names with the include/exclude lists
        names = set(self.names)
        if self.include_names is not None:
            names.intersection_update(self.include_names)
        if self.exclude_names is not None:
            names.difference_update(self.exclude_names)
            
        self.cols = [Column(name=x, index=i) for i, x in enumerate(self.names) if x in names]

    def get_cols(self, lines):
        """Initialize the header Column objects from the table ``lines``.

        Based on the previously set Header attributes find or create the column names.
        Sets ``self.cols`` with the list of Columns.  This list only includes the actual
        requested columns after filtering by the include_names and exclude_names
        attributes.  See ``self.names`` for the full list.

        :param lines: list of table lines
        :returns: None
        """

        start_line = _get_line_index(self.start_line, self.process_lines(lines))
        if start_line is None:
            # No header line so auto-generate names from n_data_cols
            if self.names is None:
                # Get the data values from the first line of table data to determine n_data_cols
                try:
                    first_data_vals = next(self.data.get_str_vals())
                except StopIteration:
                    raise InconsistentTableError('No data lines found so cannot autogenerate column names')
                n_data_cols = len(first_data_vals)
                self.names = [self.auto_format % i for i in range(1, n_data_cols+1)]

        elif self.names is None:
            # No column names supplied so read them from header line in table.
            for i, line in enumerate(self.process_lines(lines)):
                if i == start_line:
                    break
            else: # No header line matching
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

    def write(self, lines):
        if self.start_line is not None:
            for i, spacer_line in izip(range(self.start_line),
                                       itertools.cycle(self.write_spacer_lines)):
                lines.append(spacer_line)
            lines.append(self.splitter.join([x.name for x in self.cols]))

    @property
    def colnames(self):
        """Return the column names of the table"""
        return tuple(col.name for col in self.cols)

    def _get_n_data_cols(self):
        """Return the number of expected data columns from data splitting.
        This is either explicitly set (typically for fixedwidth splitters)
        or set to self.names otherwise.
        """
        if not hasattr(self, '_n_data_cols'):
            self._n_data_cols = len(self.names)
        return self._n_data_cols

    def _set_n_data_cols(self, val):
        """Return the number of expected data columns from data splitting.
        """
        self._n_data_cols = val

    n_data_cols = property(_get_n_data_cols, _set_n_data_cols)

    def get_type_map_key(self, col):
        return col.raw_type

    def get_col_type(self, col):
        try:
            type_map_key = self.get_type_map_key(col)
            return self.col_type_map[type_map_key.lower()]
        except KeyError:
            raise ValueError('Unknown data type ""%s"" for column "%s"' % (
                    col.raw_type, col.name))


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
    write_spacer_lines = ['ASCIITABLE_WRITE_SPACER_LINE']
    formats = {}
    default_formatter = str
    fill_values = []
    fill_include_names = None
    fill_exclude_names = None
    
    def __init__(self):
        self.splitter = self.__class__.splitter_class()

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
            #if input is only one <fill_spec>, then make it a list
            try:
                self.fill_values[0] + ''
                self.fill_values = [ self.fill_values ] 
            except TypeError:
                pass
            # Step 1: Set the default list of columns which are affected by fill_values
            colnames = set(self.header.colnames)
            if self.fill_include_names is not None:
                colnames.intersection_update(self.fill_include_names)
            if self.fill_exclude_names is not None:
                colnames.difference_update(self.fill_exclude_names)
        
            # Step 2a: Find out which columns are affected by this tuple
            # iterate over reversed order, so last condition is set first and overwritten by earlier conditions
            for replacement in reversed(self.fill_values):
                if len(replacement) < 2:
                    raise ValueError("Format of fill_values must be (<bad>, <fill>, <optional col1>, ...)")
                elif len(replacement) == 2:
                    affect_cols = colnames
                else:
                    affect_cols = replacement[2:]
                
                for i, key in ((i, x) for i, x in enumerate(self.header.colnames) if x in affect_cols):
                    cols[i].fill_values[replacement[0]] = str(replacement[1])

    def _set_masks(self, cols):
        """Replace string values in col.str_vales and set masks"""
        if self.fill_values:
            for col in (col for col in cols if col.fill_values):
                col.mask = [False] * len(col.str_vals)
                for i, str_val in ((i, x) for i, x in enumerate(col.str_vals) if x in col.fill_values):
                    col.str_vals[i] = col.fill_values[str_val]
                    col.mask[i] = True

    def write(self, lines):
        if hasattr(self.start_line, '__call__'):
            raise TypeError('Start_line attribute cannot be callable for write()')
        else:
            data_start_line = self.start_line or 0

        while len(lines) < data_start_line:
            lines.append(itertools.cycle(self.write_spacer_lines))

        formatters = []
        for col in self.cols:
            formatter = self.formats.get(col.name, self.default_formatter)
            if not hasattr(formatter, '__call__'):
                formatter = _format_func(formatter)
            col.formatter = formatter
            
        for vals in izip(*self.cols):
            lines.append(self.splitter.join(vals))

def _format_func(format_str):
    def func(val):
        return format_str % val
    return func


class DictLikeNumpy(dict):
    """Provide minimal compatibility with numpy rec array API for BaseOutputter
    object::
      
      table = asciitable.read('mytable.dat', numpy=False)
      table.field('x')    # List of elements in column 'x'
      table.dtype.names   # get column names in order
      table[1]            # returns row 1 as a list
      table[1][2]         # 3nd column in row 1
      table['col1'][1]    # Row 1 in column col1
      for row_vals in table:  # iterate over table rows
          print row_vals  # print list of vals in each row

    """
    # To do: - add colnames property to set colnames and dtype.names as well.
    # - ordered dict?

    class Dtype(object):
        pass

    def __init__(self, *args, **kwargs):
        self.dtype = DictLikeNumpy.Dtype()
        dict.__init__(self, *args, **kwargs)

    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item + '')
        except TypeError:
            return [dict.__getitem__(self, x)[item] for x in self.dtype.names]

    def field(self, colname):
        return self[colname]

    def __len__(self):
        return len(list(self.values())[0])

    def __iter__(self):
        self.__index = 0
        return self

    def __next__(self):
        try:
            vals = self[self.__index]
        except IndexError:
            raise StopIteration
        else:
            self.__index += 1
            return vals

    def next(self):
        return self.__next__()
        
def convert_list(python_type):
    """Return a tuple ``(converter_func, converter_type)``.  The converter
    function converts a list into a list of the given ``python_type``.  This
    argument is a function that takes a single argument and returns a single
    value of the desired type.  In general this will be one of ``int``,
    ``float`` or ``str``.  The converter type is used to track the generic data
    type (int, float, str) that is produced by the converter function.
    """
    type_map = {int: IntType,
                float: FloatType,
                str: StrType}
    converter_type = type_map.get(python_type, AllType)

    def converter(vals):
        return [python_type(x) for x in vals]
    return converter, converter_type

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
    default_converters = [convert_list(int),
                          convert_list(float),
                          convert_list(str)]

    def __call__(self, cols):
        self._convert_vals(cols)
        table = DictLikeNumpy((x.name, x.data) for x in cols)
        table.dtype.names = tuple(x.name for x in cols)
        return table

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
                except IndexError:
                    raise ValueError('Column %s failed to convert' % col.name)

class NumpyOutputter(BaseOutputter):
    """Output the table as a numpy.rec.recarray

    Missing or bad data values are handled at two levels.  The first is in
    the data reading step where if ``data.fill_values`` is set then any
    occurences of a bad value are replaced by the correspond fill value.
    At the same time a boolean list ``mask`` is created in the column object.

    The second stage is when converting to numpy arrays which by default generates
    masked arrays, if ``data.fill_values`` is set and plain arrays if it is not.
    In the rare case that plain arrays are needed set ``auto_masked`` (default = True) and
    ``default_masked`` (default = False) to control this behavior as follows:

    ===========  ==============  ===========  ============
    auto_masked  default_masked  fill_values  output
    ===========  ==============  ===========  ============
    --            True           --           masked_array
    --            False          None         array
    True          --             dict(..)     masked_array
    False         --             dict(..)     array
    ===========  ==============  ===========  ============

    To set these values use::
    
      Outputter = asciitable.NumpyOutputter()
      Outputter.default_masked = True
    
    """

    auto_masked_array = True
    default_masked_array = False

    if has_numpy:
        default_converters = [convert_numpy(numpy.int),
                              convert_numpy(numpy.float),
                              convert_numpy(numpy.str)]

    def __call__(self, cols):
        self._convert_vals(cols)
        recarr = numpy.rec.fromarrays([x.data for x in cols], names=[x.name for x in cols])
        if self.default_masked_array or (self.auto_masked_array and
                                         any(col.fill_values for col in cols)):
            maarr = recarr.view(numpy.ma.MaskedArray)
            for col in cols:
                if col.fill_values:
                    maarr[col.name] = numpy.ma.masked_where(col.mask, maarr[col.name])
            return maarr
        else:
            return recarr


class BaseReader(object):
    """Class providing methods to read an ASCII table using the specified
    header, data, inputter, and outputter instances.

    Typical usage is to instantiate a Reader() object and customize the
    ``header``, ``data``, ``inputter``, and ``outputter`` attributes.  Each
    of these is an object of the corresponding class.

    There is one method ``inconsistent_handler`` that can be used to customize the
    behavior of ``read()`` in the event that a data row doesn't match the header.
    The default behavior is to raise an InconsistentTableError.

    """
    def __init__(self):
        self.header = BaseHeader()
        self.data = BaseData()
        self.inputter = BaseInputter()
        self.outputter = BaseOutputter()
        self.meta = {}                  # Placeholder for storing table metadata 
        self.keywords = []              # Placeholder for storing table Keywords
        # Data and Header instances benefit from a little cross-coupling.  Header may need to
        # know about number of data columns for auto-column name generation and Data may
        # need to know about header (e.g. for fixed-width tables where widths are spec'd in header.
        self.data.header = self.header
        self.header.data = self.data

    def read(self, table):
        """Read the ``table`` and return the results in a format determined by
        the ``outputter`` attribute.

        The ``table`` parameter is any string or object that can be processed
        by the instance ``inputter``.  For the base Inputter class ``table`` can be
        one of:

        * File name
        * String (newline separated) with all header and data lines (must have at least 2 lines)
        * List of strings

        :param table: table input
        :returns: output table
        """
        # If ``table`` is a file then store the name in the ``data``
        # attribute. The ``table`` is a "file" if it is a string
        # without the new line specific to the OS.
        try:
            if os.linesep not in table + '':
                self.data.table_name = os.path.basename(table)
        except TypeError:
            # Not a string.
            pass
            
        # Same from __init__.  ??? Do these need to be here?
        self.data.header = self.header
        self.header.data = self.data

        self.lines = self.inputter.get_lines(table)
        self.data.get_data_lines(self.lines)
        self.header.get_cols(self.lines)
        cols = self.header.cols         # header.cols corresponds to *output* columns requested
        n_data_cols = self.header.n_data_cols # number of data cols expected from splitter
        self.data.splitter.cols = cols

        for i, str_vals in enumerate(self.data.get_str_vals()):
            if len(str_vals) != n_data_cols:                
                str_vals = self.inconsistent_handler(str_vals, n_data_cols)
                
                #if str_vals is None, we skip this row
                if str_vals is None:
                    continue
                
                #otherwise, we raise an error only if it is still inconsistent
                if len(str_vals) != n_data_cols:
                    errmsg = ('Number of header columns (%d) inconsistent with '
                              'data columns (%d) at data line %d\n'
                              'Header values: %s\n'
                              'Data values: %s' % (len(cols), len(str_vals), i,
                                                   [x.name for x in cols], str_vals))
                    raise InconsistentTableError(errmsg)

            for col in cols:
                col.str_vals.append(str_vals[col.index])

        self.data.masks(cols)
        self.table = self.outputter(cols)
        self.cols = self.header.cols

        return self.table
    
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
        #an empty list will always trigger an InconsistentTableError in read()
        return str_vals

    @property
    def comment_lines(self):
        """Return lines in the table that match header.comment regexp"""
        if not hasattr(self, 'lines'):
            raise ValueError('Table must be read prior to accessing the header_comment_lines')
        if self.header.comment:
            re_comment = re.compile(self.header.comment)
            comment_lines = [x for x in self.lines if re_comment.match(x)]
        else:
            comment_lines = []
        return comment_lines

    def write(self, table=None):
        """Write ``table`` as list of strings.

        :param table: asciitable Reader object
        :returns: list of strings corresponding to ASCII table
        """
        if table is None:
            table = self

        # link information about the columns to the writer object (i.e. self)
        self.header.cols = table.cols
        self.data.cols = self.header.cols
        self.data.masks(self.data.cols)

        # Write header and data to lines list 
        lines = []
        self.header.write(lines)
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

    def process_lines(self, lines):
        striplines = (x.strip() for x in lines)
        lines = [x for x in striplines if len(x) > 0]

        parts = []
        outlines = []
        for line in lines:
            if line.endswith(self.continuation_char):
                parts.append(line.rstrip(self.continuation_char))
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
                     'names', 'include_names', 'exclude_names',
                     'fill_values', 'fill_include_names', 'fill_exclude_names')

def _get_reader(Reader, Inputter=None, Outputter=None, numpy=True, **kwargs):
    """Initialize a table reader allowing for common customizations.  See ui.get_reader()
    for param docs.  This routine is for internal (package) use only and is useful
    because it depends only on the "core" module.
    """

    reader_kwargs = dict([k, v] for k, v in kwargs.items() if k not in extra_reader_pars)
    reader = Reader(**reader_kwargs)

    if Inputter is not None:
        reader.inputter = Inputter()

    if has_numpy and numpy:
        reader.outputter = NumpyOutputter()

    if Outputter is not None:
        reader.outputter = Outputter()

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
        reader.header.start_line = kwargs['header_start']
    if 'converters' in kwargs:
        reader.outputter.converters = kwargs['converters']
    if 'data_Splitter' in kwargs:
        reader.data.splitter = kwargs['data_Splitter']()
    if 'header_Splitter' in kwargs:
        reader.header.splitter = kwargs['header_Splitter']()
    if 'names' in kwargs:
        reader.header.names = kwargs['names']
    if 'include_names' in kwargs:
        reader.header.include_names = kwargs['include_names']
    if 'exclude_names' in kwargs:
        reader.header.exclude_names = kwargs['exclude_names']
    if 'fill_values' in kwargs:
        reader.data.fill_values = kwargs['fill_values']
    if 'fill_include_names' in kwargs:
        reader.data.fill_include_names = kwargs['fill_include_names']
    if 'fill_exclude_names' in kwargs:
        reader.data.fill_exclude_names = kwargs['fill_exclude_names']

    return reader

extra_writer_pars = ('delimiter', 'comment', 'quotechar', 'formats',
                     'names', 'include_names', 'exclude_names', 
                     'fill_values', 'fill_include_names',
                     'fill_exclude_names')

def _get_writer(Writer, **kwargs):
    """Initialize a table writer allowing for common customizations. This
    routine is for internal (package) use only and is useful because it depends
    only on the "core" module. """

    writer_kwargs = dict([k, v] for k, v in kwargs.items() if k not in extra_writer_pars)
    writer = Writer(**writer_kwargs)

    if 'delimiter' in kwargs:
        writer.header.splitter.delimiter = kwargs['delimiter']
        writer.data.splitter.delimiter = kwargs['delimiter']
    if 'write_comment' in kwargs:
        writer.header.write_comment = kwargs['write_comment']
        writer.data.write_comment = kwargs['write_comment']
    if 'quotechar' in kwargs:
        writer.header.splitter.quotechar = kwargs['quotechar']
        writer.data.splitter.quotechar = kwargs['quotechar']
    if 'formats' in kwargs:
        writer.data.formats = kwargs['formats']
    if 'names' in kwargs:
        writer.header.names = kwargs['names']
    if 'include_names' in kwargs:
        writer.header.include_names = kwargs['include_names']
    if 'exclude_names' in kwargs:
        writer.header.exclude_names = kwargs['exclude_names']
    if 'fill_values' in kwargs:
        writer.data.fill_values = kwargs['fill_values']
    if 'fill_include_names' in kwargs:
        writer.data.fill_include_names = kwargs['fill_include_names']
    if 'fill_exclude_names' in kwargs:
        writer.data.fill_exclude_names = kwargs['fill_exclude_names']
    return writer

    
