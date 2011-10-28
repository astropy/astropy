"""Asciitable: an extensible ASCII table reader and writer.

memory.py:
  Classes to read table from in-memory data structure into
  asciitable format.  This is used for writing tables.

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

import asciitable.core as core
import asciitable.basic as basic
from asciitable.core import io, next, izip, any
if core.has_numpy:
    import numpy

class Memory(core.BaseReader):
    """Read a table from a data object in memory.  Several input data formats are supported:

    **Output of asciitable.read()**::

      table = asciitable.get_reader(Reader=asciitable.Daophot)
      data = table.read('t/daophot.dat')
      mem_data_from_table = asciitable.read(table, Reader=asciitable.Memory)
      mem_data_from_data = asciitable.read(data, Reader=asciitable.Memory)

    **Numpy structured array**::

      data = numpy.zeros((2,), dtype=[('col1','i4'), ('col2','f4'), ('col3', 'a10')])
      data[:] = [(1, 2., 'Hello'), (2, 3., "World")]
      mem_data = asciitable.read(data, Reader=asciitable.Memory)
      
    **Numpy masked structured array**::

      data = numpy.ma.zeros((2,), dtype=[('col1','i4'), ('col2','f4'), ('col3', 'a10')])
      data[:] = [(1, 2., 'Hello'), (2, 3., "World")]
      data['col2'] = ma.masked
      mem_data = asciitable.read(data, Reader=asciitable.Memory)
      
      In the current version all masked values will be converted to nan.

    **Sequence of sequences**::
    
      data = [[1, 2,   3      ],
              [4, 5.2, 6.1    ],
              [8, 9,   'hello']]
      mem_data = asciitable.read(data, Reader=asciitable.Memory, names=('c1','c2','c3'))

    **Dict of sequences**::

      data = {'c1': [1, 2, 3],
              'c2': [4, 5.2, 6.1],
              'c3': [8, 9, 'hello']}
      mem_data = asciitable.read(data, Reader=asciitable.Memory, names=('c1','c2','c3'))

    """
    def __init__(self):
        self.header = MemoryHeader()
        self.data = MemoryData()
        self.inputter = MemoryInputter()
        self.outputter = core.BaseOutputter()
        self.meta = {}                  # Placeholder for storing table metadata 
        self.keywords = []              # Table Keywords

    def read(self, table):
        self.data.header = self.header
        self.header.data = self.data

        self.lines = self.inputter.get_lines(table, self.header.names)
        self.data.get_data_lines(self.lines)
        self.header.get_cols(self.lines)
        cols = self.header.cols         # header.cols corresponds to *output* columns requested
        n_data_cols = len(self.header.names) # header.names corresponds to *all* header columns in table
        self.data.splitter.cols = cols

        for i, str_vals in enumerate(self.data.get_str_vals()):
            if len(list(str_vals)) != n_data_cols:
                errmsg = ('Number of header columns (%d) inconsistent with '
                          'data columns (%d) at data line %d\n'
                          'Header values: %s\n'
                          'Data values: %s' % (len(cols), len(str_vals), i,
                                               [x.name for x in cols], str_vals))
                raise core.InconsistentTableError(errmsg)

            for col in cols:
                col.str_vals.append(str_vals[col.index])

        self.data.masks(cols)
        self.cols = cols
        if hasattr(table, 'keywords'):
            self.keywords = table.keywords

        self.outputter.default_converters = [((lambda vals: vals), core.IntType),
                                             ((lambda vals: vals), core.FloatType),
                                             ((lambda vals: vals), core.StrType)]
        self.table = self.outputter(cols)
        self.cols = self.header.cols

        return self.table

    def write(self, table=None):
        """Not available for the Memory class (raises NotImplementedError)"""
        raise NotImplementedError

MemoryReader = Memory

class MemoryInputter(core.BaseInputter):
    """Get the lines from the table input and return an iterable object that contains the data lines.

    The input table can be one of:

    * asciitable Reader object
    * NumPy structured array
    * List of lists
    * Dict of lists
    """
    def get_lines(self, table, names):
        """Get the lines from the ``table`` input.
        
        :param table: table input
        :param names: list of column names (only used for dict of lists to set column order)
        :returns: list of lines

        """
        try:  
            # If table is dict-like this will return the first key.
            # If table is list-like this will return the first row.
            first_row_or_key = next(iter(table))
        except TypeError:
            # Not iterable, is it an asciitable Reader instance?
            if isinstance(table, core.BaseReader):
                lines = table.table
            else:
                # None of the possible choices so raise exception
                raise TypeError('Input table must be iterable or else be a Reader object')
        else:
            # table is iterable, now see if it is dict-like or list-like
            try:
                # If first_row_or_key is a row (in the case that table is
                # list-like) then this will raise exception
                table[first_row_or_key]
            except (TypeError, IndexError, ValueError):
                # Table is list-like (python list-of-lists or numpy recarray)
                lines = table
            else:
                # Table is dict-like.  Turn this into a DictLikeNumpy that has
                # an API similar to a numpy recarray.
                lines = core.DictLikeNumpy(table)
                if names is None:
                    lines.dtype.names = sorted(lines.keys())
                else:
                    lines.dtype.names = names

        # ``lines`` could now be one of the following iterable objects:
        # - NumPy recarray
        # - DictLikeNumpy object
        # - Python list of lists
        return lines

def get_val_type(val):
    """Get the asciitable data type corresponding to ``val``.  Try a series
    of possibilities, organized roughly by expected frequencies of types in
    data tables.
    """
    # Try native python types
    if isinstance(val, float):
        return core.FloatType
    elif isinstance(val, int):
        return core.IntType
    elif isinstance(val, str):
        return core.StrType
    elif isinstance(val, core.long):
        return core.IntType
    elif isinstance(val, core.unicode):
        return core.StrType
        
    # Not a native Python type so try a NumPy type
    try:
        type_name = val.dtype.name
    except AttributeError:
        pass
    else:
        if 'int' in type_name:
            return core.IntType
        elif 'float' in type_name:
            return core.FloatType
        elif 'string' in type_name:
            return core.StrType

    # Nothing matched
    raise TypeError("Memory: could not infer type for data value '%s'" % val)
    
def get_lowest_type(type_set):
    """Return the lowest common denominator among a set of asciitable Types,
    in order StrType, FloatType, IntType.  
    """
    if core.StrType in type_set:
        return core.StrType
    elif core.FloatType in type_set:
        return core.FloatType
    elif core.IntType in type_set:
        return core.IntType

    raise ValueError("Type_set '%s' does not have expected values" % type_set)
        

class MemoryHeader(core.BaseHeader):
    """Memory table header reader"""
    def __init__(self):
        pass

    def get_cols(self, lines):
        """Initialize the header Column objects from the table ``lines``.

        Based on the previously set Header attributes find or create the column names.
        Sets ``self.cols`` with the list of Columns.  This list only includes the actual
        requested columns after filtering by the include_names and exclude_names
        attributes.  See ``self.names`` for the full list.

        :param lines: list of table lines
        :returns: list of table Columns
        """

        if self.names is None:
            # No column names supplied so first try to get them from NumPy structured array
            try:
                self.names = lines.dtype.names
            except AttributeError:
                # Still no col names available so auto-generate them
                try:
                    first_data_vals = next(iter(lines))
                except StopIteration:
                    raise core.InconsistentTableError(
                        'No data lines found so cannot autogenerate column names')
                n_data_cols = len(first_data_vals)
                self.names = [self.auto_format % i for i in range(1, n_data_cols+1)]

        self._set_cols_from_names()

        # ``lines`` could be one of: NumPy recarray, DictLikeNumpy obj, Python
        # list of lists. If NumPy recarray then set col.type accordingly.  In
        # the other two cases convert the data values to strings so the usual
        # data converter processing will get the correct type.
        if core.has_numpy and isinstance(lines, numpy.ndarray):
            for col in self.cols:
                type_name = lines[col.name].dtype.name
                if 'int' in type_name:
                    col.type = core.IntType
                elif 'float' in type_name:
                    col.type = core.FloatType
                elif 'str' in type_name:
                    col.type = core.StrType
        else:
            # lines is a list of lists or DictLikeNumpy.  
            col_types = {}
            col_indexes = [col.index for col in self.cols]
            for vals in lines:
                for col_index in col_indexes:
                    val = vals[col_index]
                    col_type_set = col_types.setdefault(col_index, set())
                    col_type_set.add(get_val_type(val))
            for col in self.cols:
                col.type = get_lowest_type(col_types[col.index])
            

class MemorySplitter(core.BaseSplitter):
    """Splitter for data already in memory.  It is assumed that ``lines`` are
    iterable and that each line (aka row) is also an iterable object that
    provides the column values for that row."""
    def __call__(self, lines):
        for vals in lines:
            yield vals

class MemoryData(core.BaseData):
    """Memory table data reader.  Same as the BaseData reader but with a
    special splitter and a "pass-thru" process_lines function."""

    splitter_class = MemorySplitter

    def process_lines(self, lines):
        return lines

