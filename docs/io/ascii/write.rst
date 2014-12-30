.. include:: references.txt

.. _astropy.io.ascii_write:

Writing tables
--------------

:mod:`astropy.io.ascii` is able to write ASCII tables out to a file or file-like
object using the same class structure and basic user interface as for reading
tables.

The |write| function provides a way to write a data table as a
formatted ASCII table.  For example::

  >>> import numpy as np
  >>> from astropy.io import ascii
  >>> x = np.array([1, 2, 3])
  >>> y = x ** 2
  >>> ascii.write([x, y], 'values.dat', names=['x', 'y'])

The ``values.dat`` file will then contain::

  x y
  1 1
  2 4
  3 9

Most of the input table :ref:`supported_formats` for
reading are also available for writing.  This provides a great deal of
flexibility in the format for writing.  The example below writes the data as a
LaTeX table, using the option to send the output to ``sys.stdout`` instead of a
file::

  >>> ascii.write(data, format='latex')  # doctest: +SKIP
  \begin{table}
  \begin{tabular}{cc}
  x & y \\
  1 & 1 \\
  2 & 4 \\
  3 & 9 \\
  \end{tabular}
  \end{table}

There is also a faster Cython engine for writing simple formats,
which is enabled by default for these formats (see :ref:`fast_ascii_io`).
To disable this engine, use the parameter ``fast_writer``::

   >>> ascii.write(data, 'values.csv', format='csv', fast_writer=False)  # doctest: +SKIP

Input data format
^^^^^^^^^^^^^^^^^

The input ``table`` argument to |write| can be any value that is supported for
initializing a |Table| object.  This is documented in detail in the
:ref:`construct_table` section and includes creating a table with a list of
columns, a dictionary of columns, or from `numpy` arrays (either structured or
homogeneous).  The sections below show a few examples.

Table or NumPy structured array
"""""""""""""""""""""""""""""""

An AstroPy |Table| object or a NumPy `structured array`_ (or record array) can
serve as input to the |write| function.

::

    >>> from astropy.io import ascii
    >>> from astropy.table import Table

    >>> data = Table({'a': [1, 2, 3],
    ...               'b': [4.0, 5.0, 6.0]},
    ...              names=['a', 'b'])
    >>> ascii.write(data)
    a b
    1 4.0
    2 5.0
    3 6.0

    >>> data = np.array([(1, 2., 'Hello'), (2, 3., "World")],
    ...                 dtype=('i4,f4,a10'))
    >>> ascii.write(data)
    f0 f1 f2
    1 2.0 Hello
    2 3.0 World

The output of :mod:`astropy.io.ascii.read` is a |Table| or NumPy array data
object that can be an input to the |write| function.

::

    >>> data = ascii.read('t/daophot.dat', format='daophot')  # doctest: +SKIP
    >>> ascii.write(data, 'space_delimited_table.dat')  # doctest: +SKIP

List of lists
"""""""""""""

A list of Python lists (or any iterable object) can be used as input::

    >>> x = [1, 2, 3]
    >>> y = [4, 5.2, 6.1]
    >>> z = ['hello', 'world', '!!!']
    >>> data = [x, y, z]

    >>> ascii.write(data)
    col0 col1 col2
    1 4.0 hello
    2 5.2 world
    3 6.1 !!!

The ``data`` object does not contain information about the column names so
|Table| has chosen them automatically.  To specify the names, provide the
``names`` keyword argument.  This example also shows excluding one of the columns
from the output::

    >>> ascii.write(data, names=['x', 'y', 'z'], exclude_names=['y'])
    x z
    1 hello
    2 world
    3 !!!


Dict of lists
"""""""""""""

A dictionary containing iterable objects can serve as input to |write|.  Each
dict key is taken as the column name while the value must be an iterable object
containing the corresponding column values.

Since a Python dictionary is not ordered the output column order will be
unpredictable unless the ``names`` argument is provided.

::

    >>> data = {'x': [1, 2, 3],
    ...         'y': [4, 5.2, 6.1],
    ...         'z': ['hello', 'world', '!!!']}
    >>> ascii.write(data, names=['x', 'y', 'z'])
    x y z
    1 4.0 hello
    2 5.2 world
    3 6.1 !!!


.. _io_ascii_write_parameters:

Parameters for ``write()``
^^^^^^^^^^^^^^^^^^^^^^^^^^

The |write| function accepts a number of parameters that specify the detailed output table
format.  Each of the :ref:`supported_formats` is handled by a corresponding Writer class that
can define different defaults, so the descriptions below sometimes mention "typical"
default values.  This refers to the :class:`~astropy.io.ascii.Basic` writer and other
similar Writer classes.

Some output format Writer classes, e.g. :class:`~astropy.io.ascii.Latex` or
:class:`~astropy.io.ascii.AASTex` accept additional keywords, that can
customize the output further. See the documentation of these classes for
details.

**output** : output specifier
  There are two ways to specify the output for the write operation:

  - Name of a file (string)
  - File-like object (from open(), StringIO, etc)

**table** : input table
  Any value that is supported for initializing a |Table| object (see :ref:`construct_table`).

**format** : output format (default='basic')
  This specifies the format of the ASCII table to be written, for
  example if it is a basic character delimited table, fixed format table, or a
  CDS-compatible table, etc.  The value of this parameter must
  be one of the :ref:`supported_formats`.

**delimiter** : column delimiter string
  A one-character string used to separate fields which typically defaults to the space character.
  Other common values might be "," or "|" or "\\t".

**comment** : string defining start of a comment line in output table
  For the :class:`~astropy.io.ascii.Basic` Writer this defaults to "# ".
  Which and how comments are written depends on the format chosen.
  The comments are defined as a list of strings in the input table
  `meta['comments']` element. Comments in the metadata of the given
  |Table| will normally be written before the header, although
  :class:`~astropy.io.ascii.CommentedHeader` writes table comments after the
  commented header. To disable writing comments, set ``comment=False``.

**formats**: dict of data type converters
  For each key (column name) use the given value to convert the column data to a string.
  If the format value is string-like then it is used as a Python format statement,
  e.g. '%0.2f' % value.  If it is a callable function then that function
  is called with a single argument containing the column value to be converted.
  Example::

    astropy.io.ascii.write(table, sys.stdout, formats={'XCENTER': '%12.1f',
                                                 'YCENTER': lambda x: round(x, 1)},

**names**: list of names corresponding to each data column
  Define the complete list of names for each data column.  This will override
  names determined from the data table (if available).  If not supplied then
  use names from the data table or auto-generated names.

**include_names**: list of names to include in output
  From the list of column names found from the data table or the ``names``
  parameter, select for output only columns within this list.  If not supplied
  then include all names.

**exclude_names**: list of names to exclude from output
  Exclude these names from the list of output columns.  This is applied *after*
  the ``include_names`` filtering.  If not specified then no columns are excluded.

**fill_values**: fill value specifier of lists
  This can be used to fill missing values in the table or replace values with special meaning.

  See the :ref:`replace_bad_or_missing_values` section for more information on the syntax.
  The syntax is almost the same as when reading a table.
  There is a special value ``astropy.io.ascii.masked`` that is used a say "output this string
  for all masked values in a masked table (the default is to use a ``'--'``)::

      >>> import sys
      >>> from astropy.table import Table, Column, MaskedColumn
      >>> from astropy.io import ascii
      >>> t = Table([(1, 2), (3, 4)], names=('a', 'b'), masked=True)
      >>> t['a'].mask = [True, False]
      >>> ascii.write(t, sys.stdout)
      a b
      -- 3
      2 4
      >>> ascii.write(t, sys.stdout, fill_values=[(ascii.masked, 'N/A')])
      a b
      N/A 3
      2 4

  If no ``fill_values`` is applied for masked values in ``astropy.io.ascii``, the default set
  with ``numpy.ma.masked_print_option.set_display`` applies (usually that is also ``'--'``)::

      >>> ascii.write(t, sys.stdout, fill_values=[])
      a b
      -- 3
      2 4

  Note that when writing a table all values are converted to strings, before
  any value is replaced. Because ``fill_values`` only replaces cells that
  are an exact match to the specification, you need to provide the string
  representation (stripped of whitespace) for each value. For example, in
  the following commands ``-99`` is formatted with two digits after the
  comma, so we need to replace ``-99.00`` and not ``-99``::

      >>> t = Table([(-99, 2), (3, 4)], names=('a', 'b'))
      >>> ascii.write(t, sys.stdout, fill_values = [('-99.00', 'no data')],
      ...             formats={'a': '%4.2f'})
      a b
      "no data" 3
      2.00 4

  Similarly, if you replace a value in a column that has a fixed length format,
  e.g. ``'f4.2'``, then the string you want to replace must have the same
  number of characters, in the example above ``fill_values=[(' nan',' N/A')]``
  would work.

**fill_include_names**: list of column names, which are affected by ``fill_values``.
  If not supplied, then ``fill_values`` can affect all columns.

**fill_exclude_names**: list of column names, which are not affected by ``fill_values``.
  If not supplied, then ``fill_values`` can affect all columns.

**fast_writer**: whether to use the fast Cython writer
  If this parameter is ``None`` (which it is by default), |write| will attempt
  to use the faster writer (described in :ref:`fast_ascii_io`) if possible.
  Specifying ``fast_writer=False`` disables this behavior.

**Writer** : Writer class (*deprecated* in favor of ``format``)
  This specifies the top-level format of the ASCII table to be written, for
  example if it is a basic character delimited table, fixed format table, or a
  CDS-compatible table, etc.  The value of this parameter must be a Writer
  class.  For basic usage this means one of the built-in :ref:`extension_reader_classes`.
  Note: Reader classes and Writer classes are synonymous, in other
  words Reader classes can also write, but for historical reasons they are
  often called Reader classes.
