.. include:: references.txt

.. _astropy.io.ascii_write:

Writing Tables
**************

:mod:`astropy.io.ascii` is able to write ASCII tables out to a file or file-like
object using the same class structure and basic user interface as for reading
tables.

The |write| function provides a way to write a data table as a
formatted ASCII table.

Examples
--------

..
  EXAMPLE START
  Writing ASCII Tables Using astropy.io.ascii

To write a formatted ASCII table using the |write| function::

  >>> import numpy as np
  >>> from astropy.io import ascii
  >>> x = np.array([1, 2, 3])
  >>> y = x ** 2
  >>> ascii.write([x, y], 'values.dat', names=['x', 'y'], overwrite=True)

The ``values.dat`` file will then contain::

  x y
  1 1
  2 4
  3 9

Most of the input table :ref:`supported_formats` for
reading are also available for writing. This provides a great deal of
flexibility in the format for writing. The example below writes the data as a
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

..
  EXAMPLE END

Input Data Format
=================

The input ``table`` argument to |write| can be any value that is supported for
initializing a |Table| object. This is documented in detail in the
:ref:`construct_table` section and includes creating a table with a list of
columns, a dictionary of columns, or from `numpy` arrays (either structured or
homogeneous).


Table or NumPy Structured Array
-------------------------------

An Astropy |Table| object or a NumPy `structured array`_ (or record array) can
serve as input to the |write| function.

Example
=======

..
  EXAMPLE START
  Creating a Table with a NumPy Structured Array or an Existing Table

To create a table with a ``numpy`` structured array or an existing table::

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

    >>> data = ascii.read('data/daophot.dat', format='daophot')  # doctest: +SKIP
    >>> ascii.write(data, 'space_delimited_table.dat')  # doctest: +SKIP

..
  EXAMPLE END


List of ``list`` Objects
------------------------

A list of Python ``list`` objects (or any iterable object) can be used as
input::

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
|Table| has chosen them automatically. To specify the names, provide the
``names`` keyword argument. This example also shows excluding one of the columns
from the output::

    >>> ascii.write(data, names=['x', 'y', 'z'], exclude_names=['y'])
    x z
    1 hello
    2 world
    3 !!!


``dict`` of ``list`` Objects
----------------------------

A dictionary containing iterable objects can serve as input to |write|. Each
dict key is taken as the column name while the value must be an iterable object
containing the corresponding column values.

Since a Python dictionary is not ordered, the output column order will be
unpredictable unless the ``names`` argument is provided.

Example
=======

..
  EXAMPLE START
  Writing a Table from a dict of list Objects

To write a table from a ``dict`` of ``list`` objects::

    >>> data = {'x': [1, 2, 3],
    ...         'y': [4, 5.2, 6.1],
    ...         'z': ['hello', 'world', '!!!']}
    >>> ascii.write(data, names=['x', 'y', 'z'])
    x y z
    1 4.0 hello
    2 5.2 world
    3 6.1 !!!

..
  EXAMPLE END

.. _io_ascii_write_parameters:


Parameters for ``write()``
--------------------------

The |write| function accepts a number of parameters that specify the detailed
output table format. Each of the :ref:`supported_formats` is handled by a
corresponding Writer class that can define different defaults, so the
descriptions below sometimes mention "typical" default values. This refers to
the :class:`~astropy.io.ascii.Basic` writer and other similar Writer classes.

Some output format Writer classes (e.g., :class:`~astropy.io.ascii.Latex` or
:class:`~astropy.io.ascii.AASTex`) accept additional keywords that can
customize the output further. See the documentation of these classes for
details.

**output**: output specifier
  There are two ways to specify the output for the write operation:

  - Name of a file (string)
  - File-like object (from open(), StringIO, etc.)

**table**: input table
  Any value that is supported for initializing a |Table| object (see :ref:`construct_table`).

**format**: output format (default='basic')
  This specifies the format of the ASCII table to be written, such as a basic
  character delimited table, fixed-format table, or a CDS-compatible table,
  etc. The value of this parameter must be one of the :ref:`supported_formats`.

**delimiter**: column delimiter string
  A one-character string used to separate fields which typically defaults to
  the space character. Other common values might be "," or "|" or "\\t".

**comment**: string defining start of a comment line in output table
  For the :class:`~astropy.io.ascii.Basic` Writer this defaults to "# ".
  Which comments are written and how depends on the format chosen.
  The comments are defined as a list of strings in the input table
  ``meta['comments']`` element. Comments in the metadata of the given
  |Table| will normally be written before the header, although
  :class:`~astropy.io.ascii.CommentedHeader` writes table comments after the
  commented header. To disable writing comments, set ``comment=False``.

**formats**: dict of data type converters
  For each key (column name) use the given value to convert the column data to
  a string. If the format value is string-like, then it is used as a Python
  format statement (e.g., '%0.2f' % value). If it is a callable function, then
  that function is called with a single argument containing the column value to
  be converted. Example::

    astropy.io.ascii.write(table, sys.stdout, formats={'XCENTER': '%12.1f',
                                                 'YCENTER': lambda x: round(x, 1)},

**names**: list of names corresponding to each data column
  Define the complete list of names for each data column. This will override
  names determined from the data table (if available). If not supplied then
  use names from the data table or auto-generated names.

**include_names**: list of names to include in output
  From the list of column names found from the data table or the ``names``
  parameter, select for output only columns within this list. If not supplied
  then include all names.

**exclude_names**: list of names to exclude from output
  Exclude these names from the list of output columns. This is applied *after*
  the ``include_names`` filtering. If not specified then no columns are excluded.

**fill_values**: list of fill value specifiers
  This can be used to fill missing values in the table or replace values with special meaning.

  See the :ref:`replace_bad_or_missing_values` section for more information on
  the syntax. The syntax is almost the same as when reading a table.
  There is a special value ``astropy.io.ascii.masked`` that is used to say
  "output this string for all masked values in a masked table" (the default is
  to use an empty string ``""``)::

      >>> import sys
      >>> from astropy.table import Table, Column, MaskedColumn
      >>> from astropy.io import ascii
      >>> t = Table([(1, 2), (3, 4)], names=('a', 'b'), masked=True)
      >>> t['a'].mask = [True, False]
      >>> ascii.write(t, sys.stdout)
      a b
      "" 3
      2 4
      >>> ascii.write(t, sys.stdout, fill_values=[(ascii.masked, 'N/A')])
      a b
      N/A 3
      2 4

  Note that when writing a table, all values are converted to strings before
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

  Similarly, if you replace a value in a column that has a fixed length format
  (e.g., ``'f4.2'``), then the string you want to replace must have the same
  number of characters. In the example above, ``fill_values=[(' nan',' N/A')]``
  would work.

**fill_include_names**: list of column names, which are affected by ``fill_values``
  If not supplied, then ``fill_values`` can affect all columns.

**fill_exclude_names**: list of column names, which are not affected by ``fill_values``
  If not supplied, then ``fill_values`` can affect all columns.

**fast_writer**: whether to use the fast Cython writer
  If this parameter is ``None`` (which it is by default), |write| will attempt
  to use the faster writer (described in :ref:`fast_ascii_io`) if possible.
  Specifying ``fast_writer=False`` disables this behavior.

**Writer** : Writer class (*deprecated* in favor of ``format``)
  This specifies the top-level format of the ASCII table to be written, such as
  a basic character delimited table, fixed-format table, or a CDS-compatible
  table, etc. The value of this parameter must be a Writer class. For basic
  usage this means one of the built-in :ref:`extension_reader_classes`.
  Note that Reader classes and Writer classes are synonymous; in other
  words, Reader classes can also write, but for historical reasons they are
  often called Reader classes.


.. _ecsv_format:

ECSV Format
-----------

The `Enhanced Character-Separated Values (ECSV) format
<https://github.com/astropy/astropy-APEs/blob/master/APE6.rst>`_ can be used to
write ``astropy`` `~astropy.table.Table` or `~astropy.table.QTable` datasets to
a text-only data file and then read the table back without loss of information.
The format handles the key issue of serializing column specifications and table
metadata by using a YAML-encoded data structure. The actual tabular data are
stored in a standard character separated values (CSV) format, giving
compatibility with a wide variety of non-specialized CSV table readers.

.. _ecsv_format_mixin_columns:

Mixin Columns
-------------

Starting with ``astropy`` 2.0 it is possible to store not only standard
`~astropy.table.Column` objects to ECSV but also the following
:ref:`mixin_columns`:

- `astropy.time.Time`
- `astropy.time.TimeDelta`
- `astropy.units.Quantity`
- `astropy.coordinates.Latitude`
- `astropy.coordinates.Longitude`
- `astropy.coordinates.Angle`
- `astropy.coordinates.Distance`
- `astropy.coordinates.EarthLocation`
- `astropy.coordinates.SkyCoord`

In general, a mixin column may contain multiple data components as well as
object attributes beyond the standard `~astropy.table.Column` attributes like
``format`` or ``description``. Storing such mixin columns is done by replacing
the mixin column with column(s) representing the underlying data component(s)
and then inserting metadata which informs the reader of how to reconstruct the
original column. For example, a `~astropy.coordinates.SkyCoord` mixin column in
``'spherical'`` representation would have data attributes ``ra``, ``dec``,
``distance``, along with object attributes like ``representation_type`` or
``frame``.


Example
=======

..
  EXAMPLE START
  Writing a Table with a SkyCoord Column in ECSV Format

Creating a table with a `~astropy.coordinates.SkyCoord` column can be accomplished with a mixin
column, which is supported by `ECSV <https://docs.astropy.org/en/stable/api/astropy.io.ascii.Ecsv.html>`_. To store a mixin column::

  >>> from astropy.io import ascii
  >>> from astropy.coordinates import SkyCoord
  >>> import astropy.units as u
  >>> from astropy.time import Time
  >>> from astropy.table import QTable, Column

  >>> sc = SkyCoord(ra=[1,2]*u.deg, dec=[3,4]*u.deg, distance=[5,6]*u.m,
  ...               frame='fk4', obstime=Time('2000:001'))
  >>> sc.info.description = 'flying circus'
  >>> c = Column([1,2])
  >>> q = [1,2]*u.m
  >>> q.info.format = '.2f'
  >>> t = QTable([c, q, sc], names=['c', 'q', 'sc'])

  >>> ascii.write(t, format='ecsv')   # doctest: +SKIP
  # %ECSV 0.9
  # ---
  # datatype:
  # - {name: c, datatype: int64}
  # - {name: q, unit: m, datatype: float64}
  # - {name: sc.ra, unit: deg, datatype: float64}
  # - {name: sc.dec, unit: deg, datatype: float64}
  # - {name: sc.distance, unit: m, datatype: float64}
  # meta: !!omap
  # - __serialized_columns__:
  #     q:
  #       __class__: astropy.units.quantity.Quantity
  #       value: !astropy.table.SerializedColumn {name: q}
  #     sc:
  #       __class__: astropy.coordinates.sky_coordinate.SkyCoord
  #       __info__: {description: flying circus}
  #       dec: !astropy.table.SerializedColumn
  #         __class__: astropy.coordinates.angles.Latitude
  #         value: !astropy.table.SerializedColumn {name: sc.dec}
  #       distance: !astropy.table.SerializedColumn
  #         __class__: astropy.coordinates.distances.Distance
  #         value: !astropy.table.SerializedColumn {name: sc.distance}
  #       equinox: !astropy.time.Time {format: byear_str, in_subfmt: '*', jd1: 2400000.5,
  #         jd2: 33281.92345905, out_subfmt: '*', precision: 3, scale: tai}
  #       frame: fk4
  #       obstime: !astropy.time.Time {format: yday, in_subfmt: '*', jd1: 2451544.5, jd2: 0.0,
  #         out_subfmt: '*', precision: 3, scale: utc}
  #       ra: !astropy.table.SerializedColumn
  #         __class__: astropy.coordinates.angles.Longitude
  #         value: !astropy.table.SerializedColumn {name: sc.ra}
  #         wrap_angle: !astropy.coordinates.Angle
  #           unit: !astropy.units.Unit {unit: deg}
  #           value: 360.0
  #       representation: spherical
  # schema: astropy-2.0
  c q sc.ra sc.dec sc.distance
  1 1.0 1.0 3.0 5.0
  2 2.0 2.0 4.0 6.0

The ``'__class__'`` keyword gives the fully-qualified class name and must be
one of the specifically allowed ``astropy`` classes. There is no option to add
user-specified allowed classes. The ``'__info__'`` keyword contains values for
standard `~astropy.table.Column` attributes like ``description`` or ``format``,
for any mixin columns that are represented by more than one serialized column.

..
  EXAMPLE END

.. _ecsv_format_masked_columns:

Masked Columns
--------------

By default, the ECSV format uses an empty (zero-length) string in the output
table to represent masked or missing data in `~astropy.table.MaskedColumn`
columns. In certain cases this may not be sufficient:

- String column that contains empty (zero-length) string(s) as valid data.
- Masked data values must be stored so those values can later be unmasked.

In this case, there is an available mechanism to specify that the full data
and the mask itself should be written as columns in the output table as
shown in the example below. For further context see the section on
:ref:`table_serialization_methods`.

Example
=======

..
  EXAMPLE START
  Using ECSV Format to Write Astropy Tables with Masked or Missing Data

To specify that the full data and the mask itself should be written as columns
in the output table::

  >>> from astropy.table.table_helpers import simple_table
  >>> t = simple_table(masked=True)
  >>> t['c'][0] = ""  # Valid empty string in data
  >>> t
  <Table masked=True length=3>
    a      b     c
  int64 float64 str1
  ----- ------- ----
     --     1.0
      2     2.0   --
      3      --    e

Now we tell ECSV writer to output separate data and mask columns for the
string column ``'c'``::

  >>> t['c'].info.serialize_method['ecsv'] = 'data_mask'

When this is written out, notice that the output shows all of the
data values for the ``'c'`` column (including the masked ``'d'``
value) and a new column ``'c.masked'``. It also stores metadata
that tells the ECSV reader to interpret the ``'c'`` and ``'c.masked'``
columns as components of one `~astropy.table.MaskedColumn` object:

.. doctest-skip::

  >>> ascii.write(t, format='ecsv')
  # %ECSV 0.9
  # ---
  # datatype:
  # - {name: a, datatype: int64}
  # - {name: b, datatype: float64}
  # - {name: c, datatype: string}
  # - {name: c.mask, datatype: bool}
  # meta: !!omap
  # - __serialized_columns__:
  #     c:
  #       __class__: astropy.table.column.MaskedColumn
  #       data: !astropy.table.SerializedColumn {name: c}
  #       mask: !astropy.table.SerializedColumn {name: c.mask}
  # schema: astropy-2.0
  a b c c.mask
  "" 1.0 "" False
  2 2.0 d True
  3 "" e False

When you read this back in, the empty (zero-length) string in the
first row of column ``'c'`` will be preserved. You can also write
all of the columns out as data and mask pairs using the Unified I/O
interface for tables with the ``serialize_method`` keyword argument::

  >>> t.write('out.ecsv', format='ascii.ecsv', serialize_method='data_mask')  # doctest: +SKIP

In this case, all data values (including those "under the mask" in the
original table) will be restored exactly when you read the file back.

..
  EXAMPLE END
