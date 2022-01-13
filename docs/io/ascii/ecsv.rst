.. _ecsv_format:

ECSV Format
===========

The `Enhanced Character-Separated Values (ECSV) format
<https://github.com/astropy/astropy-APEs/blob/main/APE6.rst>`_ can be used to
write ``astropy`` `~astropy.table.Table` or `~astropy.table.QTable` datasets to
a text-only human readable data file and then read the table back without loss
of information. The format stores column specifications like unit and data type
along with table metadata by using a YAML header data structure. The
actual tabular data are stored in a standard character separated values (CSV)
format, giving compatibility with a wide variety of non-specialized CSV table
readers.

.. attention::

    The ECSV format is the recommended way to store Table data in a
    human-readable ASCII file. This includes use cases from informal
    use in science research to production pipelines and data systems.

    In addition to Python, ECSV is supported in `TOPCAT
    <http://www.star.bris.ac.uk/~mbt/topcat/>`_ and in the java `STIL
    <http://www.star.bris.ac.uk/~mbt/topcat/sun253/inEcsv.html>`_ library. .

Usage
-----

When writing in the ECSV format there are only two choices for the delimiter,
either space or comma, with space being the default. Any other value of
``delimiter`` will give an error. For reading the delimiter is specified within
the file itself.

Apart from the delimiter, the only other applicable read/write arguments are
``names``, ``include_names``, and ``exclude_names``. All other arguments will be
either ignored or raise an error.

Simple Table
------------
..
  EXAMPLE START
  Writing Data Tables as ECSV: Simple Table

The following writes a table as a simple space-delimited file. The
ECSV format is auto-selected due to ``.ecsv`` suffix::

  >>> import numpy as np
  >>> from astropy.table import Table
  >>> data = Table()
  >>> data['a'] = np.array([1, 2], dtype=np.int8)
  >>> data['b'] = np.array([1, 2], dtype=np.float32)
  >>> data['c'] = np.array(['hello', 'world'])
  >>> data.write('my_data.ecsv')  # doctest: +SKIP

The contents of ``my_data.ecsv`` are shown below::

  # %ECSV 1.0
  # ---
  # datatype:
  # - {name: a, datatype: int8}
  # - {name: b, datatype: float32}
  # - {name: c, datatype: string}
  # schema: astropy-2.0
  a b c
  1 1.0 hello
  2 2.0 world

The ECSV header is the section prefixed by the ``#`` comment character. An ECSV
file must start with the ``%ECSV <version>`` line. The ``datatype`` element
defines the list of columns and the ``schema`` relates to astropy-specific
extensions that are used for writing `Mixin Columns`_.

..
  EXAMPLE END

Masked Data
-----------

You can write masked (or "missing") data in the ECSV format in two different
ways, either using an empty string to represent missing values or by splitting
the masked columns into separate data and mask columns.

Empty String
""""""""""""

The first (default) way uses an empty string as a marker in place of
masked values. This is a bit more common outside of ``astropy`` and does not
require any astropy-specific extensions.

  >>> from astropy.table import MaskedColumn
  >>> t = Table()
  >>> t['x'] = MaskedColumn([1.0, 2.0, 3.0], unit='m', dtype='float32')
  >>> t['x'][1] = np.ma.masked
  >>> t['y'] = MaskedColumn([False, True, False], dtype='bool')
  >>> t['y'][0] = np.ma.masked

  >>> t.write('my_data.ecsv', format='ascii.ecsv', overwrite=True)  # doctest: +SKIP

The contents of ``my_data.ecsv`` are shown below::

  # %ECSV 1.0
  # ---
  # datatype:
  # - {name: x, unit: m, datatype: float32}
  # - {name: y, datatype: bool}
  # schema: astropy-2.0
  x y
  1.0 ""
  "" True
  3.0 False

To read this back, you would run the following::

  >>> Table.read('my_data.ecsv')  # doctest: +SKIP
  <Table length=3>
     x      y
     m
  float32  bool
  ------- -----
      1.0    --
       --  True
      3.0 False

Data + Mask
"""""""""""

The second way is to tell the writer to break any masked column into a data
column and a mask column by supplying the ``serialize_method='data_mask'``
argument::

  >>> t.write('my_data.ecsv', serialize_method='data_mask', overwrite=True)  # doctest: +SKIP

There are two main reasons you might want to do this:

- Storing the data "under the mask" instead of replacing it with an empty string.
- Writing a string column that contains empty strings which are not masked.

The contents of ``my_data.ecsv`` are shown below. First notice that there are
two new columns ``x.mask`` and ``y.mask`` that have been added, and these explicitly
record the mask values for those columns. Next notice now that the ECSV
header is a bit more complex and includes the astropy-specific extensions that
tell the reader how to interpret the plain CSV columns ``x, x.mask, y, y.mask``
and reassemble them back into the appropriate masked columns.
::

  # %ECSV 1.0
  # ---
  # datatype:
  # - {name: x, unit: m, datatype: float32}
  # - {name: x.mask, datatype: bool}
  # - {name: y, datatype: bool}
  # - {name: y.mask, datatype: bool}
  # meta: !!omap
  # - __serialized_columns__:
  #     x:
  #       __class__: astropy.table.column.MaskedColumn
  #       data: !astropy.table.SerializedColumn {name: x}
  #       mask: !astropy.table.SerializedColumn {name: x.mask}
  #     y:
  #       __class__: astropy.table.column.MaskedColumn
  #       data: !astropy.table.SerializedColumn {name: y}
  #       mask: !astropy.table.SerializedColumn {name: y.mask}
  # schema: astropy-2.0
  x x.mask y y.mask
  1.0 False False True
  2.0 True True False
  3.0 False False False

.. note::

   For the security minded, the ``__class__`` value must within an allowed list
   of astropy classes that are trusted by the reader. You cannot use an
   arbitrary class here.

..
  EXAMPLE START
  Using ECSV Format to Write Astropy Tables with Masked or Missing Data

Per-column control
@@@@@@@@@@@@@@@@@@

In rare cases it may be necessary to specify the serialization method for each
column individually. This is shown in the example below::

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
string column ``'c'``:

.. doctest-skip::

  >>> t['c'].info.serialize_method['ecsv'] = 'data_mask'
  >>> ascii.write(t, format='ecsv')
  # %ECSV 1.0
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

When you read this back in, both the empty (zero-length) string and the masked
``'d'`` value in the column ``'c'`` will be preserved.

..
  EXAMPLE END

.. _ecsv_format_mixin_columns:

Mixin Columns
-------------

It is possible to store not only standard `~astropy.table.Column` and
`~astropy.table.MaskedColumn` objects to ECSV but also the following
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
- `astropy.table.NdarrayMixin`
- Coordinate representation types such as `astropy.coordinates.SphericalRepresentation`

In general, a mixin column may contain multiple data components as well as
object attributes beyond the standard `~astropy.table.Column` attributes like
``format`` or ``description``. Storing such mixin columns is done by replacing
the mixin column with column(s) representing the underlying data component(s)
and then inserting metadata which informs the reader of how to reconstruct the
original column. For example, a `~astropy.coordinates.SkyCoord` mixin column in
``'spherical'`` representation would have data attributes ``ra``, ``dec``,
``distance``, along with object attributes like ``representation_type`` or
``frame``.

..
  EXAMPLE START
  Writing a Table with a SkyCoord Column in ECSV Format

This example demonstrates writing a `~astropy.table.QTable` that has `~astropy.time.Time`
and `~astropy.coordinates.SkyCoord` mixin columns::

  >>> from astropy.coordinates import SkyCoord
  >>> import astropy.units as u
  >>> from astropy.table import QTable

  >>> sc = SkyCoord(ra=[1, 2] * u.deg, dec=[3, 4] * u.deg)
  >>> sc.info.description = 'flying circus'
  >>> q = [1, 2] * u.m
  >>> q.info.format = '.2f'
  >>> t = QTable()
  >>> t['c'] = [1, 2]
  >>> t['q'] = q
  >>> t['sc'] = sc

  >>> t.write('my_data.ecsv')  # doctest: +SKIP

The contents of ``my_data.ecsv`` are below::

  # %ECSV 1.0
  # ---
  # datatype:
  # - {name: c, datatype: int64}
  # - {name: q, unit: m, datatype: float64, format: .2f}
  # - {name: sc.ra, unit: deg, datatype: float64}
  # - {name: sc.dec, unit: deg, datatype: float64}
  # meta: !!omap
  # - __serialized_columns__:
  #     q:
  #       __class__: astropy.units.quantity.Quantity
  #       __info__: {format: .2f}
  #       unit: !astropy.units.Unit {unit: m}
  #       value: !astropy.table.SerializedColumn {name: q}
  #     sc:
  #       __class__: astropy.coordinates.sky_coordinate.SkyCoord
  #       __info__: {description: flying circus}
  #       dec: !astropy.table.SerializedColumn
  #         __class__: astropy.coordinates.angles.Latitude
  #         unit: &id001 !astropy.units.Unit {unit: deg}
  #         value: !astropy.table.SerializedColumn {name: sc.dec}
  #       frame: icrs
  #       ra: !astropy.table.SerializedColumn
  #         __class__: astropy.coordinates.angles.Longitude
  #         unit: *id001
  #         value: !astropy.table.SerializedColumn {name: sc.ra}
  #         wrap_angle: !astropy.coordinates.Angle
  #           unit: *id001
  #           value: 360.0
  #       representation_type: spherical
  # schema: astropy-2.0
  c q sc.ra sc.dec
  1 1.0 1.0 3.0
  2 2.0 2.0 4.0

The ``'__class__'`` keyword gives the fully-qualified class name and must be
one of the specifically allowed ``astropy`` classes. There is no option to add
user-specified allowed classes. The ``'__info__'`` keyword contains values for
standard `~astropy.table.Column` attributes like ``description`` or ``format``,
for any mixin columns that are represented by more than one serialized column.

..
  EXAMPLE END

.. _ecsv_format_masked_columns:

Multidimensional Columns
------------------------

Using ECSV it is possible to write a table that contains multidimensional
columns (both masked and unmasked). This is done by encoding each element as a
string using JSON. This functionality works for all column types that are
supported by ECSV including :ref:`mixin_columns`. This capability is added in
astropy 4.3 and ECSV version 1.0.

..
  EXAMPLE START
  Using ECSV Format to Write Astropy Tables with Multidimensional Columns

We start by defining a table with 2 rows where each element in the second column
``'b'`` is itself a 3x2 array::

  >>> t = Table()
  >>> t['a'] = ['x', 'y']
  >>> t['b'] = np.arange(12, dtype=np.float64).reshape(2, 3, 2)
  >>> t
  <Table length=2>
   a        b
  str1 float64[3,2]
  ---- ------------
     x   0.0 .. 5.0
     y  6.0 .. 11.0

  >>> t['b'][0]
  array([[0., 1.],
        [2., 3.],
        [4., 5.]])

Now we can write this to ECSV and observe how the N-d column ``'b'`` has been
written as a string with ``datatype: string``. Notice also that the column
descriptor for the column includes the new ``subtype: float64[3,2]`` attribute
specifying the type and shape of each item.

.. doctest-skip::

  >>> ascii.write(t, format='ecsv')  # doctest: +SKIP
  # %ECSV 1.0
  # ---
  # datatype:
  # - {name: a, datatype: string}
  # - {name: b, datatype: string, subtype: 'float64[3,2]'}
  # schema: astropy-2.0
  a b
  x [[0.0,1.0],[2.0,3.0],[4.0,5.0]]
  y [[6.0,7.0],[8.0,9.0],[10.0,11.0]]

When you read this back in, the sequence of JSON-encoded column items are then
decoded using JSON back into the original N-d column.

..
  EXAMPLE END

Variable-length arrays
----------------------

ECSV supports storing multidimensional columns is when the length of each array
element may vary. This data structure is supported in the `FITS standard
<https://fits.gsfc.nasa.gov/fits_standard.html>`_. While ``numpy`` does not
natively support variable-length arrays, it is possible to represent such a
structure using an object-type array of typed ``np.ndarray`` objects. This is how
the ``astropy`` FITS reader outputs a variable-length array.

This capability is added in astropy 4.3 and ECSV version 1.0.

Most commonly variable-length arrays have a 1-d array in each cell of the
column. You might a column with 1-d ``np.ndarray`` cells having lengths of 2, 5,
and 3 respectively.

The ECSV standard and ``astropy`` also supports arbitrary N-d arrays in each
cell, where all dimensions except the last one must match. For instance you
could have a column with ``np.ndarray`` cells having shapes of ``(4,4,2)``,
``(4,4,5)``, and ``(4,4,3)`` respectively.

..
  EXAMPLE START
  Using ECSV Format to Write Astropy Tables with Variable-Length Arrays

The example below shows writing a variable-length 1-d array to ECSV. Notice the
new ECSV column attribute ``subtype: 'int64[null]'``. The ``[null]`` indicates a
variable length for the one dimension. If we had been writing the N-d example
above the subtype would have been ``int64[4,4,null]``.

.. doctest-skip::

  >>> t = Table()
  >>> t['a'] = np.empty(3, dtype=object)
  >>> t['a'] = [np.array([1, 2], dtype=np.int64),
  ...           np.array([3, 4, 5], dtype=np.int64),
  ...           np.array([6, 7, 8, 9], dtype=np.int64)]
  >>> ascii.write(t, format='ecsv')
  # %ECSV 1.0
  # ---
  # datatype:
  # - {name: a, datatype: string, subtype: 'int64[null]'}
  # schema: astropy-2.0
  a
  [1,2]
  [3,4,5]
  [6,7,8,9]

..
  EXAMPLE END

Object arrays
-------------

ECSV can store object-type columns with simple Python objects consisting of
``dict``, ``list``, ``str``, ``int``, ``float``, ``bool`` and ``None`` elements.
More precisely, any object that can be serialized to `JSON
<https://www.json.org/>`__ using the standard library `json
<https://docs.python.org/3/library/json.html>`__ package is supported.

..
  EXAMPLE START
  Using ECSV Format to Write Astropy Tables with Object Arrays

The example below shows writing an object array to ECSV. Because JSON requires
a double-quote around strings, and because ECSV requires ``""`` to represent
a double-quote within a string, one tends to get double-double quotes in this
representation.

.. doctest-skip::

  >>> t = Table()
  >>> t['a'] = np.array([{'a': 1},
  ...                    {'b': [2.5, None]},
  ...                    True], dtype=object)
  >>> ascii.write(t, format='ecsv')
  # %ECSV 1.0
  # ---
  # datatype:
  # - {name: a, datatype: string, subtype: json}
  # schema: astropy-2.0
  a
  "{""a"":1}"
  "{""b"":[2.5,null]}"
  true

..
  EXAMPLE END
