.. include:: references.txt
.. |join| replace:: :func:`~astropy.table.join`
.. |Quantity| replace:: :class:`~astropy.units.Quantity`
.. |Time| replace:: :class:`~astropy.time.Time`
.. |SkyCoord| replace:: :class:`~astropy.coordinates.SkyCoord`

.. _mixin_columns:

Mixin columns
---------------

Version 1.0 of astropy introduces a new concept of the "Mixin
Column" in tables which allows integration of appropriate non-|Column| based
class objects within a |Table| object.  These mixin column objects are not
converted in any way but are used natively.

The available built-in mixin column classes are:

- |Quantity|
- |SkyCoord|
- |Time|
- `~astropy.table.NdarrayMixin`

As a first example we can create a table and add a time column::

  >>> from astropy.table import Table
  >>> from astropy.time import Time
  >>> t = Table()
  >>> t['index'] = [1, 2]
  >>> t['time'] = Time(['2001-01-02T12:34:56', '2001-02-03T00:01:02'])
  >>> print(t)
  index           time
  ----- -----------------------
      1 2001-01-02T12:34:56.000
      2 2001-02-03T00:01:02.000

The important point here is that the ``time`` column is a bona fide |Time| object::

  >>> t['time']
  <Time object: scale='utc' format='isot' value=['2001-01-02T12:34:56.000' '2001-02-03T00:01:02.000']>
  >>> t['time'].mjd
  array([ 51911.52425926,  51943.00071759])

.. _quantity_and_qtable:

Quantity and QTable
^^^^^^^^^^^^^^^^^^^

The ability to natively handle |Quantity| objects within a table makes it
easier to manipulate tabular data with units in a natural and robust way.
However, this feature introduces an ambiguity because data with a unit
(e.g. from a FITS binary table) can be represented as either a |Column| with a
``unit`` attribute or as a |Quantity| object. In order to retain complete
backward compatibility with astropy versions prior to 1.0, a minor variant of
the |Table| class called |QTable| is available.  |QTable| is exactly the same
as |Table| except that |Quantity| is the default for any data column with a
defined unit.

If you take advantage of the |Quantity| infrastructure in your analysis then
|QTable| is the preferred way to create tables with units.  If instead you use
table column units more as a descriptive label then the plain |Table| class is
probably the best class to use.

To illustrate these concepts we first create a standard |Table| where we supply as input a
|Time| object and a |Quantity| object with units of ``m / s``.  In this case
the quantity is converted to a |Column| (which has a ``unit`` attribute but
does not have all the features of a |Quantity|)::

  >>> import astropy.units as u
  >>> t = Table()
  >>> t['index'] = [1, 2]
  >>> t['time'] = Time(['2001-01-02T12:34:56', '2001-02-03T00:01:02'])
  >>> t['velocity'] = [3, 4] * u.m / u.s

  >>> print(t)
  index           time          velocity
                                 m / s
  ----- ----------------------- --------
      1 2001-01-02T12:34:56.000      3.0
      2 2001-02-03T00:01:02.000      4.0

  >>> type(t['velocity'])
  <class 'astropy.table.column.Column'>

  >>> t['velocity'].unit
  Unit("m / s")

  >>> (t['velocity'] ** 2).unit  # WRONG because Column is not smart about unit
  Unit("m / s")

So instead let's do the same thing using a quantity table |QTable|::

  >>> from astropy.table import QTable

  >>> qt = QTable()
  >>> qt['index'] = [1, 2]
  >>> qt['time'] = Time(['2001-01-02T12:34:56', '2001-02-03T00:01:02'])
  >>> qt['velocity'] = [3, 4] * u.m / u.s

The ``velocity`` column is now a |Quantity| and behaves accordingly::

  >>> type(qt['velocity'])
  <class 'astropy.units.quantity.Quantity'>

  >>> qt['velocity'].unit
  Unit("m / s")

  >>> (qt['velocity'] ** 2).unit  # GOOD!
  Unit("m2 / s2")

You can easily convert |Table| to |QTable| and vice-versa::

  >>> qt2 = QTable(t)
  >>> type(qt2['velocity'])
  <class 'astropy.units.quantity.Quantity'>

  >>> t2 = Table(qt2)
  >>> type(t2['velocity'])
  <class 'astropy.table.column.Column'>

Mixin Attributes
^^^^^^^^^^^^^^^^

The usual column attributes ``name``, ``dtype``, ``unit``, ``format``, and
``description`` are available in any mixin column via the ``info`` property::

  >>> qt['velocity'].info.name
  'velocity'

This ``info`` property is a key bit of glue that allows for a
non-Column object to behave much like a column.

The same ``info`` property is also available in standard
`~astropy.table.Column` objects.  These ``info`` attributes like
``t['a'].info.name`` simply refer to the direct `~astropy.table.Column`
attribute (e.g. ``t['a'].name``) and can be used interchangeably.
Likewise in a `~astropy.units.Quantity` object, ``info.dtype``
attribute refers to the native ``dtype`` attribute of the object.

.. Note::

   When writing generalized code that handles column objects which
   might be mixin columns, one must *always* use the ``info``
   property to access column attributes.


.. _details_and_caveats:

Details and caveats
^^^^^^^^^^^^^^^^^^^

Most common table operations behave as expected when mixin columns are part of
the table.  However, there are limitations in the current implementation.

**Adding or inserting a row**

Adding or inserting a row works as expected only for mixin classes that are
mutable (data can changed internally) and that have an ``insert()`` method.
|Quantity| supports ``insert()`` but |Time| and |SkyCoord| do not.  If we try to
insert a row into the previously defined table an exception occurs::

  >>> qt.add_row((1, '2001-02-03T00:01:02', 5 * u.m / u.s))
  Traceback (most recent call last):
    ...
  ValueError: Unable to insert row because of exception in column 'time':
  'Time' object has no attribute 'insert'

**Initializing from a list of rows or a list of dicts**

This mode of initializing a table does not work with mixin columns, so both of
the following will fail::

   >>> qt = QTable([{'a': 1 * u.m, 'b': 2},
   ...              {'a': 2 * u.m, 'b': 3}])
   Traceback (most recent call last):
    ...
   ValueError: setting an array element with a sequence.

   >>> qt = QTable(rows=[[1 * u.m, 2],
   ...                   [2 * u.m, 3]])
   Traceback (most recent call last):
    ...
   ValueError: setting an array element with a sequence.

The problem lies in knowing if and how to assemble the individual elements
for each column into an appropriate mixin column.  The current code uses
numpy to perform this function on numerical or string types, but it obviously
does not handle mixin column types like |Quantity| or |SkyCoord|.

**Masking**

Mixin columns do not support masking, but there is limited support for use of
mixins within a masked table.  In this case a ``mask`` attribute is assigned to
the mixin column object.  This ``mask`` is a special object that is a boolean
array of ``False`` corresponding to the mixin data shape.  The ``mask`` looks
like a normal numpy array but an exception will be raised if ``True`` is assigned
to any element.  The consequences of the limitation are most obvious in the
high-level table operations.

**High-level table operations**

The table below gives a summary of support for high-level operations on tables
that contain mixin columns:

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Operation
     - Support
   * - :ref:`grouped-operations`
     - Not implemented yet, but no fundamental limitation
   * - :ref:`stack-vertically`
     - Not implemented yet, pending definition of generic concatenation protocol
   * - :ref:`stack-horizontally`
     - Works if output mixin columns do not require masking
   * - :ref:`table-join`
     - Works if output mixin columns do not require masking; no mixin key
       columns allowed
   * - :ref:`unique-rows`
     - Not implemented yet, uses grouped operations

**Mixin column attributes**

For mixin columns the column attributes ``name``, ``unit``, ``dtype``,
``format``, ``description`` and ``meta`` are stored in the ``info`` attribute



dictionary called `_astropy_column_attrs`.  These attributes can be manipulated
with the functions ``col_getattr`` and ``col_setattr`` which are available in
the ``astropy.table.column`` module.  These methods are not part of
the astropy public API and are likely to change in the future.

**ASCII table writing**

Mixin columns can be written out to file using the `astropy.io.ascii` module,
but the fast C-based writers are not available.  Instead the legacy pure-Python
writers will be used.


.. _mixin_protocol:

Mixin protocol
^^^^^^^^^^^^^^

A key idea behind mixin columns is that any class which satisfies a specified
protocol can be used.  That means many user-defined class objects which handle
array-like data can be used natively within a |Table|.  The protocol is
relatively simple and requires that a class behave like a minimal numpy array
with the following properties:

- Contains array-like data
- Supports getting data as a single item, slicing, or index array access
- Has a ``shape`` attribute
- Has a ``__len__`` method for length
- Has an ``info`` class descriptor which is a subclass of the
  ``astropy.utils.data_info.MixinInfo`` class.

The `Example: ArrayWrapper`_ section shows a working minimal example of a class
which can be used as a mixin column.  A `pandas.Series
<http://pandas.pydata.org/pandas-docs/stable/generated/pandas.Series.html>`_
object can function as a mixin column as well.

Other interesting possibilities for mixin columns include:

- Columns which are dynamically computed as a function of other columns (AKA
  spreadsheet)
- Columns which are themselves a |Table|, i.e. nested tables.  A `proof of
  concept <https://github.com/astropy/astropy/pull/3963>`_ is available.

.. _arraywrapper_example:

Example: ArrayWrapper
^^^^^^^^^^^^^^^^^^^^^

The code listing below shows a example of a data container class which acts as
a mixin column class.  This class is a simple wrapper around a numpy array.  It
is used in the astropy mixin test suite and is fully compliant as a mixin
column.

::

  from astropy.utils.data_info import ParentDtypeInfo

  class ArrayWrapper(object):
      """
      Minimal mixin using a simple wrapper around a numpy array
      """
      info = ParentDtypeInfo()

      def __init__(self, data):
          self.data = np.array(data)
          if 'info' in getattr(data, '__dict__', ()):
              self.info = data.info

      def __getitem__(self, item):
          if isinstance(item, (int, np.integer)):
              out = self.data[item]
          else:
              out = self.__class__(self.data[item])
              if 'info' in self.__dict__:
                  out.info = self.info
          return out

      def __setitem__(self, item, value):
          self.data[item] = value

      def __len__(self):
          return len(self.data)

      @property
      def dtype(self):
          return self.data.dtype

      @property
      def shape(self):
          return self.data.shape

      def __repr__(self):
          return ("<{0} name='{1}' data={2}>"
                  .format(self.__class__.__name__, self.info.name, self.data))
