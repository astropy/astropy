.. |join| replace:: :func:`~astropy.table.join`

.. _mixin_columns:

Mixin Columns
*************

``astropy`` tables support the concept of "mixin columns", which
allows integration of appropriate non-|Column| based class objects within a
|Table| object. These mixin column objects are not converted in any way but are
used natively.

The available built-in mixin column classes are:

- |Quantity| and subclasses
- |SkyCoord| and coordinate frame classes
- |Time| and :class:`~astropy.time.TimeDelta`
- :class:`~astropy.coordinates.EarthLocation`
- `~astropy.table.NdarrayMixin`

Basic Example
=============

.. EXAMPLE START: Using Mixin Columns in Tables

As an example we can create a table and add a time column::

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

The important point here is that the ``time`` column is a bona fide |Time|
object::

  >>> t['time']
  <Time object: scale='utc' format='isot' value=['2001-01-02T12:34:56.000' '2001-02-03T00:01:02.000']>
  >>> t['time'].mjd  # doctest: +FLOAT_CMP
  array([51911.52425926, 51943.00071759])

.. EXAMPLE END

.. _quantity_and_qtable:

Quantity and QTable
===================

The ability to natively handle |Quantity| objects within a table makes it more
convenient to manipulate tabular data with units in a natural and robust way.
However, this feature introduces an ambiguity because data with a unit
(e.g., from a FITS binary table) can be represented as either a |Column| with a
``unit`` attribute or as a |Quantity| object. In order to cleanly resolve this
ambiguity, ``astropy`` defines a minor variant of the |Table| class called
|QTable|. The |QTable| class is exactly the same as |Table| except that
|Quantity| is the default for any data column with a defined unit.

If you take advantage of the |Quantity| infrastructure in your analysis, then
|QTable| is the preferred way to create tables with units. If instead you use
table column units more as a descriptive label, then the plain |Table| class is
probably the best class to use.

Example
-------

.. EXAMPLE START: Using Quantity Columns and QTables

To illustrate these concepts we first create a standard |Table| where we supply
as input a |Time| object and a |Quantity| object with units of ``m / s``. In
this case the quantity is converted to a |Column| (which has a ``unit``
attribute but does not have all of the features of a |Quantity|)::

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

So instead let's do the same thing using a |QTable|::

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

You can conveniently convert |Table| to |QTable| and vice-versa::

  >>> qt2 = QTable(t)
  >>> type(qt2['velocity'])
  <class 'astropy.units.quantity.Quantity'>

  >>> t2 = Table(qt2)
  >>> type(t2['velocity'])
  <class 'astropy.table.column.Column'>

.. Note::

   To summarize: the **only** difference between `~astropy.table.QTable` and
   `~astropy.table.Table` is the behavior when adding a column that has a
   specified unit. With `~astropy.table.QTable` such a column is always
   converted to a `~astropy.units.Quantity` object before being added to the
   table. Likewise if a unit is specified for an existing unit-less
   `~astropy.table.Column` in a `~astropy.table.QTable`, then the column is
   converted to `~astropy.units.Quantity`.

   The converse is that if you add a `~astropy.units.Quantity` column to an
   ordinary `~astropy.table.Table` then it gets converted to an ordinary
   `~astropy.table.Column` with the corresponding ``unit`` attribute.

.. attention::

   When a column of ``int`` ``dtype`` is converted to `~astropy.units.Quantity`,
   its ``dtype`` is converted to ``float``.

   For example, for a quality flag column of ``int``, if it is
   assigned with the :ref:`dimensionless unit <doc_dimensionless_unit>`, it will still
   be converted to ``float``. Therefore such columns typically should not be
   assigned with any unit.

.. EXAMPLE END

.. _mixin_attributes:

Mixin Attributes
================

The usual column attributes ``name``, ``dtype``, ``unit``, ``format``, and
``description`` are available in any mixin column via the ``info`` property::

  >>> qt['velocity'].info.name
  'velocity'

This ``info`` property is a key bit of glue that allows a non-|Column| object
to behave much like a |Column|.

The same ``info`` property is also available in standard
`~astropy.table.Column` objects. These ``info`` attributes like
``t['a'].info.name`` refer to the direct `~astropy.table.Column`
attribute (e.g., ``t['a'].name``) and can be used interchangeably.
Likewise in a `~astropy.units.Quantity` object, ``info.dtype``
attribute refers to the native ``dtype`` attribute of the object.

.. Note::

   When writing generalized code that handles column objects which
   might be mixin columns, you must *always* use the ``info``
   property to access column attributes.

.. _details_and_caveats:

Details and Caveats
===================

Most common table operations behave as expected when mixin columns are part of
the table. However, there are limitations in the current implementation.

**Adding or inserting a row**

Adding or inserting a row works as expected only for mixin classes that are
mutable (data can be changed internally) and that have an ``insert()`` method.
Adding rows to a |Table| with |Quantity|, |Time| or |SkyCoord| columns does
work.

**Masking**

Masking of mixin columns is enabled by the |Masked| class. See
:ref:`utils-masked` for details.

**ASCII table writing**

Tables with mixin columns can be written out to file using the
`astropy.io.ascii` module, but the fast C-based writers are not available.
Instead, the pure-Python writers will be used. For writing tables with mixin
columns it is recommended to use the :ref:`ecsv_format`. This will fully
serialize the table data and metadata, allowing full "round-trip" of the table
when it is read back.

**Binary table writing**

Tables with mixin columns can be written to binary files using FITS, HDF5 and
Parquet formats. These can be read back to recover exactly the original |Table|
including mixin columns and metadata. See :ref:`table_io` for details.

.. _mixin_protocol:

Mixin Protocol
==============

A key idea behind mixin columns is that any class which satisfies a specified
protocol can be used. That means many user-defined class objects which handle
array-like data can be used natively within a |Table|. The protocol is
relatively concise and requires that a class behave like a minimal ``numpy``
array with the following properties:

- Contains array-like data.
- Implements ``__getitem__()`` to support getting data as a
  single item, slicing, or index array access.
- Has a ``shape`` attribute.
- Has a ``__len__()`` method for length.
- Has an ``info`` class descriptor which is a subclass of the
  :class:`astropy.utils.data_info.MixinInfo` class.

The `Example: ArrayWrapper`_ section shows a minimal working example of a class
which can be used as a mixin column. A :class:`pandas.Series` object can
function as a mixin column as well.

Other interesting possibilities for mixin columns include:

- Columns which are dynamically computed as a function of other columns (AKA
  spreadsheet).
- Columns which are themselves a |Table| (i.e., nested tables). A `proof of
  concept <https://github.com/astropy/astropy/pull/3963>`_ is available.

new_like() method
-----------------

In order to support high-level operations like :func:`~astropy.table.join` and
:func:`~astropy.table.vstack`, a mixin class must provide a ``new_like()``
method in the ``info`` class descriptor. A key part of the functionality is to
ensure that the input column metadata are merged appropriately and that the
columns have consistent properties such as the shape.

A mixin class that provides ``new_like()`` must also implement
``__setitem__()`` to support setting via a single item, slicing, or index
array.

The ``new_like()`` method has the following signature::

    def new_like(self, cols, length, metadata_conflicts='warn', name=None):
        """
        Return a new instance of this class which is consistent with the
        input ``cols`` and has ``length`` rows.

        This is intended for creating an empty column object whose elements can
        be set in-place for table operations like join or vstack.

        Parameters
        ----------
        cols : list
            List of input columns
        length : int
            Length of the output column object
        metadata_conflicts : {'warn', 'error', 'silent'}
            How to handle metadata conflicts
        name : str
            Output column name

        Returns
        -------
        col : object
            New instance of this class consistent with ``cols``
        """

Examples of this are found in the `~astropy.table.column.ColumnInfo` and
`~astropy.units.quantity.QuantityInfo` classes.


.. _arraywrapper_example:

Example: ArrayWrapper
=====================

The code listing below shows an example of a data container class which acts as
a mixin column class. This class is a wrapper around a |ndarray|. It is used in
the ``astropy`` mixin test suite and is fully compliant as a mixin column.

::

  from astropy.utils.data_info import ParentDtypeInfo

  class ArrayWrapper:
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
          return f"<{self.__class__.__name__} name='{self.info.name}' data={self.data}>"

.. _table_mixin_registry:

Registering array-like objects as mixin columns
===============================================

In some cases, you may want to directly add an array-like
object as a table column while maintaining the original object properties
(instead of the default conversion of the object to a `~astropy.table.Column`).
This is done by registering the object class as a mixin column and
defining a handler which allows `~astropy.table.Table` to treat that object
class as a mixin similar to the built-in mixin columns such as `~astropy.time.Time`
or `~astropy.units.quantity.Quantity`.

This can be done for data classes that are defined in third-party packages and which
you have no control over. As an example, we define a class
that is not numpy-like and stores the data in a private attribute::

    >>> import numpy as np
    >>> class ExampleDataClass:
    ...     def __init__(self):
    ...         self._data = np.array([0, 1, 3, 4], dtype=float)

By default, this cannot be used as a table column::

    >>> t = Table()
    >>> t['data'] = ExampleDataClass()
    Traceback (most recent call last):
    ...
    TypeError: Empty table cannot have column set to scalar value

However, you can create a function (or 'handler') which takes
an instance of the data class you want to have automatically
handled and returns a mixin column::

    >>> from astropy.table.table_helpers import ArrayWrapper
    >>> def handle_example_data_class(obj):
    ...     return ArrayWrapper(obj._data)

You can then register this by providing the fully qualified name
of the class and the handler function::

    >>> from astropy.table.mixins.registry import register_mixin_handler
    >>> register_mixin_handler('__main__.ExampleDataClass', handle_example_data_class)
    >>> t['data'] = ExampleDataClass()
    >>> t
    <Table length=4>
      data
    float64
    -------
        0.0
        1.0
        3.0
        4.0

Because we defined the data class as part of the example
above, the fully qualified name starts with ``__main__``,
but for a class in a third-party package, this might look
like ``package.Class`` for example.
