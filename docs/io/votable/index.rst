.. doctest-skip-all

.. include:: references.txt

.. _astropy-io-votable:

*******************************************
VOTable XML Handling (`astropy.io.votable`)
*******************************************

Introduction
============

The `astropy.io.votable` sub-package converts VOTable XML files to and
from ``numpy`` record arrays. This subpackage was originally developed
as ``vo.table``.

Getting Started
===============

This section provides a quick introduction of using :mod:`astropy.io.votable`. The
goal is to demonstrate the package's basic features without getting into too
much detail.

.. note::

    If you want to read or write a single table in VOTable format, the
    recommended method is via the high-level :ref:`table_io`. In particular
    see the :ref:`Unified I/O VOTables <table_io_votable>` section.

Reading a VOTable File
----------------------

To read in a VOTable file, pass a file path to
`~astropy.io.votable.parse`::

    from astropy.io.votable import parse
    votable = parse("votable.xml")

``votable`` is a `~astropy.io.votable.tree.VOTableFile` object, which
can be used to retrieve and manipulate the data and save it back out
to disk.

VOTable files are made up of nested ``RESOURCE`` elements, each of
which may contain one or more ``TABLE`` elements. The ``TABLE``
elements contain the arrays of data.

To get at the ``TABLE`` elements, you can write a loop over the
resources in the ``VOTABLE`` file::

    for resource in votable.resources:
        for table in resource.tables:
            # ... do something with the table ...
            pass

However, if the nested structure of the resources is not important,
you can use `~astropy.io.votable.tree.VOTableFile.iter_tables` to
return a flat list of all tables::

    for table in votable.iter_tables():
        # ... do something with the table ...
        pass

Finally, if you expect only one table in the file, it might be most convenient
to use `~astropy.io.votable.tree.VOTableFile.get_first_table`::

  table = votable.get_first_table()

Alternatively, there is a convenience method to parse a VOTable file and
return the first table all in one step::

  from astropy.io.votable import parse_single_table
  table = parse_single_table("votable.xml")

From a `~astropy.io.votable.tree.TableElement` object, you can get the data itself
in the ``array`` member variable::

  data = table.array

This data is a ``numpy`` record array.

The columns get their names from both the ``ID`` and ``name``
attributes of the ``FIELD`` elements in the ``VOTABLE`` file.

Examples
^^^^^^^^

..
  EXAMPLE START
  Reading a VOTable File with astropy.io.votable

Suppose we had a ``FIELD`` specified as follows:

.. code-block:: xml

   <FIELD ID="Dec" name="dec_targ" datatype="char" ucd="POS_EQ_DEC_MAIN"
          unit="deg">
    <DESCRIPTION>
     representing the ICRS declination of the center of the image.
    </DESCRIPTION>
   </FIELD>

.. note::

    The mapping from VOTable ``name`` and ``ID`` attributes to ``numpy``
    dtype ``names`` and ``titles`` is highly confusing.

    In VOTable, ``ID`` is guaranteed to be unique, but is not
    required. ``name`` is not guaranteed to be unique, but is
    required.

    In ``numpy`` record dtypes, ``names`` are required to be unique and
    are required. ``titles`` are not required, and are not required
    to be unique.

    Therefore, VOTable's ``ID`` most closely maps to ``numpy``'s
    ``names``, and VOTable's ``name`` most closely maps to ``numpy``'s
    ``titles``. However, in some cases where a VOTable ``ID`` is not
    provided, a ``numpy`` ``name`` will be generated based on the VOTable
    ``name``. Unfortunately, VOTable fields do not have an attribute
    that is both unique and required, which would be the most
    convenient mechanism to uniquely identify a column.

    When converting from an `astropy.io.votable.tree.TableElement` object to
    an `astropy.table.Table` object, you can specify whether to give
    preference to ``name`` or ``ID`` attributes when naming the
    columns. By default, ``ID`` is given preference. To give
    ``name`` preference, pass the keyword argument
    ``use_names_over_ids=True``::

      >>> votable.get_first_table().to_table(use_names_over_ids=True)

This column of data can be extracted from the record array using::

  >>> table.array['dec_targ']
  array([17.15153360566, 17.15153360566, 17.15153360566, 17.1516686826,
         17.1516686826, 17.1516686826, 17.1536197136, 17.1536197136,
         17.1536197136, 17.15375479055, 17.15375479055, 17.15375479055,
         17.1553884541, 17.15539736932, 17.15539752176,
         17.25736014763,
         # ...
         17.2765703], dtype=object)

or equivalently::

  >>> table.array['Dec']
  array([17.15153360566, 17.15153360566, 17.15153360566, 17.1516686826,
         17.1516686826, 17.1516686826, 17.1536197136, 17.1536197136,
         17.1536197136, 17.15375479055, 17.15375479055, 17.15375479055,
         17.1553884541, 17.15539736932, 17.15539752176,
         17.25736014763,
         # ...
         17.2765703], dtype=object)

..
  EXAMPLE END

Building a New Table from Scratch
---------------------------------

It is also possible to build a new table, define some field datatypes,
and populate it with data.

Example
^^^^^^^

..
  EXAMPLE START
  Building a New Table from a VOTable File

To build a new table from a VOTable file::

  from astropy.io.votable.tree import VOTableFile, Resource, TableElement, Field

  # Create a new VOTable file...
  votable = VOTableFile()

  # ...with one resource...
  resource = Resource()
  votable.resources.append(resource)

  # ... with one table
  table = TableElement(votable)
  resource.tables.append(table)

  # Define some fields
  table.fields.extend([
          Field(votable, name="filename", datatype="char", arraysize="*"),
          Field(votable, name="matrix", datatype="double", arraysize="2x2")])

  # Now, use those field definitions to create the numpy record arrays, with
  # the given number of rows
  table.create_arrays(2)

  # Now table.array can be filled with data
  table.array[0] = ('test1.xml', [[1, 0], [0, 1]])
  table.array[1] = ('test2.xml', [[0.5, 0.3], [0.2, 0.1]])

  # Now write the whole thing to a file.
  # Note, we have to use the top-level votable file object
  votable.to_xml("new_votable.xml")

..
  EXAMPLE END

Outputting a VOTable File
-------------------------

This section describes writing table data in the VOTable format using the
`~astropy.io.votable` package directly. For some cases, however, the high-level
:ref:`table_io` will often suffice and is somewhat more convenient to use. See
the :ref:`Unified I/O VOTable <table_io_votable>` section for details.

To save a VOTable file, call the
`~astropy.io.votable.tree.VOTableFile.to_xml` method. It accepts
either a string or Unicode path, or a Python file-like object::

  votable.to_xml('output.xml')

There are a number of data storage formats supported by
`astropy.io.votable`. The ``TABLEDATA`` format is XML-based and
stores values as strings representing numbers. The ``BINARY`` format
is more compact, and stores numbers in base64-encoded binary. VOTable
version 1.3 adds the ``BINARY2`` format, which allows for masking of
any data type, including integers and bit fields which cannot be
masked in the older ``BINARY`` format. The storage format can be set
on a per-table basis using the `~astropy.io.votable.tree.TableElement.format`
attribute, or globally using the
`~astropy.io.votable.tree.VOTableFile.set_all_tables_format` method::

  votable.get_first_table().format = 'binary'
  votable.set_all_tables_format('binary')
  votable.to_xml('binary.xml')

Using `astropy.io.votable`
==========================

Standard Compliance
-------------------

`astropy.io.votable.tree.TableElement` supports the `VOTable Format Definition
Version 1.1
<https://www.ivoa.net/documents/REC/VOTable/VOTable-20040811.html>`_,
`Version 1.2
<https://www.ivoa.net/documents/VOTable/20091130/REC-VOTable-1.2.html>`_,
`Version 1.3
<https://www.ivoa.net/documents/VOTable/20130920/REC-VOTable-1.3-20130920.html>`_,
and `Version 1.4
<https://www.ivoa.net/documents/VOTable/20191021/REC-VOTable-1.4-20191021.html>`_.
Some flexibility is provided to support the 1.0 draft version and
other nonstandard usage in the wild, see :ref:`verifying-votables` for more
details.

.. note::

  Each warning and VOTABLE-specific exception emitted has a number and
  is documented in more detail in :ref:`warnings` and
  :ref:`exceptions`.

Output always conforms to the 1.1, 1.2, 1.3, or 1.4 spec, depending on the
input.

.. _verifying-votables:

Verifying VOTables
^^^^^^^^^^^^^^^^^^

Many VOTable files in the wild do not conform to the VOTable specification. You
can set what should happen when a violation is encountered with the ``verify``
keyword, which can take three values:

    * ``'ignore'`` - Attempt to parse the VOTable silently. This is the default
      setting.
    * ``'warn'`` - Attempt to parse the VOTable, but raise appropriate
      :ref:`warnings`. It is possible to limit the number of warnings of the
      same type to a maximum value using the
      `astropy.io.votable.exceptions.conf.max_warnings
      <astropy.io.votable.exceptions.Conf.max_warnings>` item in the
      :ref:`astropy_config`.
    * ``'exception'`` - Do not parse the VOTable and raise an exception.

The ``verify`` keyword can be used with the :func:`~astropy.io.votable.parse`
or :func:`~astropy.io.votable.parse_single_table` functions::

  from astropy.io.votable import parse
  votable = parse("votable.xml", verify='warn')

It is possible to change the default ``verify`` value through the
`astropy.io.votable.conf.verify <astropy.io.votable.Conf.verify>` item in the
:ref:`astropy_config`.

Note that ``'ignore'`` or ``'warn'``  mean that ``astropy`` will attempt to
parse the VOTable, but if the specification has been violated then success
cannot be guaranteed.

It is good practice to report any errors to the author of the application that
generated the VOTable file to bring the file into compliance with the
specification.

Missing Values
--------------

Any value in the table may be "missing". `astropy.io.votable` stores
a  ``numpy`` masked array in each `~astropy.io.votable.tree.TableElement`
instance. This behaves like an ordinary ``numpy`` masked array, except
for variable-length fields. For those fields, the datatype of the
column is "object" and another ``numpy`` masked array is stored there.
Therefore, operations on variable-length columns will not work — this
is because variable-length columns are not directly supported
by ``numpy`` masked arrays.

Datatype Mappings
-----------------

The datatype specified by a ``FIELD`` element is mapped to a ``numpy``
type according to the following table:

  ================================ =========================
  VOTABLE type                     NumPy type
  ================================ =========================
  boolean                          b1
  -------------------------------- -------------------------
  bit                              b1
  -------------------------------- -------------------------
  unsignedByte                     u1
  -------------------------------- -------------------------
  char (*variable length*)         O - A ``bytes()`` object.
  -------------------------------- -------------------------
  char (*fixed length*)            S
  -------------------------------- -------------------------
  unicodeChar (*variable length*)  O - A `str` object
  -------------------------------- -------------------------
  unicodeChar (*fixed length*)     U
  -------------------------------- -------------------------
  short                            i2
  -------------------------------- -------------------------
  int                              i4
  -------------------------------- -------------------------
  long                             i8
  -------------------------------- -------------------------
  float                            f4
  -------------------------------- -------------------------
  double                           f8
  -------------------------------- -------------------------
  floatComplex                     c8
  -------------------------------- -------------------------
  doubleComplex                    c16
  ================================ =========================

If the field is a fixed-size array, the data is stored as a ``numpy``
fixed-size array.

If the field is a variable-size array (that is, ``arraysize`` contains
a '*'), the cell will contain a Python list of ``numpy`` values. Each
value may be either an array or scalar depending on the ``arraysize``
specifier.

Examining Field Types
---------------------

To look up more information about a field in a table, you can use the
`~astropy.io.votable.tree.TableElement.get_field_by_id` method, which returns
the `~astropy.io.votable.tree.Field` object with the given ID.

Example
^^^^^^^

..
  EXAMPLE START
  Examining Field Types in VOTables with astropy.io.votable

To look up more information about a field::

  >>> field = table.get_field_by_id('Dec')
  >>> field.datatype
  'char'
  >>> field.unit
  'deg'

.. note::
   Field descriptors should not be mutated. To change the set of
   columns, convert the Table to an `astropy.table.Table`, make the
   changes, and then convert it back.

..
  EXAMPLE END

.. _votable-serialization:

Data Serialization Formats
--------------------------

VOTable supports a number of different serialization formats.

- `TABLEDATA
  <http://www.ivoa.net/documents/VOTable/20130920/REC-VOTable-1.3-20130920.html#ToC36>`__
  stores the data in pure XML, where the numerical values are written
  as human-readable strings.

- `BINARY
  <http://www.ivoa.net/documents/VOTable/20130920/REC-VOTable-1.3-20130920.html#ToC38>`__
  is a binary representation of the data, stored in the XML as an
  opaque ``base64``-encoded blob.

- `BINARY2
  <http://www.ivoa.net/documents/VOTable/20130920/REC-VOTable-1.3-20130920.html#ToC39>`__
  was added in VOTable 1.3, and is identical to "BINARY", except that
  it explicitly records the position of missing values rather than
  identifying them by a special value.

- `FITS
  <http://www.ivoa.net/documents/VOTable/20130920/REC-VOTable-1.3-20130920.html#ToC37>`__
  stores the data in an external FITS file. This serialization is not
  supported by the `astropy.io.votable` writer, since it requires
  writing multiple files.

- ``PARQUET``
  stores the data in an external PARQUET file, similar to FITS serialization.
  Reading and writing is fully supported by the `astropy.io.votable` writer and
  the `astropy.io.votable.parse` reader. The parquet file can be
  referenced with either absolute and relative paths. The parquet
  serialization can be used as part of the unified Table I/O (see next
  section), by setting the ``format`` argument to ``'votable.parquet'``.

The serialization format can be selected in two ways:

    1) By setting the ``format`` attribute of a
    `astropy.io.votable.tree.TableElement` object::

        votable.get_first_table().format = "binary"
        votable.to_xml("new_votable.xml")

    2) By overriding the format of all tables using the
    ``tabledata_format`` keyword argument when writing out a VOTable
    file::

        votable.to_xml("new_votable.xml", tabledata_format="binary")

Converting to/from an `astropy.table.Table`
-------------------------------------------

The VOTable standard does not map conceptually to an
`astropy.table.Table`. However, a single table within the ``VOTable``
file may be converted to and from an `astropy.table.Table`::

  from astropy.io.votable import parse_single_table
  table = parse_single_table("votable.xml").to_table()

As a convenience, there is also a function to create an entire VOTable
file with just a single table::

  from astropy.io.votable import from_table, writeto
  votable = from_table(table)
  writeto(votable, "output.xml")

.. note::

  By default, ``to_table`` will use the ``ID`` attribute from the files to
  create the column names for the `~astropy.table.Table` object. However,
  it may be that you want to use the ``name`` attributes instead. For this,
  set the ``use_names_over_ids`` keyword to `True`. Note that since field
  ``names`` are not guaranteed to be unique in the VOTable specification,
  but column names are required to be unique in ``numpy`` structured arrays (and
  thus `astropy.table.Table` objects), the names may be renamed by appending
  numbers to the end in some cases.

Performance Considerations
--------------------------

File reads will be moderately faster if the ``TABLE`` element includes
an nrows_ attribute. If the number of rows is not specified, the
record array must be resized repeatedly during load.

.. _nrows: http://www.ivoa.net/documents/REC/VOTable/VOTable-20040811.html#ToC10

.. _votable_mivot:

Reading and writing VO model annotations
========================================

Introduction
------------
Model Instances in VOTables (`MIVOT <https://ivoa.net/documents/MIVOT/20230620/REC-mivot-1.0.pdf>`_)
defines a syntax to map VOTable data to any model serialised in VO-DML (Virtual Observatory Data Modeling Language).
This annotation schema operates as a bridge between data and the models. It associates both column/param metadata and data
from the VOTable to the data model elements (class, attributes, types, etc.). It also brings up VOTable data or
metadata that were possibly missing in the table, e.g., coordinate system description, or curation tracing.
The data model elements are grouped in an independent annotation block complying with the MIVOT XML schema which
is added as an extra resource above the table element.
The MIVOT syntax allows to describe a data structure as a hierarchy of classes.
It is also able to represent relations and compositions between them. It can moreover build up data model objects by
aggregating instances from different tables of the VOTable.

Astropy implementation
----------------------
The purpose of Astropy is not to process VO annotations.
It is just to allow related packages to get and set MIVOT blocks from/into VOTables.
For this reason, in this implementation MIVOT annotations are both imported and exported as strings.
The current implementation prevents client code from injecting into VOTables strings
that are not MIVOT serializations.

MivotBlock implementation:

- MIVOT blocks are handled by the :class:`astropy.io.votable.tree.MivotBlock` class.
- A MivotBlock instance can only be carried by a resource with "type=meta".
- This instance holds the XML mapping block as a string.
- MivotBlock objects are instanced by the Resource parser.
- The MivotBlock class has its own logic that operates both parsing and IO functionalities.

Example
^^^^^^^

.. code-block:: xml

       <VOTABLE xmlns="http://www.ivoa.net/xml/VOTable/v1.3"
         xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" version="1.3">
         <RESOURCE>
           <RESOURCE type="meta">
             <VODML xmlns="http://www.ivoa.net/xml/mivot">
              ...
             </VODML>
           </RESOURCE>
           <TABLE name="myDataTable">
            ....
           </TABLE>
         </RESOURCE>
    </VOTABLE>

Reading a VOTable containing a MIVOT block
------------------------------------------

To read in a VOTable file containing or not a MIVOT Resource, pass a file path to`~astropy.io.votable.parse`:

.. code-block:: python

   >>> from astropy.io.votable import parse
   >>> from astropy.utils.data import get_pkg_data_filename
   >>> votable = parse(get_pkg_data_filename("data/test.order.xml", package="astropy.io.votable.tests"))
   <VODML xmlns="http://www.ivoa.net/xml/mivot">
   </VODML>

   <VODML xmlns="http://www.ivoa.net/xml/mivot">
   </VODML>

The parse function will call the MIVOT parser if it detects a MIVOT block.

Building a Resource containing a MIVOT block
--------------------------------------------

Construct the MIVOT block by passing the XML block as a parameter:

.. code-block:: python

   >>> from astropy.io.votable import tree
   >>> from astropy.io.votable.tree import MivotBlock, Resource, VOTableFile
   >>> mivot_block = MivotBlock("""
   <VODML xmlns="http://www.ivoa.net/xml/mivot" >
      <REPORT status="OK">Unit test mapping block</REPORT>
      <GLOBALS>  </GLOBALS>
   </VODML>
   """)

Build a new resource:

.. code-block:: python

   >>> mivot_resource = Resource()

Give it the type meta:

.. code-block:: python

   >>> mivot_resource.type = "meta"

Then add it the MIVOT block:

.. code-block:: python

   >>> mivot_resource.mivot_block = mivot_block

Now you have a MIVOT resource that you can add to an object Resource creating a new Resource:

.. code-block:: python

   >>> votable = VOTableFile()
   >>> r1 = Resource()
   >>> r1.type = "results"
   >>> r1.resources.append(mivot_resource)

You can add an `astropy.io.votable.tree.TableElement` to the resource:

.. code-block:: python

   >>> table = tree.TableElement(votable)
   >>> r1.tables.append(t1)
   >>> votable.resources.append(r1)
   >>> for resource in votable.resources:
   ...     print(resource.mivot_block.content)
   <VODML xmlns="http://www.ivoa.net/xml/mivot" >
      <REPORT status="OK">Unit test mapping block</REPORT>
      <GLOBALS>  </GLOBALS>
   </VODML>


See Also
========

- `VOTable Format Definition Version 1.1
  <https://www.ivoa.net/documents/REC/VOTable/VOTable-20040811.html>`_

- `VOTable Format Definition Version 1.2
  <https://www.ivoa.net/documents/VOTable/20091130/REC-VOTable-1.2.html>`_

- `VOTable Format Definition Version 1.3
  <https://www.ivoa.net/documents/VOTable/20130920/REC-VOTable-1.3-20130920.html>`_

- `VOTable Format Definition Version 1.4
  <https://www.ivoa.net/documents/VOTable/20191021/REC-VOTable-1.4-20191021.html>`_

- `MIVOT Recommendation Version 1.0
  <https://ivoa.net/documents/MIVOT/20230620/REC-mivot-1.0.pdf>`_

.. note that if this section gets too long, it should be moved to a separate
   doc page - see the top of performance.inc.rst for the instructions on how to do
   that
.. include:: performance.inc.rst

Reference/API
=============

.. toctree::
   :maxdepth: 2

   ref_api
   api_exceptions
