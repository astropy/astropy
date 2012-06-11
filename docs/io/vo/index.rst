.. include:: references.txt

VOTable XML handling (`astropy.io.vo`)
======================================

Introduction
------------

The `astropy.io.vo` subpackage converts VOTable XML files to and from
Numpy record arrays.

Getting Started
---------------

Reading a VOTable file
^^^^^^^^^^^^^^^^^^^^^^

To read in a VOTable file, pass a file path to
`astropy.io.vo.table.parse`::

  from astropy.io.vo.table import parse
  votable = parse("votable.xml")

``votable`` is a `~astropy.io.vo.tree.VOTableFile` object, which can
be used to retrieve and manipulate the data and save it back out to
disk.

VOTable files are made up of nested ``RESOURCE`` elements, each of
which may contain one or more ``TABLE`` elements.  The ``TABLE``
elements contain the arrays of data.

To get at the ``TABLE`` elements, one can write a loop over the
resources in the ``VOTABLE`` file::

  for resource in votable.resources:
    for table in resource.tables:
      # ... do something with the table ...
      pass

However, if the nested structure of the resources is not important,
one can use `~astropy.io.vo.tree.VOTable.iter_tables` to return a flat
list of all tables::

  for table in votable.iter_tables():
    # ... do something with the table ...
    pass

Finally, if there is expected to be only one table in the file, it
might be simplest to just use
`~astropy.io.vo.tree.VOTable.get_first_table`::

  table = votable.get_first_table()

Even easier, there is a convenience method to parse a VOTable file and
return the first table all in one step::

  from astropy.io.vo.table import parse_single_table
  table = parse_single_table("votable.xml")

From a `~astropy.io.vo.tree.Table` object, one can get the data itself
in the ``array`` member variable::

  data = table.array

This data is a Numpy record array.  The columns get their names from
both the ``ID`` and ``name`` attributes of the ``FIELD`` elements in
the ``VOTABLE`` file.  For example, suppose we had a ``FIELD``
specified as follows:

.. code-block:: xml

   <FIELD ID="Dec" name="dec_targ" datatype="char" ucd="POS_EQ_DEC_MAIN"
          unit="deg">
    <DESCRIPTION>
     representing the ICRS declination of the center of the image.
    </DESCRIPTION>
   </FIELD>

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

Building a new table from scratch
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is also possible to build a new table, define some field datatypes
and populate it with data::

  from astropy.io.vo.tree import VOTableFile, Resource, Table, Field

  # Create a new VOTable file...
  votable = VOTableFile()

  # ...with one resource...
  resource = Resource()
  votable.resources.append(resource)

  # ... with one table
  table = Table(votable)
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

Outputting a VOTable file
^^^^^^^^^^^^^^^^^^^^^^^^^

To save a VOTable file, simply call the
`~astropy.io.vo.tree.VOTableFile.to_xml` method.  It accepts either a
string or unicode path, or a Python file-like object::

  votable.to_xml('output.xml')

There are currently two data storage formats supported by
`astropy.io.vo`.  The ``TABLEDATA`` format is XML-based and stores
values as strings representing numbers.  The ``BINARY`` format is more
compact, and stores numbers in base64-encoded binary.  The storage
format can be set on a per-table basis using the
`~astropy.io.vo.tree.Table.format` attribute, or globally using the
`~astropy.io.vo.tree.VOTableFile.set_all_tables_format` method::

  votable.get_first_table().format = 'binary'
  votable.set_all_tables_format('binary')
  votable.to_xml('binary.xml')

Using `io.vo`
-------------

Standard compliance
^^^^^^^^^^^^^^^^^^^

`astropy.io.vo.table` supports the `VOTable Format Definition Version
1.1
<http://www.ivoa.net/Documents/REC/VOTable/VOTable-20040811.html>`_
and `Version 1.2
<http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html>`_.
Some flexibility is provided to support the 1.0 draft version and
other non-standard usage in the wild.  To support these cases, set the
keyword argument ``pedantic`` to ``False`` when parsing.

.. note::

  Each warning and VOTABLE-specific exception emitted has a number and
  is documented in more detail in :ref:`warnings` and
  :ref:`exceptions`.

Output always conforms to the 1.1 or 1.2 spec, depending on the input.

.. _pedantic-mode:

Pedantic mode
_____________

Many VOTABLE files in the wild do not conform to the VOTABLE
specification.  If reading one of these files causes exceptions, you
may turn off pedantic mode in `astropy.io.vo` by passing
``pedantic=False`` to the `~astropy.io.vo.table.parse` or
`~astropy.io.vo.table.parse_single_table` functions::

  from astropy.io.vo.table import parse
  votable = parse("votable.xml", pedantic=False)

Note, however, that it is good practice to report these errors to the
author of the application that generated the VOTABLE file to bring the
file into compliance with the specification.

Even with ``pedantic`` turned off, many warnings may still be omitted.
These warnings are all of the type
`~astropy.io.vo.exceptions.VOTableSpecWarning` and can be turned off
using the standard Python `warnings` module.

Missing values
^^^^^^^^^^^^^^

Any value in the table may be "missing".  `astropy.io.vo.table` stores
a parallel array in each `~astropy.io.vo.tree.Table` instance called
`~astropy.io.vo.tree.Table.mask` to keep track of missing values.
This array is ``False`` anywhere the value is missing.

.. note::
   In the future, the ``array`` and ``mask`` members will likely be
   combined into a single masked record array.  There are
   implementation bugs in current versions of Numpy that prevent this
   at the moment.

Datatype mappings
^^^^^^^^^^^^^^^^^

The datatype specified by a ``FIELD`` element is mapped to a Numpy
type according to the following table:

  ================================ ========================================================================
  VOTABLE type                     Numpy type
  ================================ ========================================================================
  boolean                          b1
  -------------------------------- ------------------------------------------------------------------------
  bit                              b1
  -------------------------------- ------------------------------------------------------------------------
  unsignedByte                     u1
  -------------------------------- ------------------------------------------------------------------------
  char (*variable length*)         O - In Python 2.x, a `str` object; in 3.x, a `bytes` object.
  -------------------------------- ------------------------------------------------------------------------
  char (*fixed length*)            S
  -------------------------------- ------------------------------------------------------------------------
  unicodeChar (*variable length*)  O - In Python 2.x, a `unicode` object, in utf-16; in 3.x a `str` object
  -------------------------------- ------------------------------------------------------------------------
  unicodeChar (*fixed length*)     U
  -------------------------------- ------------------------------------------------------------------------
  short                            i2
  -------------------------------- ------------------------------------------------------------------------
  int                              i4
  -------------------------------- ------------------------------------------------------------------------
  long                             i8
  -------------------------------- ------------------------------------------------------------------------
  float                            f4
  -------------------------------- ------------------------------------------------------------------------
  double                           f8
  -------------------------------- ------------------------------------------------------------------------
  floatComplex                     c8
  -------------------------------- ------------------------------------------------------------------------
  doubleComplex                    c16
  ================================ ========================================================================

If the field is a fixed size array, the data is stored as a Numpy
fixed-size array.

If the field is a variable size array (that is ``arraysize`` contains
a '*'), the cell will contain a Python list of Numpy values.  Each
value may be either an array or scalar depending on the ``arraysize``
specifier.

Examining field types
^^^^^^^^^^^^^^^^^^^^^

To look up more information about a field in a table, one can use the
`~astropy.io.vo.tree.Table.get_field_or_param_by_id` method, which
returns the `~astropy.io.vo.tree.Field` object with the given ID.  For
example::

  >>> field = table.get_field_or_param_by_id('Dec')
  >>> field.datatype
  'char'
  >>> field.unit
  'deg'

.. note::
   Field descriptors should not be mutated -- they will have no effect
   on the record arrays storing the data.  This shortcoming will be
   addressed in a future version of `astropy.io.vo`.

Performance considerations
^^^^^^^^^^^^^^^^^^^^^^^^^^

File reads will be moderately faster if the ``TABLE`` element includes
an nrows_ attribute.  If the number of rows is not specified, the
record array must be resized repeatedly during load.

.. _nrows: http://www.ivoa.net/Documents/REC/VOTable/VOTable-20040811.html#ToC10

See Also
--------

- `VOTable Format Definition Version 1.1
  <http://www.ivoa.net/Documents/REC/VOTable/VOTable-20040811.html>`_

- `VOTable Format Definition Version 1.2
  <http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html>`_

Reference/API
-------------

.. automodapi:: astropy.io.vo
   :no-main-section:
   :subsections: table, tree, converters, ucd, unit, util, validator, xmlutil
   :no-inheritance-diagram:

astropy.io.vo.exceptions Module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 1

   api_exceptions.rst
