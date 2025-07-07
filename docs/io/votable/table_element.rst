.. doctest-skip-all

Table element
-------------

VOTable files can contain ``RESOURCE`` elements, each of
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

Datatype Mappings
^^^^^^^^^^^^^^^^^

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
^^^^^^^^^^^^^^^^^^^^^

To look up more information about a field in a table, you can use the
`~astropy.io.votable.tree.TableElement.get_field_by_id` method, which returns
the `~astropy.io.votable.tree.Field` object with the given ID.

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

Building a New Table from Scratch
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is also possible to build a new table, define some field datatypes,
and populate it with data.

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

Missing Values
^^^^^^^^^^^^^^

Any value in the table may be "missing". `astropy.io.votable` stores
a  ``numpy`` masked array in each `~astropy.io.votable.tree.TableElement`
instance. This behaves like an ordinary ``numpy`` masked array, except
for variable-length fields. For those fields, the datatype of the
column is "object" and another ``numpy`` masked array is stored there.
Therefore, operations on variable-length columns will not work — this
is because variable-length columns are not directly supported
by ``numpy`` masked arrays.
