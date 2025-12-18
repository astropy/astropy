
.. _table_io_votable:

VO Tables
---------

Reading/writing from/to `VO table <http://www.ivoa.net/documents/VOTable/>`_
files is supported with ``format='votable'``. In most cases, existing VO
tables should be automatically identified as such based on the header of the
file, but if not, or if writing to disk, then the format should be explicitly
specified.

Examples
^^^^^^^^

..
  EXAMPLE START
  Reading from and Writing to VO Tables

.. testsetup::
    >>> # Shortcut to set up a single table in a VO file
    >>> from astropy.table import Table
    >>> from astropy.io import votable
    >>> t1 = Table({'object': ['M31', 'Vega'], 'J': [10.0, 5.0]})
    >>> t1.meta['ID'] = 'twomass'
    >>> vot = votable.from_table(t1)
    >>> vot.to_xml('aj285677t3_votable.xml')
    >>> # Make and append a second table
    >>> from astropy.io.votable.tree import TableElement, Field
    >>> table = TableElement(votable)
    >>> table.fields.extend([Field(votable, name="filename", datatype="char", arraysize="*"),
    ...                      Field(votable, name="data", datatype="float", arraysize="2x2")])
    >>> table.create_arrays(1)
    >>> table.ID = "spitzer"
    >>> table.array[0] = ("test1.xml", [[1, 0], [0, 1]])
    >>> vot.resources[0].tables.append(table)
    >>> vot.to_xml('catalog.xml')

If a VO table file contains only a single table, then it can be read in with::

    >>> t = Table.read('aj285677t3_votable.xml')

If more than one table is present in the file, an error will be raised,
unless the table ID is specified via the ``table_id=`` argument::

    >>> t = Table.read('catalog.xml')
    Traceback (most recent call last):
    ...
    ValueError: Multiple tables found: table id should be set via the table_id= argument. The available tables are twomass, spitzer, or integers less than 2.

    >>> t = Table.read('catalog.xml', table_id='twomass')

To write to a new file, the ID of the table should also be specified (unless
``t.meta['ID']`` is defined)::

    >>> t.write('new_catalog.xml', table_id='updated_table', format='votable')

When writing, the ``compression=True`` argument can be used to force
compression of the data on disk, and the ``overwrite=True`` argument can be
used to overwrite an existing file.

.. testcleanup::

   >>> import os
   >>> os.remove('aj285677t3_votable.xml')
   >>> os.remove('catalog.xml')
   >>> os.remove('new_catalog.xml')

..
  EXAMPLE END
