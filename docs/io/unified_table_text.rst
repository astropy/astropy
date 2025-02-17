.. _table_io_ascii:

ASCII Formats
-------------

The :meth:`~astropy.table.Table.read` and
:meth:`~astropy.table.Table.write` methods can be used to read and write formats
supported by `astropy.io.ascii`.

Use ``format='ascii'`` in order to interface to the generic
:func:`~astropy.io.ascii.read` and :func:`~astropy.io.ascii.write`
functions from `astropy.io.ascii`. When reading a table, this means
that all supported ASCII table formats will be tried in order to successfully
parse the input.

Examples
^^^^^^^^

..
  EXAMPLE START
  Reading and Writing ASCII Formats

To read and write formats supported by `astropy.io.ascii`:

.. doctest-skip::

  >>> t = Table.read('astropy/io/ascii/tests/t/latex1.tex', format='ascii')
  >>> print(t)
  cola colb colc
  ---- ---- ----
     a    1    2
     b    3    4

When writing a table with ``format='ascii'`` the output is a basic
character-delimited file with a single header line containing the
column names.

All additional arguments are passed to the `astropy.io.ascii`
:func:`~astropy.io.ascii.read` and :func:`~astropy.io.ascii.write`
functions. Further details are available in the sections on
:ref:`io_ascii_read_parameters` and :ref:`io_ascii_write_parameters`. For
example, to change the column delimiter and the output format for the ``colc``
column use:

.. doctest-skip::

  >>> t.write(sys.stdout, format='ascii', delimiter='|', formats={'colc': '%0.2f'})
  cola|colb|colc
  a|1|2.00
  b|3|4.00


.. note::

   When specifying an ASCII table format using the unified interface, the
   format name is prefixed with ``ascii`` in order to identify the format as
   ASCII-based. Compare the table above to the `astropy.io.ascii` list of
   :ref:`supported formats <supported_formats>` where the prefix is not
   needed. Therefore the following are equivalent:

.. doctest-skip::

     >>> dat = ascii.read('file.dat', format='daophot')
     >>> dat = Table.read('file.dat', format='ascii.daophot')

.. attention:: **ECSV is recommended**

   For writing and reading tables to ASCII in a way that fully reproduces the
   table data, types, and metadata (i.e., the table will "round-trip"), we
   highly recommend using the :ref:`ecsv_format`. This writes the actual data
   in a space-delimited format (the ``basic`` format) that any ASCII table
   reader can parse, but also includes metadata encoded in a comment block that
   allows full reconstruction of the original columns. This includes support
   for :ref:`ecsv_format_mixin_columns` (such as
   `~astropy.coordinates.SkyCoord` or `~astropy.time.Time`) and
   :ref:`ecsv_format_masked_columns`.

..
  EXAMPLE END
