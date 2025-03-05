.. doctest-skip-all

.. include:: references.txt

.. _astropy-io-votable:

***************************************
VOTable Handling (`astropy.io.votable`)
***************************************

Introduction
============

The `astropy.io.votable` sub-package converts VOTable XML files to and
from ``numpy`` record arrays.

.. note::

    If you want to read or write a single table in VOTable format, the
    recommended method is via the :ref:`table_io` interface. In particular,
    see the :ref:`Unified I/O VO Tables <table_io_votable>` section.

Getting Started
===============

Reading a VOTable File
----------------------

To read a VOTable file, pass a file path to
`~astropy.io.votable.parse`::

    from astropy.io.votable import parse
    votable = parse("votable.xml")

``votable`` is a `~astropy.io.votable.tree.VOTableFile` object, which
can be used to retrieve and manipulate the data and save it back out
to disk.

Writing a VOTable File
----------------------

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

The VOTable elements
====================

VOTables are built from nested elements. Let's for example build a
votable containing an ``INFO`` element::

  >>> from astropy.io.votable.tree import VOTableFile, Info
  >>> vot = VOTableFile()
  >>> vot.infos.append(Info(name="date_obs", value="2025-01-01"))

These elements can be:

- `~astropy.io.votable.tree.CooSys`
- `~astropy.io.votable.tree.TimeSys`
- `~astropy.io.votable.tree.Info`
- `~astropy.io.votable.tree.Param`
- `~astropy.io.votable.tree.Group`
- `~astropy.io.votable.tree.Resource`
- `~astropy.io.votable.tree.Link`
- `~astropy.io.votable.tree.TableElement`
- `~astropy.io.votable.tree.Field`
- `~astropy.io.votable.tree.Values`
- `~astropy.io.votable.tree.MivotBlock`

Here are some detailed explanations on some of these elements:

.. toctree::
   :maxdepth: 1

   table_element
   coosys_element
   mivot_blocks

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
`Version 1.4
<https://www.ivoa.net/documents/VOTable/20191021/REC-VOTable-1.4-20191021.html>`_,
and `Version 1.5
<https://ivoa.net/documents/VOTable/20250116/REC-VOTable-1.5.html>`_,
Some flexibility is provided to support the 1.0 draft version and
other nonstandard usage in the wild, see :ref:`verifying-votables` for more
details.

.. note::

  Each warning and VOTABLE-specific exception emitted has a number and
  is documented in more detail in :ref:`warnings` and
  :ref:`exceptions`.

Output always conforms to the 1.1, 1.2, 1.3, 1.4, 1.5 spec, depending on the
input.

.. _verifying-votables:

Verifying VOTables
------------------

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


Data Origin
===========

.. _astropy-io-votable-dataorigin:

.. include:: dataorigin.rst


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

- `VOTable Format Definition Version 1.5
  <https://ivoa.net/documents/VOTable/20250116/REC-VOTable-1.5.html>`_

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
