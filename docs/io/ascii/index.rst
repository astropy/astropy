.. include:: references.txt

.. _io-ascii:

*********************************
ASCII Tables (`astropy.io.ascii`)
*********************************

Introduction
============

`astropy.io.ascii` provides methods for reading and writing a wide range of
ASCII data table formats via built-in :ref:`extension_reader_classes`. The
emphasis is on flexibility and convenience of use, although readers can
optionally use a less flexible C-based engine for reading and writing for
improved performance.

The following shows a few of the ASCII formats that are available, while the
section on `Supported formats`_ contains the full list.

* :class:`~astropy.io.ascii.Basic`: basic table with customizable delimiters and header configurations
* :class:`~astropy.io.ascii.Cds`: `CDS format table <http://vizier.u-strasbg.fr/doc/catstd.htx>`_ (also Vizier and ApJ machine readable tables)
* :class:`~astropy.io.ascii.Daophot`: table from the IRAF DAOphot package
* :class:`~astropy.io.ascii.Ecsv`: :ref:`ecsv_format` for lossless round-trip of data tables
* :class:`~astropy.io.ascii.FixedWidth`: table with fixed-width columns (see also :ref:`fixed_width_gallery`)
* :class:`~astropy.io.ascii.Ipac`: `IPAC format table <https://irsa.ipac.caltech.edu/applications/DDGEN/Doc/ipac_tbl.html>`_
* :class:`~astropy.io.ascii.HTML`: HTML format table contained in a <table> tag
* :class:`~astropy.io.ascii.Latex`: LaTeX table with datavalue in the ``tabular`` environment
* :class:`~astropy.io.ascii.Rdb`: tab-separated values with an extra line after the column definition line
* :class:`~astropy.io.ascii.SExtractor`: `SExtractor format table <http://www.astromatic.net/software/sextractor>`_

The strength of `astropy.io.ascii` is the support for astronomy-specific
formats (often with metadata) and specialized data types such as
:ref:`SkyCoord <astropy-coordinates-high-level>`, :ref:`Time
<astropy-time>`, and :ref:`Quantity <quantity>`. For reading or writing large
data tables in a generic format such as CSV, using the :ref:`Table - Pandas
interface <pandas>` is a recommended option to consider.

.. note::

    It is also possible (and encouraged) to use the functionality from
    :mod:`astropy.io.ascii` through a higher level interface in the
    :ref:`Data Tables <astropy-table>` package. See :ref:`table_io` for more details.

Getting Started
===============

Reading Tables
--------------

The majority of commonly encountered ASCII tables can be read with the
|read| function. Assume you have a file named ``sources.dat`` with the
following contents::

  obsid redshift  X      Y     object
  3102  0.32      4167  4085   Q1250+568-A
  877   0.22      4378  3892   "Source 82"

This table can be read with the following::

  >>> from astropy.io import ascii
  >>> data = ascii.read("sources.dat")  # doctest: +SKIP
  >>> print(data)                       # doctest: +SKIP
  obsid redshift  X    Y      object
  ----- -------- ---- ---- -----------
   3102     0.32 4167 4085 Q1250+568-A
    877     0.22 4378 3892   Source 82

The first argument to the |read| function can be the name of a file, a string
representation of a table, or a list of table lines. The return value
(``data`` in this case) is a :ref:`Table <astropy-table>` object.

By default, |read| will try to :ref:`guess the table format <guess_formats>`
by trying all of the `supported formats`_.

.. Warning::

   Guessing the file format is often slow for large files because the reader
   tries parsing the file with every allowed format until one succeeds.
   For large files it is recommended to disable guessing with ``guess=False``.

If guessing the format does not work, as in the case for unusually formatted
tables, you may need to give `astropy.io.ascii` additional hints about
the format.

Examples
^^^^^^^^

..
  EXAMPLE START
  Reading ASCII Tables Using astropy.io.ascii

For unusually formatted tables, give additional hints about the format::

   >>> lines = ['objID                   & osrcid            & xsrcid       ',
   ...          '----------------------- & ----------------- & -------------',
   ...          '              277955213 & S000.7044P00.7513 & XS04861B6_005',
   ...          '              889974380 & S002.9051P14.7003 & XS03957B7_004']
   >>> data = ascii.read(lines, data_start=2, delimiter='&')
   >>> print(data)
     objID         osrcid          xsrcid
   --------- ----------------- -------------
   277955213 S000.7044P00.7513 XS04861B6_005
   889974380 S002.9051P14.7003 XS03957B7_004

If the format of a file is known (e.g., it is a fixed-width table or an IPAC
table), then it is more efficient and reliable to provide a value for the
``format`` argument from one of the values in the `supported formats`_. For
example::

   >>> data = ascii.read(lines, format='fixed_width_two_line', delimiter='&')

For simpler formats such as CSV, |read| will automatically try reading with the
Cython/C parsing engine, which is significantly faster than the ordinary Python
implementation (described in :ref:`fast_ascii_io`). If the fast engine fails,
|read| will fall back on the Python reader by default. The argument
``fast_reader`` can be specified to control this behavior. For example, to
disable the fast engine::

   >>> data = ascii.read(lines, format='csv', fast_reader=False)

For reading very large tables see the section on :ref:`chunk_reading` or
use `pandas <https://pandas.pydata.org/>`_ (see Note below).

.. Note::

   Reading a table which contains unicode characters is supported with the
   pure Python readers by specifying the ``encoding`` parameter. The fast
   C-readers do not support unicode. For large data files containing unicode,
   we recommend reading the file using `pandas <https://pandas.pydata.org/>`_
   and converting to a :ref:`Table <astropy-table>` via the :ref:`Table -
   Pandas interface <pandas>`.

..
  EXAMPLE END

Writing Tables
--------------

The |write| function provides a way to write a data table as a formatted ASCII
table.

Examples
^^^^^^^^

..
  EXAMPLE START
  Writing Data Tables as Formatted ASCII Tables

The following writes a table as a simple space-delimited file::

  >>> import numpy as np
  >>> from astropy.table import Table, Column, MaskedColumn
  >>> x = np.array([1, 2, 3])
  >>> y = x ** 2
  >>> data = Table([x, y], names=['x', 'y'])
  >>> ascii.write(data, 'values.dat', overwrite=True)

The ``values.dat`` file will then contain::

  x y
  1 1
  2 4
  3 9

Most of the input Reader formats supported by `astropy.io.ascii` for reading are
also supported for writing. This provides a great deal of flexibility in the
format for writing. The example below writes the data as a LaTeX table, using
the option to send the output to ``sys.stdout`` instead of a file::

  >>> import sys
  >>> ascii.write(data, sys.stdout, format='latex')
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

Finally, one can write data in the :ref:`ecsv_format` which allows for
preserving table metadata such as column data types and units. In this way, a
data table (including one with masked entries) can be stored and read back as
ASCII with no loss of information.

  >>> t = Table(masked=True)
  >>> t['x'] = MaskedColumn([1.0, 2.0], unit='m', dtype='float32')
  >>> t['x'][1] = np.ma.masked
  >>> t['y'] = MaskedColumn([False, True], dtype='bool')

  >>> import io
  >>> fh = io.StringIO()
  >>> t.write(fh, format='ascii.ecsv')  # doctest: +SKIP
  >>> table_string = fh.getvalue()      # doctest: +SKIP
  >>> print(table_string)               # doctest: +SKIP
  # %ECSV 0.9
  # ---
  # datatype:
  # - {name: x, unit: m, datatype: float32}
  # - {name: y, datatype: bool}
  x y
  1.0 False
  "" True

  >>> Table.read(table_string, format='ascii')  # doctest: +SKIP
  <Table masked=True length=2>
     x      y
     m
  float32  bool
  ------- -----
      1.0 False
       --  True

.. Note::

   For most supported formats one can write a masked table and then
   read it back without losing information about the masked table
   entries. This is accomplished by using a blank string entry to
   indicate a masked (missing) value. See the :ref:`replace_bad_or_missing_values`
   section for more information.

..
  EXAMPLE END

.. _supported_formats:

Supported Formats
=================

A full list of the supported ``format`` values and corresponding format types
for ASCII tables is given below. The ``Write`` column indicates which formats
support write functionality, and the ``Fast`` column indicates which formats
are compatible with the fast Cython/C engine for reading and writing.

========================= ===== ==== ============================================================================================
           Format         Write Fast                                          Description
========================= ===== ==== ============================================================================================
``aastex``                  Yes      :class:`~astropy.io.ascii.AASTex`: AASTeX deluxetable used for AAS journals
``basic``                   Yes  Yes :class:`~astropy.io.ascii.Basic`: Basic table with custom delimiters
``cds``                              :class:`~astropy.io.ascii.Cds`: CDS format table
``commented_header``        Yes  Yes :class:`~astropy.io.ascii.CommentedHeader`: Column names in a commented line
``csv``                     Yes  Yes :class:`~astropy.io.ascii.Csv`: Basic table with comma-separated values
``daophot``                          :class:`~astropy.io.ascii.Daophot`: IRAF DAOphot format table
``ecsv``                    Yes      :class:`~astropy.io.ascii.Ecsv`: Enhanced CSV format
``fixed_width``             Yes      :class:`~astropy.io.ascii.FixedWidth`: Fixed width
``fixed_width_no_header``   Yes      :class:`~astropy.io.ascii.FixedWidthNoHeader`: Fixed-width with no header
``fixed_width_two_line``    Yes      :class:`~astropy.io.ascii.FixedWidthTwoLine`: Fixed-width with second header line
``html``                    Yes      :class:`~astropy.io.ascii.HTML`: HTML format table
``ipac``                    Yes      :class:`~astropy.io.ascii.Ipac`: IPAC format table
``latex``                   Yes      :class:`~astropy.io.ascii.Latex`: LaTeX table
``no_header``               Yes  Yes :class:`~astropy.io.ascii.NoHeader`: Basic table with no headers
``rdb``                     Yes  Yes :class:`~astropy.io.ascii.Rdb`: Tab-separated with a type definition header line
``rst``                     Yes      :class:`~astropy.io.ascii.RST`: reStructuredText simple format table
``sextractor``                       :class:`~astropy.io.ascii.SExtractor`: SExtractor format table
``tab``                     Yes  Yes :class:`~astropy.io.ascii.Tab`: Basic table with tab-separated values
========================= ===== ==== ============================================================================================


.. attention:: **ECSV is recommended**

   For writing and reading tables to ASCII in a way that fully reproduces the
   table data, types and metadata (i.e., the table will "round-trip"), we highly
   recommend using the :ref:`ecsv_format`. This writes the actual data in a
   simple space-delimited format (the ``basic`` format) that any ASCII table
   reader can parse, but also includes metadata encoded in a comment block that
   allows full reconstruction of the original columns. This includes support
   for :ref:`ecsv_format_mixin_columns` (such as
   `~astropy.coordinates.SkyCoord` or `~astropy.time.Time`) and
   :ref:`ecsv_format_masked_columns`.


Using `astropy.io.ascii`
========================

The details of using `astropy.io.ascii` are provided in the following sections:

Reading tables
---------------

.. toctree::
   :maxdepth: 2

   read

Writing tables
---------------

.. toctree::
   :maxdepth: 2

   write

Fixed-Width Gallery
--------------------

.. toctree::
   :maxdepth: 2

   fixed_width_gallery

Fast ASCII Engine
-----------------

.. toctree::
   :maxdepth: 2

   fast_ascii_io

Base Class Elements
-------------------

.. toctree::
   :maxdepth: 2

   base_classes

Extension Reader Classes
------------------------

.. toctree::
   :maxdepth: 2

   extension_classes

.. note that if this section gets too long, it should be moved to a separate
   doc page - see the top of performance.inc.rst for the instructions on how to do
   that
.. include:: performance.inc.rst

Reference/API
=============

.. automodapi:: astropy.io.ascii
