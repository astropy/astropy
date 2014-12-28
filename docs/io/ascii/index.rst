.. include:: references.txt

.. _io-ascii:

*********************************
ASCII Tables (`astropy.io.ascii`)
*********************************

Introduction
============

`astropy.io.ascii` provides methods for reading and writing a wide range of ASCII data table
formats via built-in :ref:`extension_reader_classes`.  The emphasis is on flexibility and ease of use,
although readers can optionally use a less flexible C/Cython engine for reading and writing for
improved performance.

The following shows a few of the ASCII formats that are available, while the section on
`Supported formats`_ contains the full list.

* :class:`~astropy.io.ascii.Basic`: basic table with customizable delimiters and header configurations
* :class:`~astropy.io.ascii.Cds`: `CDS format table <http://vizier.u-strasbg.fr/doc/catstd.htx>`_ (also Vizier and ApJ machine readable tables)
* :class:`~astropy.io.ascii.Daophot`: table from the IRAF DAOphot package
* :class:`~astropy.io.ascii.Ecsv`: Enhanced CSV format
* :class:`~astropy.io.ascii.FixedWidth`: table with fixed-width columns (see also :ref:`fixed_width_gallery`)
* :class:`~astropy.io.ascii.Ipac`: `IPAC format table <http://irsa.ipac.caltech.edu/applications/DDGEN/Doc/ipac_tbl.html>`_
* :class:`~astropy.io.ascii.HTML`: HTML format table contained in a <table> tag
* :class:`~astropy.io.ascii.Latex`: LaTeX table with datavalue in the ``tabular`` environment
* :class:`~astropy.io.ascii.Rdb`: tab-separated values with an extra line after the column definition line
* :class:`~astropy.io.ascii.SExtractor`: `SExtractor format table <http://www.astromatic.net/software/sextractor>`_

The :mod:`astropy.io.ascii` package is built on a modular and extensible class
structure with independent :ref:`base_class_elements` so that new formats can
be easily accommodated.

.. note::

    It is also possible to use the functionality from
    :mod:`astropy.io.ascii` through a higher-level interface in the
    :mod:`astropy.table` package. See :ref:`table_io` for more details.

Getting Started
===============

Reading Tables
--------------

The majority of commonly encountered ASCII tables can be easily read with the |read|
function.  Assume you have a file named ``sources.dat`` with the following contents::

  obsid redshift  X      Y     object
  3102  0.32      4167  4085   Q1250+568-A
  877   0.22      4378  3892   "Source 82"

This table can be read with the following::

  >>> from astropy.io import ascii
  >>> data = ascii.read("sources.dat")  # doctest: +SKIP
  >>> print data  # doctest: +SKIP
  obsid redshift  X    Y      object
  ----- -------- ---- ---- -----------
   3102     0.32 4167 4085 Q1250+568-A
    877     0.22 4378 3892   Source 82

The first argument to the |read| function can be the name of a file, a string
representation of a table, or a list of table lines.  By default |read| will
try to `guess the table format <#guess-table-format>`_ by trying all the
`supported formats`_.  If this does not work (for unusually formatted tables) then
one needs give astropy.io.ascii additional hints about the format, for
example::

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

If the format of a file is known (e.g. it is a fixed width table or an IPAC table),
then it is more efficient and reliable to provide a value for the ``format`` argument from one
of the values in the `supported formats`_.  For example::

   >>> data = ascii.read(lines, format='fixed_width_two_line', delimiter='&')

For simpler formats such as CSV, |read| will automatically try reading with the
Cython/C parsing engine, which is significantly faster than the ordinary Python
implementation (described in :ref:`fast_ascii_io`). If the fast engine fails,
|read| will fall back on the Python reader by default. The argument
``fast_reader`` can be specified to control this behavior. For example, to
disable the fast engine::

   >>> data = ascii.read(lines, format='csv', fast_reader=False)

Writing Tables
--------------

The |write| function provides a way to write a data table as a formatted ASCII
table.  For example the following writes a table as a simple space-delimited
file::

  >>> import numpy as np
  >>> from astropy.table import Table, Column
  >>> x = np.array([1, 2, 3])
  >>> y = x ** 2
  >>> data = Table([x, y], names=['x', 'y'])
  >>> ascii.write(data, 'values.dat')

The ``values.dat`` file will then contain::

  x y
  1 1
  2 4
  3 9

Most of the input Reader formats supported by `astropy.io.ascii` for reading are
also supported for writing.  This provides a great deal of flexibility in the
format for writing.  The example below writes the data as a LaTeX table, using
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

Finally, one can write data in the ECSV table format which allows preserving
table meta-data such as column data types and units.  In this way a data table
can be stored and read back as ASCII with no loss of information.

  >>> t = Table()
  >>> t['x'] = Column([1.0, 2.0], unit='m', dtype='float32')
  >>> t['y'] = Column([False, True], dtype='bool')

  >>> from StringIO import StringIO
  >>> fh = StringIO()
  >>> t.write(fh, format='ascii.ecsv')
  >>> table_string = fh.getvalue()
  >>> print(table_string)
  # %ECSV 1.0
  # ---
  # columns:
  # - {name: x, unit: m, type: float32}
  # - {name: y, type: bool}
  x y
  1.0 False
  2.0 True

  >>> Table.read(table_string, format='ascii')
  <Table rows=2 names=('x','y') units=('m',None)>
  array([(1.0, False), (2.0, True)],
        dtype=[('x', '<f4'), ('y', '?')])


.. _supported_formats:

Supported formats
=================

A full list of the supported ``format`` values and corresponding format types for ASCII
tables is given below.  The ``Write`` column indicates which formats support write
functionality, and the ``Fast`` column indicates which formats are compatible with
the fast Cython/C engine for reading and writing.

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
``fixed_width_no_header``   Yes      :class:`~astropy.io.ascii.FixedWidthNoHeader`: Fixed width with no header
``fixed_width_two_line``    Yes      :class:`~astropy.io.ascii.FixedWidthTwoLine`: Fixed width with second header line
``html``                    Yes      :class:`~astropy.io.ascii.HTML`: HTML format table
``ipac``                    Yes      :class:`~astropy.io.ascii.Ipac`: IPAC format table
``latex``                   Yes      :class:`~astropy.io.ascii.Latex`: LaTeX table
``no_header``               Yes  Yes :class:`~astropy.io.ascii.NoHeader`: Basic table with no headers
``rdb``                     Yes  Yes :class:`~astropy.io.ascii.Rdb`: Tab-separated with a type definition header line
``sextractor``                       :class:`~astropy.io.ascii.SExtractor`: SExtractor format table
``tab``                     Yes  Yes :class:`~astropy.io.ascii.Tab`: Basic table with tab-separated values
========================= ===== ==== ============================================================================================


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

Fixed-width Gallery
--------------------

.. toctree::
   :maxdepth: 2

   fixed_width_gallery

Fast ASCII Engine
-----------------

.. toctree::
   :maxdepth: 2

   fast_ascii_io

Base class elements
-------------------

.. toctree::
   :maxdepth: 2

   base_classes

Extension Reader classes
------------------------

.. toctree::
   :maxdepth: 2

   extension_classes


Reference/API
=============

.. automodapi:: astropy.io.ascii
