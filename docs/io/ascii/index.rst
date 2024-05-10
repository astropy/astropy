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
improved performance. This subpackage was originally developed as ``asciitable``.

The following shows a few of the ASCII formats that are available, while the
section on `Supported formats`_ contains the full list.

* :class:`~astropy.io.ascii.Basic`: basic table with customizable delimiters and header configurations
* :class:`~astropy.io.ascii.Cds`: `CDS format table <https://vizier.unistra.fr/doc/catstd.htx>`_ (also Vizier)
* :class:`~astropy.io.ascii.Daophot`: table from the IRAF DAOphot package
* :class:`~astropy.io.ascii.Ecsv`: :ref:`ecsv_format` for lossless round-trip of data tables (**recommended**)
* :class:`~astropy.io.ascii.FixedWidth`: table with fixed-width columns (see also :ref:`fixed_width_gallery`)
* :class:`~astropy.io.ascii.Ipac`: `IPAC format table <https://irsa.ipac.caltech.edu/applications/DDGEN/Doc/ipac_tbl.html>`_
* :class:`~astropy.io.ascii.HTML`: HTML format table contained in a <table> tag
* :class:`~astropy.io.ascii.Latex`: LaTeX table with datavalue in the ``tabular`` environment
* :class:`~astropy.io.ascii.Mrt`: AAS `Machine-Readable Tables (MRT) <https://journals.aas.org/mrt-standards/>`_)
* :class:`~astropy.io.ascii.SExtractor`: `SExtractor format table <https://sextractor.readthedocs.io/en/latest/>`_

The strength of `astropy.io.ascii` is the support for astronomy-specific
formats (often with metadata) and specialized data types such as
:ref:`SkyCoord <astropy-coordinates-high-level>`, :ref:`Time
<astropy-time>`, and :ref:`Quantity <quantity>`. For reading or writing large
data tables in a generic format such as CSV, using the :ref:`Table - Pandas
interface <pandas>` is an option to consider.

.. note::

    It is strongly encouraged to use the functionality from
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

To specify specific data types for one or more columns, use the ``converters``
argument (see :ref:`io-ascii-read-converters` for details). For instance if the
``obsid`` is actually a string identifier (instead of an integer) you can read
the table with the code below. This also illustrates using the preferred
:ref:`Table interface <table_io>` for reading::

  >>> from astropy.table import Table
  >>> sources = """
  ... target observatory obsid
  ... TW_Hya Chandra     22178
  ... MP_Mus XMM         0406030101"""
  >>> data = Table.read(sources, format='ascii', converters={'obsid': str})
  >>> data
  <Table length=2>
  target observatory   obsid
   str6      str7      str10
  ------ ----------- ----------
  TW_Hya     Chandra      22178
  MP_Mus         XMM 0406030101

Writing Tables
--------------

The |write| function provides a way to write a data table as a formatted ASCII
table.  Most of the input table :ref:`supported_formats` for reading are also
available for writing. This provides a great deal of flexibility in the format
for writing.

..
  EXAMPLE START
  Writing Data Tables as Formatted ASCII Tables

The following shows how to write a formatted ASCII table using the |write|
function::

  >>> import numpy as np
  >>> from astropy.io import ascii
  >>> from astropy.table import Table
  >>> data = Table()
  >>> data['x'] = np.array([1, 2, 3], dtype=np.int32)
  >>> data['y'] = data['x'] ** 2
  >>> ascii.write(data, 'values.dat', overwrite=True)

The ``values.dat`` file will then contain::

  x y
  1 1
  2 4
  3 9

It is also possible and encouraged to use the write functionality from
:mod:`astropy.io.ascii` through a higher level interface in the :ref:`Data
Tables <astropy-table>` package (see :ref:`table_io` for more details). For
example::

  >>> data.write('values.dat', format='ascii', overwrite=True)

.. attention:: **ECSV is recommended**

   For a reproducible ASCII version of your table, we recommend using the
   :ref:`ecsv_format`. This stores all the table meta-data (in particular the
   column types and units) to a comment section at the beginning while
   maintaining compatibility with most plain CSV readers. It also allows storing
   richer data like `~astropy.coordinates.SkyCoord` or multidimensional or
   variable-length columns. ECSV is also supported in Java by |STIL| and
   |TOPCAT| (see :ref:`ecsv_format`).

To write our simple example table to ECSV we use::

  >>> data.write('values.ecsv', overwrite=True)  # doctest: +SKIP

The ``.ecsv`` extension is recognized and implies using ECSV (equivalent to
``format='ascii.ecsv'``). The ``values.ecsv`` file will then contain::

  # %ECSV 1.0
  # ---
  # datatype:
  # - {name: x, datatype: int32}
  # - {name: y, datatype: int32}
  # schema: astropy-2.0
  x y
  1 1
  2 4
  3 9

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
``cds``                     Yes      :class:`~astropy.io.ascii.Cds`: CDS format table
``commented_header``        Yes  Yes :class:`~astropy.io.ascii.CommentedHeader`: Column names in a commented line
``csv``                     Yes  Yes :class:`~astropy.io.ascii.Csv`: Basic table with comma-separated values
``daophot``                          :class:`~astropy.io.ascii.Daophot`: IRAF DAOphot format table
``ecsv``                    Yes      :class:`~astropy.io.ascii.Ecsv`: Enhanced CSV format (**recommended**)
``fixed_width``             Yes      :class:`~astropy.io.ascii.FixedWidth`: Fixed width
``fixed_width_no_header``   Yes      :class:`~astropy.io.ascii.FixedWidthNoHeader`: Fixed-width with no header
``fixed_width_two_line``    Yes      :class:`~astropy.io.ascii.FixedWidthTwoLine`: Fixed-width with second header line
``html``                    Yes      :class:`~astropy.io.ascii.HTML`: HTML format table
``ipac``                    Yes      :class:`~astropy.io.ascii.Ipac`: IPAC format table
``latex``                   Yes      :class:`~astropy.io.ascii.Latex`: LaTeX table
``mrt``                     Yes      :class:`~astropy.io.ascii.Mrt`: AAS Machine-Readable Table format
``no_header``               Yes  Yes :class:`~astropy.io.ascii.NoHeader`: Basic table with no headers
``qdp``                     Yes      :class:`~astropy.io.ascii.QDP`: Quick and Dandy Plotter files
``rdb``                     Yes  Yes :class:`~astropy.io.ascii.Rdb`: Tab-separated with a type definition header line
``rst``                     Yes      :class:`~astropy.io.ascii.RST`: reStructuredText simple format table
``sextractor``                       :class:`~astropy.io.ascii.SExtractor`: SExtractor format table
``tab``                     Yes  Yes :class:`~astropy.io.ascii.Tab`: Basic table with tab-separated values
========================= ===== ==== ============================================================================================

Getting Help
============

Some formats have additional options that can be set to control the behavior of the
reader or writer. For more information on these options, you can either see the
documentation for the specific format class (e.g. :class:`~astropy.io.ascii.HTML`) or
use the ``help`` function of the ``read`` or ``write`` functions. For example:

.. doctest-skip::

  >>> ascii.read.help()  # Common help for all formats
  >>> ascii.read.help("html")  # Common help plus "html" format-specific args
  >>> ascii.write.help("latex")  # Common help plus "html" format-specific args

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

ECSV Format
-----------

.. toctree::
   :maxdepth: 2

   ecsv

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

.. toctree::
   :maxdepth: 2

   ref_api
