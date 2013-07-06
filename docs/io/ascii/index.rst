.. include:: references.txt

.. _io-ascii:

*********************************
ASCII Tables (`astropy.io.ascii`)
*********************************

Introduction
============

`astropy.io.ascii` provides methods for reading and writing a wide range of ASCII data table
formats via built-in :ref:`extension_reader_classes`.  The emphasis is on flexibility and ease of use.

The following shows a few of the ASCII formats that are available, while the section on
`Supported formats`_ contains the full list.

* :class:`~astropy.io.ascii.basic.Basic`: basic table with customizable delimiters and header configurations
* :class:`~astropy.io.ascii.cds.Cds`: `CDS format table <http://vizier.u-strasbg.fr/doc/catstd.htx>`_ (also Vizier and ApJ machine readable tables)
* :class:`~astropy.io.ascii.daophot.Daophot`: table from the IRAF DAOphot package
* :class:`~astropy.io.ascii.fixedwidth.FixedWidth`: table with fixed-width columns (see also :ref:`fixed_width_gallery`)
* :class:`~astropy.io.ascii.ipac.Ipac`: `IPAC format table <http://irsa.ipac.caltech.edu/applications/DDGEN/Doc/ipac_tbl.html>`_
* :class:`~astropy.io.ascii.latex.Latex`: LaTeX table with datavalue in the `tabular` environment
* :class:`~astropy.io.ascii.basic.Rdb`: tab-separated values with an extra line after the column definition line
* :class:`~astropy.io.ascii.sextractor.SExtractor`: `SExtractor format table <http://www.astromatic.net/software/sextractor>`_

The :mod:`astropy.io.ascii` package is built on a modular and extensible class
structure with independent :ref:`base_class_elements` so that new formats can
be easily accomodated.

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
  >>> data = ascii.read("sources.dat")
  >>> print data
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
                '----------------------- & ----------------- & -------------',
                '              277955213 & S000.7044P00.7513 & XS04861B6_005',
                '              889974380 & S002.9051P14.7003 & XS03957B7_004']
   >>> data = ascii.read(lines, data_start=2, delimiter='&')
   >>> print data
     objID         osrcid          xsrcid
   --------- ----------------- -------------
   277955213 S000.7044P00.7513 XS04861B6_005
   889974380 S002.9051P14.7003 XS03957B7_004

If the format of a file is known (e.g. it is a fixed width table or an IPAC table),
then it is more efficient and reliable to provide a value for the ``format`` argument from one
of the values in the `supported formats`_.  For example::

   >>> data = ascii.read(lines, format='fixed_width_two_line', delimiter='&')

Writing Tables
--------------

The |write| function provides a way to write a data table as a formatted ASCII
table.  For example the following writes a table as a simple space-delimited
file::

  >>> x = np.array([1, 2, 3])
  >>> y = x ** 2
  >>> data = Table([x, y], names=['x', 'y'])
  >>> ascii.write(data, 'values.dat')

The ``values.dat`` file will then contain::

  x y
  1 1
  2 4
  3 9

All of the input Reader formats supported by `astropy.io.ascii` for reading are
also supported for writing.  This provides a great deal of flexibility in the
format for writing.  The example below writes the data as a LaTeX table, using
the option to send the output to ``sys.stdout`` instead of a file::

  >>> ascii.write(data, sys.stdout, format='latex')
  \begin{table}
  \begin{tabular}{cc}
  x & y \\
  1 & 1 \\
  2 & 4 \\
  3 & 9 \\
  \end{tabular}
  \end{table}

.. _supported_formats:

Supported formats
=================

A full list of the supported ``format`` values and corresponding format types for ASCII
tables is given below.  The ``Write`` column indicates which formats support write
functionality.

========================= ===== ============================================================================================
           Format               Write                                          Description
========================= ===== ============================================================================================
``aastex``                  Yes :class:`~astropy.io.ascii.latex.AASTex`: AASTeX deluxetable used for AAS journals
``basic``                   Yes :class:`~astropy.io.ascii.basic.Basic`: Basic table with custom delimiters
``cds``                         :class:`~astropy.io.ascii.cds.Cds`: CDS format table
``commented_header``        Yes :class:`~astropy.io.ascii.basic.CommentedHeader`: Column names in a commented line
``daophot``                     :class:`~astropy.io.ascii.daophot.Daophot`: IRAF DAOphot format table
``fixed_width``             Yes :class:`~astropy.io.ascii.fixedwidth.FixedWidth`: Fixed width
``fixed_width_no_header``   Yes :class:`~astropy.io.ascii.fixedwidth.FixedWidthNoHeader`: Fixed width with no header
``fixed_width_two_line``    Yes :class:`~astropy.io.ascii.fixedwidth.FixedWidthTwoLine`: Fixed width with second header line
``ipac``                    Yes :class:`~astropy.io.ascii.ipac.Ipac`: IPAC format table
``latex``                   Yes :class:`~astropy.io.ascii.latex.Latex`: LaTeX table
``no_header``               Yes :class:`~astropy.io.ascii.basic.NoHeader`: Basic table with no headers
``rdb``                     Yes :class:`~astropy.io.ascii.basic.Rdb`: Tab-separated with a type definition header line
``sextractor``                  :class:`~astropy.io.ascii.sextractor.SExtractor`: SExtractor format table
``tab``                     Yes :class:`~astropy.io.ascii.basic.Tab`: Basic table with tab-separated values
========================= ===== ============================================================================================


Using `io.ascii`
================

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


