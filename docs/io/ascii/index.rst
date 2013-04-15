.. include:: references.txt

.. _io-ascii:

*********************************
ASCII Tables (`astropy.io.ascii`)
*********************************

Introduction
============

`astropy.io.ascii` provides methods for reading and writing a wide range of ASCII data table
formats via built-in :ref:`extension_reader_classes`.  The emphasis is on flexibility and ease of use.

The following formats are supported:

* :class:`~astropy.io.ascii.latex.AASTex`: AASTeX `deluxetable` used for AAS journals
* :class:`~astropy.io.ascii.basic.Basic`: basic table with customizable delimiters and header configurations
* :class:`~astropy.io.ascii.cds.Cds`: `CDS format table <http://vizier.u-strasbg.fr/doc/catstd.htx>`_ (also Vizier and ApJ machine readable tables)
* :class:`~astropy.io.ascii.basic.CommentedHeader`: column names given in a line that begins with the comment character
* :class:`~astropy.io.ascii.daophot.Daophot`: table from the IRAF DAOphot package
* :class:`~astropy.io.ascii.fixedwidth.FixedWidth`: table with fixed-width columns (see also :ref:`fixed_width_gallery`)
* :class:`~astropy.io.ascii.ipac.Ipac`: `IPAC format table <http://irsa.ipac.caltech.edu/applications/DDGEN/Doc/ipac_tbl.html>`_
* :class:`~astropy.io.ascii.latex.Latex`: LaTeX table with datavalue in the `tabular` environment
* :class:`~astropy.io.ascii.basic.NoHeader`: basic table with no header where columns are auto-named
* :class:`~astropy.io.ascii.basic.Rdb`: tab-separated values with an extra line after the column definition line
* :class:`~astropy.io.ascii.sextractor.SExtractor`: `SExtractor format table <http://www.astromatic.net/software/sextractor>`_
* :class:`~astropy.io.ascii.basic.Tab`: tab-separated values

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
supported formats.  If this does not work (for unusually formatted tables) then
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

  >>> ascii.write(data, sys.stdout, Writer=ascii.Latex)
  \begin{table}
  \begin{tabular}{cc}
  x & y \\
  1 & 1 \\
  2 & 4 \\
  3 & 9 \\
  \end{tabular}
  \end{table}

Using `io.ascii`
================

The details of using `astropy.io.ascii` are provided in the following sections:

.. toctree::
   :maxdepth: 2

   read
   write
   fixed_width_gallery
   base_classes
   extension_classes


Reference/API
=============

.. automodapi:: astropy.io.ascii


