.. include:: references.txt

ASCII Tables (`astropy.io.ascii`)
===================================

Introduction
--------------

`astropy.io.ascii` provides methods for reading and writing a wide range of ASCII data table
formats via built-in :ref:`extension_reader_classes`.  The ephemsis is on flexibility and ease of use.

The following formats are supported:

* :class:`~astropy.io.ascii.Basic`: basic table with customizable delimiters and header configurations
* :class:`~astropy.io.ascii.Cds`: `CDS format table <http://vizier.u-strasbg.fr/doc/catstd.htx>`_ (also Vizier and ApJ machine readable tables)
* :class:`~astropy.io.ascii.CommentedHeader`: column names given in a line that begins with the comment character
* :class:`~astropy.io.ascii.Daophot`: table from the IRAF DAOphot package
* :class:`~astropy.io.ascii.FixedWidth`: table with fixed-width columns (:ref:`fixed_width_gallery`)
* :class:`~astropy.io.ascii.Ipac`: `IPAC format table <http://irsa.ipac.caltech.edu/applications/DDGEN/Doc/ipac_tbl.html>`_
* :class:`~astropy.io.ascii.Latex`, :class:`~astropy.io.ascii.AASTex`: LaTeX tables (plain and AASTex)
* :class:`~astropy.io.ascii.Memory`: table already in memory (list of lists, dict of lists, etc)
* :class:`~astropy.io.ascii.NoHeader`: basic table with no header where columns are auto-named
* :class:`~astropy.io.ascii.Rdb`: tab-separated values with an extra line after the column definition line
* :class:`~astropy.io.ascii.Tab`: tab-separated values

The :mod:`astropy.io.ascii` package is built on a modular and extensible class
structure with independent :ref:`base_class_elements` so that new formats can
be easily accomodated.

Getting Started
------------------

Reading Tables
^^^^^^^^^^^^^^^

The majority of commonly encountered ASCII tables can be easily read with the |read|
function.  Assume you have a file named ``sources.dat`` with the following contents::

  obsid redshift  X      Y     object     
  3102  0.32      4167  4085   Q1250+568-A
  877   0.22      4378  3892   "Source 82"

This table can be read with the following::

  import astropy.io.ascii as ascii
  data = ascii.read("sources.dat")
  print data

The first argument to the |read| function can be the name of a file, a string
representation of a table, or a list of table lines.  By default |read| will
try to `guess the table format <#guess-table-format>`_ by trying all the
supported formats.  If this does not work (for unusually formatted tables) then
one needs give astropy.io.ascii additional hints about the format, for
example::

   lines = ['objID                   & osrcid            & xsrcid       ',  
            '----------------------- & ----------------- & -------------',
            '              277955213 & S000.7044P00.7513 & XS04861B6_005',
            '              889974380 & S002.9051P14.7003 & XS03957B7_004']  
   data = ascii.read(lines, data_start=2, delimiter='&')
   print data


Writing Tables
^^^^^^^^^^^^^^^

The |write| function provides a way to write a data table as a formatted ASCII table.  For example::

  >>> from astropy.table import Table
  >>> x = np.array([1, 2, 3])
  >>> y = x ** 2
  >>> data = Table([x, y], names=['x', 'y'])
  >>> astropy.io.ascii.write(data, 'values.dat')

The ``values.dat`` file will then contain::

  x y
  1 1
  2 4
  3 9

All of the input (read) formats supported by `astropy.io.ascii` are also
supported for writing.  This provides a great deal of flexibility in the format
for writing.  The example below writes the data as a LaTeX table, using the
option to send the output to ``sys.stdout`` instead of a file::

  >>> ascii.write(data, sys.stdout, Writer=ascii.Latex)
  \begin{table}
  \begin{tabular}{cc}
  x & y \\
  1 & 1 \\
  2 & 4 \\
  3 & 9 \\
  \end{tabular}
  \end{table}

Using `astropy.io.ascii`
-------------------------

The details of using `astropy.io.ascii` are provided in the following sections:

.. toctree::
   :maxdepth: 2

   read
   write
   fixed_width_gallery
   base_classes
   extension_classes

.. automodapi:: astropy.io.ascii


