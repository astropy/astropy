.. include:: references.txt

.. _extension_reader_classes:

Extension Reader Classes
************************

The following classes extend the base :class:`~astropy.io.ascii.BaseReader` functionality to handle reading and writing
different table formats. Some, such as the :class:`~astropy.io.ascii.Basic` Reader class
are fairly general and include a number of configurable attributes. Others
such as :class:`~astropy.io.ascii.Cds` or :class:`~astropy.io.ascii.Daophot` are specialized to read certain
well-defined but idiosyncratic formats.

* :class:`~astropy.io.ascii.AASTex`: AASTeX `deluxetable <https://fits.gsfc.nasa.gov/standard30/deluxetable.sty>`_ used for AAS journals.
* :class:`~astropy.io.ascii.Basic`: basic table with customizable delimiters and header configurations.
* :class:`~astropy.io.ascii.Cds`: `CDS format table <http://vizier.u-strasbg.fr/doc/catstd.htx>`_ (also Vizier and ApJ machine readable tables).
* :class:`~astropy.io.ascii.CommentedHeader`: column names given in a line that begins with the comment character.
* :class:`~astropy.io.ascii.Csv`: comma-separated values.
* :class:`~astropy.io.ascii.Daophot`: table from the IRAF DAOphot package.
* :class:`~astropy.io.ascii.FixedWidth`: table with fixed-width columns (see also :ref:`fixed_width_gallery`).
* :class:`~astropy.io.ascii.FixedWidthNoHeader`: table with fixed-width columns and no header.
* :class:`~astropy.io.ascii.FixedWidthTwoLine`: table with fixed-width columns and a two-line header.
* :class:`~astropy.io.ascii.HTML`: HTML format table contained in a <table> tag.
* :class:`~astropy.io.ascii.Ipac`: `IPAC format table <https://irsa.ipac.caltech.edu/applications/DDGEN/Doc/ipac_tbl.html>`_.
* :class:`~astropy.io.ascii.Latex`: LaTeX table with datavalue in the ``tabular`` environment.
* :class:`~astropy.io.ascii.Mrt`: `AAS Machine-Readable Table format <https://journals.aas.org/mrt-standards/>`_.
* :class:`~astropy.io.ascii.NoHeader`: basic table with no header where columns are auto-named.
* :class:`~astropy.io.ascii.Rdb`: tab-separated values with an extra line after the column definition line.
* :class:`~astropy.io.ascii.RST`: `reStructuredText simple format table <https://docutils.sourceforge.io/docs/ref/rst/restructuredtext.html#simple-tables>`_.
* :class:`~astropy.io.ascii.SExtractor`: `SExtractor format table <https://sextractor.readthedocs.io/en/latest/>`_.
* :class:`~astropy.io.ascii.Tab`: tab-separated values.


A gallery of complicated cases
==============================

The purpose of this section is to provide a few example how we can
deal with tables with a format that fails the default readers or that
cannot be written with the default writers from the list above. We
walk through a few examples to show what steps we can take to solve
our table/reading writing problem, but of course this gallery does not
include every possible way in which a table might fail to read or
write.

Making a table easier to read
-----------------------------

Sometimes the Reader/Writer classes are not enough, even with all the
flexibility to customize, the delimiters, column names, format
converters and other properties. To read just a single table that has
a format close to, but not identical with, any of the formats in the
list above, the fastest solution is often to open that one table file
in a text editor to modify it until it does conform to a format that
any of the readers above can already read.

Badly formatted header line
^^^^^^^^^^^^^^^^^^^^^^^^^^^
The following table will fail to parse (raising an
`~astropy.io.InconsistentTableError') because the header line looks as
if there were three columns, while in fact, there are only two::

  Name spectral type
  Vega A0
  Altair A7

Opening this file in a text editor to fix the format is easy::

  Name "spectral type"
  Vega A0
  Altair A7

or
::

  Name spectral_type
  Vega A0
  Altair A7

would both work, e.g.

..
  EXAMPLE START
  Make a table easier to read

::

  >>> from astropy.io import ascii
  >>> table = """
  ... Star "spectral type"
  ... Vega A0
  ... Altair A7
  ... """
  >>> ascii.read(table)
  <Table length=2>
   Star  spectral type
   str6       str2
  ------ -------------
    Vega            A0
  Altair            A7

Similarly, the header could be commented out and the names of the table column
can be set manually::

  >>> from astropy.io import ascii
  >>> table = """
  ... #Star spectral type
  ... Vega A0
  ... Altair A7
  ... """
  >>> ascii.read(table, names=["Star", "spectral type"], format='no_header')
  <Table length=2>
   Star  spectral type
   str6       str2
  ------ -------------
    Vega            A0
  Altair            A7

This last experiment might show us that reading works just fine, if we just
find a way to ignore the badly formatted header line. That can be done witout
any modification of the table itself, but just using the ``data_start``
parameter::

   >>> table = """
   ... Star spectral type
   ... Vega A0
   ... Altair A7
   ... """
   >>> ascii.read(table, names=["Star", "spectral type"], data_start=1)
   <Table length=2>
    Star  spectral type
    str6       str2
   ------ -------------
     Vega            A0
   Altair            A7

..
  EXAMPLE END
