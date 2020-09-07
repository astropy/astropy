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
* :class:`~astropy.io.ascii.NoHeader`: basic table with no header where columns are auto-named.
* :class:`~astropy.io.ascii.Rdb`: tab-separated values with an extra line after the column definition line.
* :class:`~astropy.io.ascii.RST`: `reStructuredText simple format table <https://docutils.sourceforge.io/docs/ref/rst/restructuredtext.html#simple-tables>`_.
* :class:`~astropy.io.ascii.SExtractor`: `SExtractor format table <http://www.astromatic.net/software/sextractor>`_.
* :class:`~astropy.io.ascii.Tab`: tab-separated values.
