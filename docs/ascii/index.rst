.. include:: references.txt

================
astropy.io.ascii
================
An extensible ASCII table reader and writer for Python 2 and 3.

The :mod:`astropy.io.ascii` packge can read and write a wide range of ASCII table
formats via built-in :ref:`extension_reader_classes`:

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

Contents:

.. include:: toc.txt

