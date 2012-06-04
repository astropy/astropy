.. include:: references.txt

.. _extension_reader_classes:

Extension Reader classes
==========================

The following classes extend the base :class:`~astropy.io.ascii.core.BaseReader` functionality to handle reading and writing
different table formats.  Some, such as the :class:`~astropy.io.ascii.Basic` Reader class
are fairly general and include a number of configurable attributes.  Others
such as :class:`~astropy.io.ascii.Cds` or :class:`~astropy.io.ascii.Daophot` are specialized to read certain
well-defined but idiosyncratic formats.

* :class:`~astropy.io.ascii.AASTex`: AASTeX `deluxetable` used for AAS journals
* :class:`~astropy.io.ascii.Basic`: basic table with customizable delimiters and header configurations
* :class:`~astropy.io.ascii.Cds`: `CDS format table <http://vizier.u-strasbg.fr/doc/catstd.htx>`_ (also Vizier and ApJ machine readable tables)
* :class:`~astropy.io.ascii.CommentedHeader`: column names given in a line that begins with the comment character
* :class:`~astropy.io.ascii.Daophot`: table from the IRAF DAOphot package
* :class:`~astropy.io.ascii.FixedWidth`: table with fixed-width columns (see also :ref:`fixed_width_gallery`)
* :class:`~astropy.io.ascii.FixedWidthNoHeader`: table with fixed-width columns and no header
* :class:`~astropy.io.ascii.FixedWidthTwoLine`: table with fixed-width columns and a two-line header
* :class:`~astropy.io.ascii.Ipac`: `IPAC format table <http://irsa.ipac.caltech.edu/applications/DDGEN/Doc/ipac_tbl.html>`_
* :class:`~astropy.io.ascii.Latex`: LaTeX table with datavalue in the `tabular` environment
* :class:`~astropy.io.ascii.NoHeader`: basic table with no header where columns are auto-named
* :class:`~astropy.io.ascii.Rdb`: tab-separated values with an extra line after the column definition line
* :class:`~astropy.io.ascii.Tab`: tab-separated values

.. autoclass:: astropy.io.ascii.AASTex
   :show-inheritance:
   :members:
   :undoc-members:

.. autoclass:: astropy.io.ascii.Basic
   :show-inheritance:
   :members:
   :undoc-members:

.. autoclass:: astropy.io.ascii.Cds
   :show-inheritance:
   :members:
   :undoc-members:

.. autoclass:: astropy.io.ascii.CommentedHeader
   :show-inheritance:
   :members:
   :undoc-members:

.. autoclass:: astropy.io.ascii.Daophot
   :show-inheritance:
   :members:
   :undoc-members:

.. autoclass:: astropy.io.ascii.FixedWidth
   :show-inheritance:
   :members:
   :undoc-members:

.. autoclass:: astropy.io.ascii.FixedWidthNoHeader
   :show-inheritance:
   :members:
   :undoc-members:

.. autoclass:: astropy.io.ascii.FixedWidthTwoLine
   :show-inheritance:
   :members:
   :undoc-members:

.. autoclass:: astropy.io.ascii.Ipac
   :show-inheritance:
   :members:
   :undoc-members:

.. autoclass:: astropy.io.ascii.Latex
   :show-inheritance:
   :members:
   :undoc-members:

.. autoclass:: astropy.io.ascii.Memory
   :show-inheritance:
   :members:
   :undoc-members:

.. autoclass:: astropy.io.ascii.NoHeader
   :show-inheritance:
   :members:
   :undoc-members:

.. autoclass:: astropy.io.ascii.Rdb
   :show-inheritance:
   :members:
   :undoc-members:

.. autoclass:: astropy.io.ascii.Tab
   :show-inheritance:
   :members:
   :undoc-members:

