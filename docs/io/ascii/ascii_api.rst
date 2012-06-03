.. include:: references.txt

astropy.io.ascii API
===========================

.. automodule:: astropy.io.ascii

Functions
--------------

.. autofunction:: read

.. autofunction:: get_reader

.. autofunction:: write

.. autofunction:: get_writer

.. autofunction:: convert_list

.. autofunction:: convert_numpy

.. autofunction:: set_guess

Core Classes
--------------
.. autoclass:: BaseReader
   :show-inheritance:
   :members:
   :inherited-members:
   :undoc-members:

.. autoclass:: BaseData
   :show-inheritance:
   :members:
   :inherited-members:
   :undoc-members:

.. autoclass:: BaseHeader
   :show-inheritance:
   :members:
   :inherited-members:
   :undoc-members:

.. autoclass:: BaseInputter
   :show-inheritance:
   :members:
   :inherited-members:
   :undoc-members:

.. autoclass:: BaseOutputter
   :show-inheritance:
   :members:
   :inherited-members:
   :undoc-members:

.. autoclass:: BaseSplitter
   :show-inheritance:
   :members:
   :inherited-members:
   :undoc-members:

.. autoclass:: Column
   :show-inheritance:
   :members:
   :inherited-members:
   :undoc-members:

.. autoclass:: DefaultSplitter
   :show-inheritance:
   :members:
   :inherited-members:
   :undoc-members:

.. autoclass:: DictLikeNumpy
   :show-inheritance:
   :members:
   :undoc-members:

.. autoclass:: InconsistentTableError
   :show-inheritance:

.. autoclass:: NumpyOutputter
   :show-inheritance:
   :members:
   :inherited-members:
   :undoc-members:


.. _extension_reader_classes:

Extension Reader Classes
-------------------------

The following classes extend the base Reader functionality to handle different
table formats.  Some, such as the :class:`Basic` Reader class are fairly
general and include a number of configurable attributes.  Others such as
:class:`Cds` or :class:`Daophot` are specialized to read certain well-defined
but idiosyncratic formats.

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

.. autoclass:: AASTex
   :show-inheritance:
   :members:
   :undoc-members:

.. autoclass:: Basic
   :show-inheritance:
   :members:
   :undoc-members:

.. autoclass:: Cds
   :show-inheritance:
   :members:
   :undoc-members:

.. autoclass:: CommentedHeader
   :show-inheritance:
   :members:
   :undoc-members:

.. autoclass:: Daophot
   :show-inheritance:
   :members:
   :undoc-members:

.. autoclass:: FixedWidth
   :show-inheritance:
   :members:
   :undoc-members:

.. autoclass:: FixedWidthNoHeader
   :show-inheritance:
   :members:
   :undoc-members:

.. autoclass:: FixedWidthTwoLine
   :show-inheritance:
   :members:
   :undoc-members:

.. autoclass:: Ipac
   :show-inheritance:
   :members:
   :undoc-members:

.. autoclass:: Latex
   :show-inheritance:
   :members:
   :undoc-members:

.. autoclass:: Memory
   :show-inheritance:
   :members:
   :undoc-members:

.. autoclass:: NoHeader
   :show-inheritance:
   :members:
   :undoc-members:

.. autoclass:: Rdb
   :show-inheritance:
   :members:
   :undoc-members:

.. autoclass:: Tab
   :show-inheritance:
   :members:
   :undoc-members:

Other extension classes
-----------------------
These classes provide support for extension readers.

.. autoclass:: astropy.io.ascii.cds.CdsData
   :show-inheritance:
   :members:
   :undoc-members:

.. autoclass:: astropy.io.ascii.cds.CdsHeader
   :show-inheritance:
   :members:
   :undoc-members:

.. autoclass:: astropy.io.ascii.basic.CommentedHeaderHeader
   :show-inheritance:
   :members:
   :undoc-members:

.. autoclass:: ContinuationLinesInputter
   :show-inheritance:
   :members:
   :undoc-members:

.. autoclass:: astropy.io.ascii.daophot.DaophotHeader
   :show-inheritance:
   :members:
   :undoc-members:

.. autoclass:: FixedWidthSplitter
   :show-inheritance:
   :members:
   :undoc-members:

.. autoclass:: FixedWidthHeader
   :show-inheritance:
   :members:
   :undoc-members:

.. autoclass:: FixedWidthData
   :show-inheritance:
   :members:
   :undoc-members:

.. autoclass:: astropy.io.ascii.ipac.IpacData
   :show-inheritance:
   :members:
   :undoc-members:

.. autoclass:: astropy.io.ascii.ipac.IpacHeader
   :show-inheritance:
   :members:
   :undoc-members:

.. autoclass:: astropy.io.ascii.latex.LatexHeader
   :show-inheritance:
   :members:
   :undoc-members:

.. autoclass:: astropy.io.ascii.latex.LatexData
   :show-inheritance:
   :members:
   :undoc-members:

.. autoclass:: astropy.io.ascii.latex.LatexSplitter
   :show-inheritance:
   :members:
   :undoc-members:

.. autoclass:: astropy.io.ascii.latex.AASTexHeader
   :show-inheritance:
   :members:
   :undoc-members:

.. autoclass:: astropy.io.ascii.latex.AASTexData
   :show-inheritance:
   :members:
   :undoc-members:

.. autoclass:: astropy.io.ascii.latex.AASTexHeaderSplitter
   :show-inheritance:
   :members:
   :undoc-members:


.. toctree::
   :maxdepth: 1

