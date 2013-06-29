.. _read_write_tables:

Reading and writing Table objects
===================================

Astropy provides a unified interface for reading and writing data
in different formats.  For many common cases this will 
simplify the process of file I/O and reduce the need to master
the separate details of all the I/O packages within Astropy.  For details and 
examples of using this interface see the :ref:`table_io` 
section.

The :class:`~astropy.table.table.Table` class includes two methods,
:meth:`~astropy.table.table.Table.read` and
:meth:`~astropy.table.table.Table.write`, that make it possible to read from
and write to files. A number of formats are automatically supported (see
:ref:`built_in_readers_writers`) and new file formats and extensions can be
registered with the :class:`~astropy.table.table.Table` class (see
:ref:`io_registry`).
::

    >>> from astropy.table import Table

the :meth:`~astropy.table.table.Table.read` method should be used as::

    >>> t = Table.read(filename, format='format')

where ``'format'`` is the format of the file to read in, e.g.::

    >>> t = Table.read('photometry.dat', format='ascii.daophot')

For certain file formats, the format can be automatically detected, for
example from the filename extension::

    >>> t = Table.read('table.tex')

Similarly, for writing, the format can be explicitly specified, e.g.::

    >>> t.write(filename, format='fits')

but as for the :meth:`~astropy.table.table.Table.read` method, the format may
be automatically identified in some cases.

Any additional arguments specified will depend on the format (see e.g. see
`Built-in readers/writers`_)

Built-in readers/writers
^^^^^^^^^^^^^^^^^^^^^^^^

ASCII formats
"""""""""""""

The :meth:`~astropy.table.table.Table.read` and
:meth:`~astropy.table.table.Table.write` methods can be used to read and write formats
supported by `astropy.io.ascii`.

Use ``format='ascii'`` in order to interface to the generic
:func:`~astropy.io.ascii.ui.read` and :func:`~astropy.io.ascii.ui.write`
functions from `astropy.io.ascii`.  When reading a table this means
that all supported ASCII table formats will be tried in order to successfully
parse the input.  For example::

  >>> t = Table.read('astropy/io/ascii/tests/t/latex1.tex', format='ascii')
  >>> print t
  cola colb colc
  ---- ---- ----
     a    1    2
     b    3    4

When writing a table with ``format='ascii'`` the output is a basic
character-delimited file with a single header line containing the
column names.

All additional arguments are passed to the `astropy.io.ascii`
:func:`~astropy.io.ascii.ui.read` and
:func:`~astropy.io.ascii.ui.write` functions. For example, to change
column  delimiter and the output format for the ``colc`` column use::

  >>> t.write(sys.stdout, format='ascii', delimiter='|', formats={'colc': '%0.2f'})
  cola|colb|colc
  a|1|2.00
  b|3|4.00

A full list of the supported ``format`` values and corresponding format types
for ASCII tables is given below.  The ``Suffix`` column indicates the filename
suffix where the format will be auto-detected, while the ``Write`` column
indicates which support write functionality.

=========================== ====== ===== ============================================================================================
           Format           Suffix Write                                          Description                                          
=========================== ====== ===== ============================================================================================
ascii                                Yes ASCII table in any supported format (uses guessing)                                           
ascii.aastex                         Yes :class:`~astropy.io.ascii.latex.AASTex`: AASTeX deluxetable used for AAS journals             
ascii.basic                          Yes :class:`~astropy.io.ascii.basic.Basic`: Basic table with custom delimiters                    
ascii.cds                                :class:`~astropy.io.ascii.cds.Cds`: CDS format table                                          
ascii.commented_header               Yes :class:`~astropy.io.ascii.basic.CommentedHeader`: Column names in a commented line            
ascii.daophot                            :class:`~astropy.io.ascii.daophot.Daophot`: IRAF DAOphot format table                         
ascii.fixed_width                    Yes :class:`~astropy.io.ascii.fixedwidth.FixedWidth`: Fixed width                                 
ascii.fixed_width_no_header          Yes :class:`~astropy.io.ascii.fixedwidth.FixedWidthNoHeader`: Fixed width with no header          
ascii.fixed_width_two_line           Yes :class:`~astropy.io.ascii.fixedwidth.FixedWidthTwoLine`: Fixed width with second header line
ascii.ipac                           Yes :class:`~astropy.io.ascii.ipac.Ipac`: IPAC format table                                       
ascii.latex                   .tex   Yes :class:`~astropy.io.ascii.latex.Latex`: LaTeX table                                           
ascii.no_header                      Yes :class:`~astropy.io.ascii.basic.NoHeader`: Basic table with no headers                        
ascii.rdb                     .rdb   Yes :class:`~astropy.io.ascii.basic.Rdb`: Tab-separated with a type definition header line        
ascii.sextractor                         :class:`~astropy.io.ascii.sextractor.SExtractor`: SExtractor format table                     
ascii.tab                            Yes :class:`~astropy.io.ascii.basic.Tab`: Basic table with tab-separated values                   
=========================== ====== ===== ============================================================================================

HDF5
""""

Reading/writing from/to `HDF5 <http://www.hdfgroup.org/HDF5/>`_ files is
supported with ``format='hdf5'`` (this requires `h5py
<http://code.google.com/p/h5py/>`_ to be installed). However, the ``.hdf5``
file extension is automatically recognized when writing files, and HDF5 files
are automatically identified (even with a different extension) when reading
in (using the first few bytes of the file to identify the format), so in most
cases you will not need to explicitly specify ``format='hdf5'``.

Since HDF5 files can contain multiple tables, the full path to the table
should be specified via the ``path=`` argument when reading and writing.
For example, to read a table called ``data`` from an HDF5 file named
``observations.hdf5``, you can do::

    >>> t = Table.read('observations.hdf5', path='data')

To read a table nested in a group in the HDF5 file, you can do::

    >>> t = Table.read('observations.hdf5', path='group/data')

To write a table to a new file, the path should also be specified::

    >>> t.write('new_file.hdf5', path='updated_data')

It is also possible to write a table to an existing file using ``append=True``::

    >>> t.write('observations.hdf5', path='updated_data', append=True)

Finally, when writing to HDF5 files, the ``compression=`` argument can be
used to ensure that the data is compressed on disk::

    >>> t.write('new_file.hdf5', path='updated_data', compression=True)

As with other formats, the ``overwrite=True`` argument is supported for
overwriting existing files.

VO Tables
"""""""""

Reading/writing from/to `VO table <http://www.ivoa.net/Documents/VOTable/>`_
files is supported with ``format='votable'``. In most cases, existing VO
tables should be automatically identified as such based on the header of the
file, but if not, or if writing to disk, then the format should be explicitly
specified.

If a VO table file only contains a single table, then it can be read in with::

    >>> t = Table.read('aj285677t3_votable.xml')

If more that one table are present in the file, an error will be raised,
unless the table ID is specified via the ``table_id=`` argument::

    >>> t = Table.read('catalog.xml')
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
      File "/Volumes/Raptor/Library/Python/2.7/lib/python/site-packages/astropy/table/table.py", line 1559, in read
        table = reader(*args, **kwargs)
      File "/Volumes/Raptor/Library/Python/2.7/lib/python/site-packages/astropy/io/votable/connect.py", line 44, in read_table_votable
        raise ValueError("Multiple tables found: table id should be set via the id= argument. The available tables are " + ', '.join(tables.keys()))
    ValueError: Multiple tables found: table id should be set via the table_id= argument. The available tables are twomass, spitzer

    >>> t = Table.read('catalog.xml', table_id='twomass')

To write to a new file, the ID of the table should also be specified (unless
``t.meta['ID']`` is defined)::

    >>> t.write('new_catalog.xml', table_id='updated_table', format='votable')

When writing, the ``compression=True`` argument can be used to force
compression of the data on disk, and the ``overwrite=True`` argument can be
used to overwrite an existing file.

Other
"""""

In future, FITS tables will also be supported via the
:class:`~astropy.table.table.Table` class. For now, these can be read and
written directly with `astropy.io.fits`.
