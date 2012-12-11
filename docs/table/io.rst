.. _table_io:

Reading and writing Table objects
---------------------------------

Introduction
^^^^^^^^^^^^

The :class:`~astropy.table.table.Table` class includes two methods,
:meth:`~astropy.table.table.Table.read` and
:meth:`~astropy.table.table.Table.write`, that make it possible to read from
and write to files. A number of formats are automatically supported (see
`Built-in readers/writers`_) and new file formats and extensions can be
registered with the :class:`~astropy.table.table.Table` class (see `Creating a
custom reader/writer`_). After importing the :class:`~astropy.table.table.Table` class::

    >>> from astropy.table import Table

the :meth:`~astropy.table.table.Table.read` method should be used as::

    >>> t = Table.read(filename, format='format')

where ``'format'`` is the format of the file to read in, e.g.::

    >>> t = Table.read('photometry.dat', format='daophot')

For certain file formats, the format can be automatically detected, for
example from the filename extension::

    >>> t = Table.read('table.tex')

Similarly, for writing, the format can be explicitly specified::

    >>> t.write(filename, format='format')

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
supported by `astropy.io.ascii`:

IPAC
++++

`IPAC tables <http://irsa.ipac.caltech.edu/applications/DDGEN/Doc/ipac_tbl.html>`_
can be read with ``format='ipac'``::

  >>> t = Table.read('2mass.tbl', format='ipac')

Note that there are different conventions for characters occuring below the
position of the ``|`` symbol in IPAC tables. By default, any character
below a ``|`` will be ignored (since this is the current standard),
but if you need to read files that assume characters below the ``|``
symbols belong to the column before or after the ``|``, you can specify
``definition='left'`` or ``definition='right'`` respectively when reading
the table (the default is ``definition='ignore'``). The following examples demonstrate the different conventions:

* ``definition='ignore'``::

    |   ra  |  dec  |
    | float | float |
      1.2345  6.7890

* ``definition='left'``::

    |   ra  |  dec  |
    | float | float |
       1.2345  6.7890

* ``definition='right'``::

    |   ra  |  dec  |
    | float | float |
    1.2345  6.7890


Advanced information is available in the :class:`~astropy.io.ascii.ipac.Ipac`
class (any arguments apart from the filename and ``format`` are passed to
this class when ``format='ipac'``).

CDS/Machine Readable
++++++++++++++++++++

`CDS/Machine readable tables <http://vizier.u-strasbg.fr/doc/catstd.htx>`_ can be read with ``format='cds'``::

    >>> t = Table.read('aj285677t3.txt', format='cds')

If the table definition is given in a separate ``ReadMe`` file, this can be
specified with::

    >>> t = Table.read('aj285677t3.txt', format='cds', readme="ReadMe")

Advanced information is available in the :class:`~astropy.io.ascii.cds.Cds`
class (any arguments apart from the filename and ``format`` are passed to
this class when ``format='cds'``).

DAOPhot
+++++++

`DAOPhot <http://stsdas.stsci.edu/cgi-bin/gethelp.cgi?daophot.hlp>`_ tables
can be read with ``format='daophot'``::

  >>> t = Table.read('photometry.dat', format='daophot')

Advanced information is available in the
:class:`~astropy.io.ascii.daophot.Daophot` class (any arguments apart from
the filename and ``format`` are passed to this class when
``format='daophot'``).

LaTeX
+++++

`LaTeX <http://www.latex-project.org/>`_ tables can be read and written with
``format='latex'``. Provided the ``.tex``` extension is used, the format does
not need to be explicitly specified::

      >>> t = Table.read('paper_table.tex')
      >>> t.write('new_paper_table.tex')

If a different extension is used, the format should be specified::

      >>> t.write('new_paper_table.inc', format='latex')

Advanced information is available in the
:class:`~astropy.io.ascii.latex.Latex` class (any arguments apart from the
filename and ``format`` are passed to this class  when ``format='latex'``).

RDB
+++

`RDB <http://hea-www.harvard.edu/MST/simul/software/docs/rdb.html>`_ tables
can be read and written with ``format='rdb'`` Provided the ``.rdb`` extension
is used, the format does not need to be explicitly specified::

      >>> t = Table.read('discovery_data.rdb')
      >>> t.write('updated_data.rdb')

If a different extension is used, the format should be specified::

      >>> t.write('updated_data.txt', format='rdb')

Advanced information is available in the :class:`~astropy.io.ascii.basic.Rdb`
class (any arguments apart from the filename and ``format`` are passed to
this class when ``format='rdb'``).

Arbitrary ASCII formats
+++++++++++++++++++++++

``format='ascii'`` can be used to interface to the bare
:func:`~astropy.io.ascii.ui.read` and :func:`~astropy.io.ascii.ui.write`
functions from `astropy.io.ascii`, e.g.::

       >>> t = Table.read('table.tex', format='ascii')

All additional arguments are passed to the `astropy.io.ascii`
:func:`~astropy.io.ascii.ui.read` and
:func:`~astropy.io.ascii.ui.write`. For example, in the following case::

       >>> t = Table.read('photometry.dat', format='ascii', data_start=2, delimiter='|')

the ``data_start`` and ``delimiter`` arguments are passed to the
:func:`~astropy.io.ascii.ui.read` function from `astropy.io.ascii` (and
similarly for writing).

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

Creating a custom reader/writer
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following example demonstrates how to create a reader for the
Table class. First, we can create a highly simplistic FITS reader
which just reads the data as a structured array::

    from astropy.table import Table

    def fits_reader(filename, hdu=1):
        from astropy.io import fits
        data = fits.open(filename)[hdu].data
        return Table(data)

and then register it with astropy.table::

    from astropy.table import io_registry
    io_registry.register_reader('fits', fits_reader)

Reader functions can take any arguments except ``format`` (since this
is reserved for the ``Table.read`` method) and should return a
``Table`` object.

We can then read in a FITS table with::

    t = Table.read('catalog.fits', format='fits')

In practice, it would be nice to have the ``read`` method automatically
identify that this file was a FITS file, so we can construct a function that
can recognize FITS files, which we refer to here as an *identifier*
function. An identifier function should take three arguments: the first
should be a string which indicates whether the identifier is being called
from ``read`` or ``write``, and the second and third are the positional and
keyword arguments passed to ``Table.read`` respectively (and are therefore a
list and a dictionary). We can write a simplistic function that only looks
at filenames (but in practice, this function could even look at the first
few bytes of the file for example). The only requirement is that it return a
boolean indicating whether the input matches that expected for the format::

    def fits_identify(origin, args, kwargs):
        return isinstance(args[0], basestring) and \
               args[0].lower().split('.')[-1] in ['fits', 'fit']

We then register this identifier function with ``astropy.table``::

    io_registry.register_identifier('fits', fits_identify)

And we can then do::

    t = Table.read('catalog.fits')

If multiple formats match the current input, then an exception is
raised, and similarly if no format matches the current input. In that
case, the format should be explicitly given with the ``format=``
keyword argument.

Similarly, it is possible to create custom writers. To go with our simplistic FITS reader above, we can write a simplistic FITS writer::

   def fits_writer(table, filename, clobber=False):
       import numpy as np
       from astropy.io import fits
       fits.writeto(filename, np.array(table), clobber=clobber)

We then register the writer::

   io_registry.register_writer('fits', fits_writer)

And we can then write the file out to a FITS file::

   t.write('catalog_new.fits', format='fits')

If we have registered the identifier as above, we can simply do::

   t.write('catalog_new.fits')
