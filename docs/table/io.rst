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
custom reader/writer`_). The :meth:`~astropy.table.table.Table.read` method should be used as::

    t = Table.read(filename, format='format')

where ``'format'`` is the format of the file to read in, e.g.::

    t = Table.read('photometry.dat', format='daophot')

For certain file formats, the format can be automatically detected, for
example from the filename extension::

    t = Table.read('table.tex')

Similarly, for writing, the format can be explicitly specified::

    t.write(filename, format='format')

and as for the :meth:`~astropy.table.table.Table.read` method, the format can
be automatically identified from the extension.

Any additional arguments specified will be passed to the format-specific
read/write functions (see e.g. see `Built-in readers/writers`_) - for
example, in the following case::

    t = Table.read('photometry.dat', format='ascii', data_start=2, delimiter='|')

the ``data_start`` and ``delimiter`` arguments are passed to :func:`astropy.io.ascii.ui.read`.

Built-in readers/writers
^^^^^^^^^^^^^^^^^^^^^^^^

ASCII formats
"""""""""""""

The :meth:`~astropy.table.table.Table.read` and
:meth:`~astropy.table.table.Table.write` methods can be used to read and write formats
supported by `astropy.io.ascii`:

    * ``format='ascii'`` can be used to interface to the bare
      :func:`~astropy.io.ascii.ui.read` and :func:`~astropy.io.ascii.ui.write`
      functions from `astropy.io.ascii`, e.g.::

         t = Table.read('table.tex', format='ascii')

         from astropy.io.ascii import Rdb
         t.write('table.rdb', format='ascii', Writer=Rdb)

    * ``format='ipac'`` can be used to read IPAC tables

    * ``format='cds'`` can be used to read CDS/Machine readable tables

    * ``format='daophot'`` can be used to read DAOPhot tables

    * ``format='latex'`` can be used to read and write Latex tables

    * ``format='rdb'`` can be used to read and write RDB tables

The following file extensions are automatically recognized and do not require the ``format=`` argument to be specified:

    * ``.tex`` is recognized as ``format='latex'``
    * ``.rdb`` is recognized as ``format='rdb'``

HDF5
""""

Reading/writing from/to HDF5 files is supported with ``format='hdf5'``. The
``.hdf5`` file extension is automatically recognized when writing files, and
HDF5 files are automatically identified (even with a different extension) when
reading in (using the first few bytes of the file to identify the format).

For information on available arguments, see
:func:`~astropy.io.misc.hdf5.read_table_hdf5` and
:func:`~astropy.io.misc.hdf5.write_table_hdf5`.

Other
"""""

In future, FITS and VO tables will also be supported.

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
